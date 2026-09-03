from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

if __package__:
    from .ply_io import read_ply, write_ply
else:
    sys.path.append(str(Path(__file__).resolve().parent))
    from ply_io import read_ply, write_ply


def drop_field(data: np.ndarray, field: str) -> np.ndarray:
    names = data.dtype.names or ()
    if field not in names:
        return np.asarray(data)
    keep_names = [name for name in names if name != field]
    dtype = np.dtype([(name, data.dtype.fields[name][0]) for name in keep_names])
    out = np.empty(len(data), dtype=dtype)
    for name in keep_names:
        out[name] = data[name]
    return out


def find_frame_id_field(data: np.ndarray) -> str:
    names = data.dtype.names or ()
    if "frame_index" in names:
        return "frame_index"
    if "frameindex" in names:
        return "frameindex"
    # Backward compatibility for older preprocessed streams.
    if "reflectance" in names:
        return "reflectance"
    # TMC13 may write the reflectance attribute using its internal PLY name.
    if "refc" in names:
        return "refc"
    raise ValueError("Decoded PLY has no frame_index/reflectance/refc property; cannot recover frame ids")


def parse_frame_map(frame_map: str) -> list[int]:
    txt = frame_map.strip()
    if not txt:
        raise ValueError("--frame-map cannot be empty")
    return [int(tok) for tok in txt.split(",") if tok.strip()]


def resolve_scale(args: argparse.Namespace, meta: dict, frame_count: int) -> int:
    if args.frame_id_scale is not None:
        return int(args.frame_id_scale)
    if "frame_id_scale" in meta:
        return int(meta["frame_id_scale"])
    return 1


def run_postprocess(args: argparse.Namespace) -> None:
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    use_frame_map = bool(args.frame_map)
    meta: dict = {}
    expected = np.asarray([], dtype=np.int64)

    if use_frame_map:
        abs_ids = parse_frame_map(args.frame_map)
        frame_count = len(abs_ids)
        if frame_count < 1:
            raise ValueError("--frame-map must contain at least one frame id")
    else:
        if not args.preprocessed_dir:
            raise ValueError("--preprocessed-dir is required when --frame-map is not provided")
        pre_dir = Path(args.preprocessed_dir)
        with (pre_dir / "concat_metadata.json").open("r", encoding="utf-8") as fp:
            meta = json.load(fp)
        frame_count = int(args.frame_count or meta["frame_count"])
        expected = np.asarray(meta.get("frame_sizes", []), dtype=np.int64)
        abs_ids = list(range(frame_count))

    frame_id_scale = resolve_scale(args, meta, frame_count)
    if frame_id_scale < 1:
        raise ValueError("frame_id_scale must be >= 1")

    decoded_path = Path(args.decoded_ply)
    print(f"[concat-frame-index-postprocess] loading decoded PLY: {decoded_path}")
    decoded = read_ply(decoded_path, mmap=True)
    frame_id_field = find_frame_id_field(decoded.data)
    print(f"[concat-frame-index-postprocess] using decoded PLY field '{frame_id_field}' as frame id carrier")

    carrier_values = np.asarray(decoded.data[frame_id_field], dtype=np.uint32)
    if len(carrier_values) and int(carrier_values.max()) > np.iinfo(np.uint8).max:
        raise ValueError(
            f"Decoded {frame_id_field} contains values above uint8 range; "
            "expected uchar frame_index in [0, 255]"
        )

    if frame_id_scale > 1:
        frame_ids = np.rint(carrier_values.astype(np.float64) / frame_id_scale).astype(np.uint8)
    else:
        frame_ids = carrier_values.astype(np.uint8)

    min_id = int(frame_ids.min()) if len(frame_ids) else 0
    max_id = int(frame_ids.max()) if len(frame_ids) else 0
    
    print(min_id, max_id)
    
    if max_id >= frame_count:
        raise ValueError(
            f"Decoded local frame id {max_id}, but frame_count={frame_count}; "
            "frame-id scale or frame-map mismatch"
        )

    report = {
        "frames": [],
        "total_decoded": int(len(decoded.data)),
        "frame_id_carrier": f"uchar {frame_id_field}",
        "frame_id_scale": frame_id_scale,
        "frame_index_min": min_id,
        "frame_index_max": max_id,
        "frame_map": abs_ids if use_frame_map else [],
        "per_frame_counts_ok": None,
        "report_tag": args.report_tag,
    }

    counts = np.bincount(frame_ids, minlength=frame_count)
    if (not use_frame_map) and len(expected) == frame_count:
        report["per_frame_counts_ok"] = bool(np.array_equal(counts[:frame_count], expected))

    for local_id, abs_id in enumerate(abs_ids):
        mask = frame_ids == local_id
        frame_data = drop_field(np.asarray(decoded.data[mask]), frame_id_field)

        out_path = out_dir / f"frame_{abs_id:04d}.ply"
        write_ply(out_path, frame_data)
        report["frames"].append(
            {
                "local_frame": int(local_id),
                "frame": int(abs_id),
                "gaussians": int(len(frame_data)),
                "path": str(out_path),
            }
        )
        print(
            f"[concat-frame-index-postprocess] local {local_id:04d} -> frame {abs_id:04d}: "
            f"{len(frame_data):,} Gaussians -> {out_path.name}"
        )

    report_name = "concat_frame_index_postprocess_report.json"
    if args.report_tag:
        report_name = f"concat_frame_index_postprocess_report_{args.report_tag}.json"

    with (out_dir / report_name).open("w", encoding="utf-8") as fp:
        json.dump(report, fp, indent=2)

    if report["per_frame_counts_ok"] is False:
        expected_counts = expected[:frame_count].astype(int).tolist()
        decoded_counts = counts[:frame_count].astype(int).tolist()
        raise ValueError(
            "per-frame counts do not match original frame sizes; "
            f"expected={expected_counts}, decoded={decoded_counts}"
        )

    print(f"[concat-frame-index-postprocess] done -> {out_dir}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Split a GPCC-decoded concatenated PLY into per-frame PLYs using uchar frame ids. "
            "Supports contiguous mode (--preprocessed-dir) and GoF mode (--frame-map)."
        )
    )
    parser.add_argument("--preprocessed-dir", default="",
                        help="Path containing concat_metadata.json (required in non-GoF mode)")
    parser.add_argument("--decoded-ply", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--frame-count", type=int, default=None)
    parser.add_argument(
        "--frame-id-scale",
        type=int,
        default=1,
        help="local frame id scale used in the decoded carrier; default=1 (no scaling)",
    )
    parser.add_argument("--frame-map", default="",
                        help="Comma-separated absolute frame ids for local ids 0..g-1 (GoF mode)")
    parser.add_argument("--report-tag", default="",
                        help="Suffix for report filename (useful when splitting multiple GOFs)")
    args = parser.parse_args()

    if args.frame_id_scale is not None and args.frame_id_scale < 1:
        raise ValueError("--frame-id-scale must be >= 1 when provided")

    run_postprocess(args)


if __name__ == "__main__":
    main()
