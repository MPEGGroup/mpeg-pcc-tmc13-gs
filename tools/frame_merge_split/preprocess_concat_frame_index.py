from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

if __package__:
    from .ply_io import drop_field, prune_by_score, read_ply, write_ply
else:
    sys.path.append(str(Path(__file__).resolve().parent))
    from ply_io import drop_field, prune_by_score, read_ply, write_ply


def glob_source_files(pattern: str, frame_count: int | None = None) -> list[Path]:
    files = sorted(Path("/").glob(pattern[1:]) if pattern.startswith("/") else Path().glob(pattern))
    if frame_count is not None:
        files = files[:frame_count]
    return files


def add_frame_index_ids(data: np.ndarray, frame_ids: np.ndarray, frame_id_scale: int = 1) -> np.ndarray:
    if data.dtype.names is None:
        raise ValueError("Structured PLY data required")
    if "frame_index" in data.dtype.names:
        raise ValueError("Input PLY already has a frame_index property; refusing to overwrite it")
    max_frame_id = int(frame_ids.max()) if len(frame_ids) else 0
    max_carrier = max_frame_id * int(frame_id_scale)
    if max_carrier > np.iinfo(np.uint8).max:
        raise ValueError("Frame ids exceed uint8 frame_index range (max 255)")

    dtype = np.dtype(data.dtype.descr + [("frame_index", "u1")])
    out = np.empty(len(data), dtype=dtype)
    for name in data.dtype.names:
        out[name] = data[name]
    out["frame_index"] = (frame_ids.astype(np.uint16, copy=False) * int(frame_id_scale)).astype(np.uint8)
    return out


def contiguous_groups(frame_count: int, group_size: int) -> list[list[int]]:
    if group_size <= 0:
        return [list(range(frame_count))]
    groups: list[list[int]] = []
    for start in range(0, frame_count, group_size):
        end = min(start + group_size, frame_count)
        groups.append(list(range(start, end)))
    return groups


def write_grouped_contiguous(args: argparse.Namespace, source_files: list[Path]) -> None:
    out_dir = Path(args.out_dir)
    frame_count = len(source_files)
    groups = contiguous_groups(frame_count, int(args.group_frame_count))

    manifest = {
        "strategy": "contiguous",
        "group_frame_count": int(args.group_frame_count),
        "prune": float(args.prune),
        "score_field": str(args.score_field),
        "frame_count": int(frame_count),
        "frame_id_carrier": "uchar frame_index (per-group LOCAL id)",
        "groups": [],
    }

    for gi, members in enumerate(groups):
        g = len(members)
        if args.frame_id_scale > 0:
            scale = int(args.frame_id_scale)
        else:
            scale = max((2 ** 8 - 1) // max(g - 1, 1), 1)

        parts, ids = [], []
        for local_id, fidx in enumerate(members):
            data = read_ply(source_files[fidx], mmap=False).data
            data = prune_by_score(data, float(args.prune), str(args.score_field))
            data = drop_field(data, str(args.score_field))
            parts.append(data)
            ids.append(np.full(len(data), local_id, dtype=np.uint16))

        concat = np.concatenate(parts)
        local_ids = np.concatenate(ids)
        concat = add_frame_index_ids(concat, local_ids, scale)

        name = f"GOF_g{gi:02d}.ply"
        out_path = out_dir / name
        write_ply(out_path, concat)
        (out_dir / f"GOF_g{gi:02d}.frames").write_text(",".join(str(m) for m in members), encoding="utf-8")

        sizes = [int(len(p)) for p in parts]
        manifest["groups"].append(
            {
                "name": name,
                "frames": [int(m) for m in members],
                "frame_sizes": sizes,
                "frame_id_scale": int(scale),
                "total_points": int(len(concat)),
            }
        )
        print(
            f"[concat-frame-index-preprocess] {name}: frames={members} scale={scale} "
            f"({len(concat):,} pts -> {out_path.stat().st_size / 1e9:.2f} GB)"
        )

    with (out_dir / "grouping_manifest.json").open("w", encoding="utf-8") as fp:
        json.dump(manifest, fp, indent=2)

    print(f"[concat-frame-index-preprocess] grouped contiguous done -> {out_dir}")


def run_preprocess(args: argparse.Namespace) -> None:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    source_files = glob_source_files(args.input_glob, args.frame_count)
    if not source_files:
        raise FileNotFoundError(f"No source files matched: {args.input_glob}")
    print(f"[concat-frame-index-preprocess] loading {len(source_files)} frames ...")

    plys = [read_ply(f, mmap=False) for f in source_files]
    frames = [p.data for p in plys]

    frame_count = len(frames)
    frame_sizes = [len(f) for f in frames]
    total_points = sum(frame_sizes)
    print(f"[concat-frame-index-preprocess] total Gaussians: {total_points:,}  ({frame_count} frames)")

    if frame_count > 256:
        raise ValueError("frame_count exceeds uint8 frame_index capacity (max 256 frames)")

    if args.group_frame_count > 0:
        write_grouped_contiguous(args, source_files)
        return

    # Default simulation uses unscaled local ids.
    # Use --frame-id-scale 0 to enable automatic spreading across uint8.
    frame_id_scale = args.frame_id_scale
    if frame_id_scale <= 0:
        frame_id_scale = max((2 ** 8 - 1) // max(frame_count - 1, 1), 1)
    print(f"[concat-frame-index-preprocess] frame_id_scale = {frame_id_scale}")

    concat_data = np.concatenate(frames)
    frame_ids = np.repeat(np.arange(frame_count, dtype=np.uint8), frame_sizes)
    concat_data = add_frame_index_ids(concat_data, frame_ids, frame_id_scale)

    # Point order is irrelevant: every point carries its own frame id in the
    # frame_index attribute, and the GPCC octree re-derives its coding order from
    # the voxel positions, so no global (Morton) sort of the input is needed.
    concat_ply_path = out_dir / "concat.ply"
    write_ply(concat_ply_path, concat_data)
    print(
        f"[concat-frame-index-preprocess] wrote {concat_ply_path}  "
        f"({concat_ply_path.stat().st_size / 1e9:.2f} GB)"
    )

    metadata = {
        "frame_count": frame_count,
        "frame_sizes": frame_sizes,
        "total_points": total_points,
        "source_files": [str(f) for f in source_files],
        "frame_id_carrier": "uchar frame_index",
        "frame_id_scale": int(frame_id_scale),
    }
    with (out_dir / "concat_metadata.json").open("w", encoding="utf-8") as fp:
        json.dump(metadata, fp, indent=2)

    print(f"[concat-frame-index-preprocess] done -> {out_dir}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Contiguous preprocessor for ConcatGPCC. "
            "Either write one concat.ply (default) or contiguous GOF_gXX.ply groups."
        )
    )
    parser.add_argument("--input-glob", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--frame-count", type=int, default=None)
    parser.add_argument(
        "--group-frame-count",
        type=int,
        default=0,
        help="frames per contiguous group; 0 = single concat.ply (no grouping)",
    )
    parser.add_argument(
        "--prune",
        type=float,
        default=0.0,
        help="fraction of least-important Gaussians to prune per frame in grouped mode",
    )
    parser.add_argument(
        "--score-field",
        default="score",
        help="importance-score column used for pruning and dropped in grouped mode",
    )
    parser.add_argument(
        "--frame-id-scale",
        type=int,
        default=1,
        help="uchar multiplier per frame id; default=1 (no scaling), 0 = auto spread across uint8",
    )
    args = parser.parse_args()
    if args.frame_id_scale < 0:
        raise ValueError("--frame-id-scale must be >= 0 (0 = auto)")
    if args.group_frame_count < 0:
        raise ValueError("--group-frame-count must be >= 0")
    if not (0.0 <= args.prune < 1.0):
        raise ValueError("--prune must be in [0, 1)")
    run_preprocess(args)


if __name__ == "__main__":
    main()
