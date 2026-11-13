"""Embed COLMAP camera positions into per-frame Gaussian-splat PLYs.

The MPEG GSC metric renders each frame from camera viewpoints stored as
`camera_position` comments in the PLY header.  This is a thin batch wrapper
around the official `cameraPosition` tool (mpeg-gsc-tools / gsTools): for every
`frame_XXXX.ply` it locates the matching COLMAP reconstruction and writes
`frame_XXXX<camera-output-suffix>` for camera-only mode, or
`frame_XXXX<camera-score-output-suffix>` for camera+score mode.

This is a data-preparation utility only.  The ConcatGPCC pipeline assumes the
PLYs are already embedded (`*.camerapos.ply`) and does not call it.

Expected dataset layout (per sequence):
    <ply-dir>/frame_XXXX.ply
    <colmap-root>/frame_XXXX/sparse/{cameras,images}.{bin,txt}
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path


def embed_one(camerapos_bin: Path, raw_ply: Path, cameras: Path,
              images: Path, out_ply: Path) -> None:
    cmd = [
        str(camerapos_bin),
        "camera",
        f"--input={raw_ply}",
        f"--camera={cameras}",
        f"--image={images}",
        f"--output={out_ply}",
        "-v", "1",
    ]
    subprocess.run(cmd, check=True)


def embed_score(camerapos_bin: Path, raw_ply: Path, score_file: Path,
                out_ply: Path, score_name: str,
                percentile: float, log_scale: bool) -> None:
    cmd = [
        str(camerapos_bin),
        "attribute",
        f"--input={raw_ply}",
        f"--name={score_name}",
        f"--values={score_file}",
        f"--output={out_ply}",
        "--transform",
        f"--percentile={percentile}",
        "-v", "1",
    ]
    if log_scale:
        cmd.append("--log-scale")
    subprocess.run(cmd, check=True)


def run(args: argparse.Namespace) -> int:
    ply_dir = Path(args.ply_dir)
    out_dir = Path(args.output_dir) if args.output_dir else ply_dir
    colmap_root = Path(args.colmap_root) if args.colmap_root else ply_dir.parent / "colmap_data"
    camerapos_bin = Path(args.camerapos_bin)

    if not camerapos_bin.is_file():
        sys.exit(f"ERROR: cameraPosition binary not found: {camerapos_bin}")
    out_dir.mkdir(parents=True, exist_ok=True)

    protected_suffixes = {".with_score.tmp.ply"}
    if args.camera_output_suffix:
        protected_suffixes.add(args.camera_output_suffix)
    if args.camera_score_output_suffix:
        protected_suffixes.add(args.camera_score_output_suffix)
    if args.scores_output_suffix:
        protected_suffixes.add(args.scores_output_suffix)

    raw_files = sorted(
        p for p in ply_dir.glob("frame_*.ply")
        if not any(p.name.endswith(sfx) for sfx in protected_suffixes)
    )
    if args.frame_count is not None:
        raw_files = raw_files[:args.frame_count]
    if not raw_files:
        sys.exit(f"ERROR: no raw frame_*.ply found in {ply_dir}")

    embedded = skipped = missing = 0
    missing_scores = 0
    for raw in raw_files:
        stem = raw.stem  # frame_XXXX
        if args.scores_only:
            out_ply = out_dir / f"{stem}{args.scores_output_suffix}"
            if out_ply.exists() and not args.overwrite:
                skipped += 1
                continue
            score_file = Path(args.scores_root) / f"{stem}{args.scores_suffix}"
            if not score_file.is_file():
                print(f"[embed-camera] WARN missing score file for {stem}: {score_file}", file=sys.stderr)
                missing_scores += 1
                continue
            print(f"[embed-camera] {stem} -> {out_ply.name}")
            embed_score(
                camerapos_bin,
                raw,
                score_file,
                out_ply,
                args.score_name,
                args.score_percentile,
                args.score_log_scale,
            )
            embedded += 1
            continue

        out_suffix = (
            args.camera_score_output_suffix
            if args.embed_scores
            else args.camera_output_suffix
        )
        out_ply = out_dir / f"{stem}{out_suffix}"
        if out_ply.exists() and not args.overwrite:
            skipped += 1
            continue
        sparse = colmap_root / stem / "sparse"
        binary_model = sparse / "cameras.bin", sparse / "images.bin"
        text_model = sparse / "cameras.txt", sparse / "images.txt"
        if all(path.is_file() for path in binary_model):
            cameras, images = binary_model
        elif all(path.is_file() for path in text_model):
            cameras, images = text_model
        else:
            print(
                f"[embed-camera] WARN missing complete COLMAP .bin or .txt model "
                f"for {stem}: {sparse}",
                file=sys.stderr,
            )
            missing += 1
            continue
        input_for_camera = raw
        scored_tmp: Path | None = None
        if args.embed_scores:
            score_file = Path(args.scores_root) / f"{stem}{args.scores_suffix}"
            if not score_file.is_file():
                print(f"[embed-camera] WARN missing score file for {stem}: {score_file}", file=sys.stderr)
                missing_scores += 1
                continue
            with tempfile.NamedTemporaryFile(
                prefix=f"{stem}.",
                suffix=".with_score.tmp.ply",
                dir=str(out_dir),
                delete=False,
            ) as tmp:
                scored_tmp = Path(tmp.name)
            try:
                embed_score(
                    camerapos_bin,
                    raw,
                    score_file,
                    scored_tmp,
                    args.score_name,
                    args.score_percentile,
                    args.score_log_scale,
                )
                input_for_camera = scored_tmp
            except Exception:
                if scored_tmp.exists():
                    scored_tmp.unlink()
                raise
        print(f"[embed-camera] {stem} -> {out_ply.name}")
        try:
            embed_one(camerapos_bin, input_for_camera, cameras, images, out_ply)
        finally:
            if scored_tmp and scored_tmp.exists():
                scored_tmp.unlink()
        embedded += 1

    print(
        f"[embed-camera] done: embedded={embedded} skipped={skipped} "
        f"missing_colmap={missing} missing_scores={missing_scores} -> {out_dir}"
    )
    return 1 if (missing or missing_scores) else 0


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Batch-embed COLMAP camera positions into frame_*.ply (writes frame_*.camerapos.ply)."
    )
    parser.add_argument("--ply-dir", required=True, help="Dir containing frame_XXXX.ply")
    parser.add_argument("--camerapos-bin", required=True,
                        help="Path to the gsTools cameraPosition binary")
    parser.add_argument("--colmap-root", default="",
                        help=(
                            "Dir with <stem>/sparse/{cameras,images}.{bin,txt} "
                            "(default: <ply-dir>/../colmap_data)"
                        ))
    parser.add_argument("--output-dir", default="",
                        help="Where to write *.camerapos.ply (default: --ply-dir)")
    parser.add_argument("--camera-output-suffix", default=".camerapos.ply",
                        help="Output suffix for camera-only mode (default: .camerapos.ply)")
    parser.add_argument("--camera-score-output-suffix", default=".camerapos.ply",
                        help="Output suffix for camera+score mode (default: .camerapos.ply)")
    parser.add_argument("--frame-count", type=int, default=None,
                        help="Embed only the first N frames (default: all)")
    parser.add_argument("--overwrite", action="store_true",
                        help="Re-embed even if the .camerapos.ply already exists")
    parser.add_argument("--embed-scores", action="store_true",
                        help="Embed a score attribute before camera embedding")
    parser.add_argument("--scores-only", action="store_true",
                        help="Embed scores only (writes <frame><scores-output-suffix>)")
    parser.add_argument("--scores-root", default="",
                        help=(
                            "Dir containing score files named <frame><scores-suffix> "
                            "(default: <ply-dir>/../scores)"
                        ))
    parser.add_argument("--scores-suffix", default=".impScores.txt",
                        help="Score file suffix appended to frame stem (default: .impScores.txt)")
    parser.add_argument("--score-name", default="score",
                        help="Attribute name to write for scores (default: score)")
    parser.add_argument("--scores-output-suffix", default=".with_score.ply",
                        help="Output suffix when --scores-only is set (default: .with_score.ply)")
    parser.add_argument("--score-percentile", type=float, default=1.0,
                        help="Percentile for score transform (default: 1.0)")
    parser.add_argument("--score-log-scale", action="store_true", default=True,
                        help="Use --log-scale in score transform (default: enabled)")
    parser.add_argument("--no-score-log-scale", dest="score_log_scale", action="store_false",
                        help="Disable --log-scale in score transform")
    args = parser.parse_args()
    if args.scores_only:
        args.embed_scores = True
    if not args.scores_root:
        args.scores_root = str(Path(args.ply_dir).parent / "scores")
    sys.exit(run(args))


if __name__ == "__main__":
    main()
    
    
'''
./batch_embed_data.sh \
    --dataset-path /media/gs_storage/DATASET/MPEG/CTC_156 \
    --colmap-subfolder colmap_data_eval \
    --config ./testset_inputs/dataset.config.ctc.sh \
    --gstools-bin /mnt/data01/mpeg_integration/mpeg-gsc-tools/gsTools/build/msvc/Release/bin/gstools \
    --embed-scores \
    --scores-subfolder scores \
    --scores-suffix .impScores.txt \
    --scores-output-suffix .jitter_camerapos.ply \
    --score-name score \
    --score-percentile 1.0 \
    --overwrite
'''