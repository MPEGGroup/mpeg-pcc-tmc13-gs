#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

usage() {
    cat <<'EOF'
Usage:
    batch_importance_scores.sh \
        --dataset-path <path> \
        --config <dataset.config.*.sh> \
        [--bin <path-to-importance-scores>] \
        [--plys-subfolder <name>] \
        [--scores-subfolder <name>] \
        [--ply-pattern <printf-pattern>] \
        [--scores-suffix <suffix>] \
        [--start-frame <N>] \
        [--frame-count <N>] \
        [--viewpoint <file>] \
        [--downscale-factor <float>] \
        [--downscale-min-width <N>] \
        [--threads <N>] \
        [--overwrite] \
        [--verbose] \
        [--dry-run]

Description:
    Reads the sequence list from a dataset config file (ADD_SEQ <sequence> <width>
    <height>) and runs importance-scores on every PLY frame of every sequence.

Expected per-sequence layout:
    <dataset-path>/<sequence>/<plys-subfolder>/frame_XXXX.camerapos.ply

Output (one text file per frame):
    <dataset-path>/<sequence>/<scores-subfolder>/frame_XXXX.impScores.txt

Note:
    The source PLY files must carry the "comment camera_position" header entries
    (i.e. the *.camerapos.ply files), unless --viewpoint is provided.

Example:
    ./batch_importance_scores.sh \
        --dataset-path /media/gs_storage/DATASET/MPEG/CTC_156 \
        --config ./testset_inputs/dataset.config.ctc.sh
EOF
}

DATASET_PATH=""
CONFIG_FILE=""
BIN="$REPO_ROOT/build/Release/bin/importance-scores"
PLYS_SUBFOLDER="plys"
SCORES_SUBFOLDER="scores"
PLY_PATTERN="frame_%04d.camerapos.ply"
SCORES_SUFFIX=".impScores.txt"
START_FRAME=""
FRAME_COUNT=""
VIEWPOINT=""
DOWNSCALE_FACTOR=""
DOWNSCALE_MIN_WIDTH=""
THREADS=""
OVERWRITE=0
VERBOSE=0
DRY_RUN=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dataset-path)        DATASET_PATH="$2"; shift 2 ;;
        --config)              CONFIG_FILE="$2"; shift 2 ;;
        --bin)                 BIN="$2"; shift 2 ;;
        --plys-subfolder)      PLYS_SUBFOLDER="$2"; shift 2 ;;
        --scores-subfolder)    SCORES_SUBFOLDER="$2"; shift 2 ;;
        --ply-pattern)         PLY_PATTERN="$2"; shift 2 ;;
        --scores-suffix)       SCORES_SUFFIX="$2"; shift 2 ;;
        --start-frame)         START_FRAME="$2"; shift 2 ;;
        --frame-count)         FRAME_COUNT="$2"; shift 2 ;;
        --viewpoint)           VIEWPOINT="$2"; shift 2 ;;
        --downscale-factor)    DOWNSCALE_FACTOR="$2"; shift 2 ;;
        --downscale-min-width) DOWNSCALE_MIN_WIDTH="$2"; shift 2 ;;
        --threads)             THREADS="$2"; shift 2 ;;
        --overwrite)           OVERWRITE=1; shift ;;
        --verbose)             VERBOSE=1; shift ;;
        --dry-run)             DRY_RUN=1; shift ;;
        -h|--help)             usage; exit 0 ;;
        *)
            echo "ERROR: Unknown argument: $1" >&2
            usage
            exit 1
            ;;
    esac
done

if [[ -z "$DATASET_PATH" || -z "$CONFIG_FILE" ]]; then
    echo "ERROR: --dataset-path and --config are required." >&2
    usage
    exit 1
fi

if [[ ! -d "$DATASET_PATH" ]]; then
    echo "ERROR: dataset path does not exist: $DATASET_PATH" >&2
    exit 1
fi

if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "ERROR: config file does not exist: $CONFIG_FILE" >&2
    exit 1
fi

if [[ ! -x "$BIN" ]]; then
    echo "ERROR: importance-scores binary not found or not executable: $BIN" >&2
    exit 1
fi

if [[ -n "$VIEWPOINT" && ! -f "$VIEWPOINT" ]]; then
    echo "ERROR: viewpoint file does not exist: $VIEWPOINT" >&2
    exit 1
fi

# Glob and index-extraction patterns derived from the printf-style ply pattern.
PLY_GLOB="${PLY_PATTERN//%04d/[0-9][0-9][0-9][0-9]}"
PLY_PREFIX="${PLY_PATTERN%%%04d*}"
PLY_SUFFIX="${PLY_PATTERN#*%04d}"

if [[ "$PLY_GLOB" == "$PLY_PATTERN" ]]; then
    echo "ERROR: --ply-pattern must contain %04d: $PLY_PATTERN" >&2
    exit 1
fi

SEQUENCES=()
ADD_SEQ() { SEQUENCES+=("$1"); }

# shellcheck source=/dev/null
source "$CONFIG_FILE"

if (( ${#SEQUENCES[@]} == 0 )); then
    echo "ERROR: no sequences were registered by config: $CONFIG_FILE" >&2
    exit 1
fi

echo "[imp-scores] dataset:   $DATASET_PATH"
echo "[imp-scores] config:    $CONFIG_FILE"
echo "[imp-scores] binary:    $BIN"
echo "[imp-scores] pattern:   $PLY_PATTERN"
echo "[imp-scores] sequences: ${#SEQUENCES[@]}"

processed=0
skipped=0
failed=0

for seq in "${SEQUENCES[@]}"; do
    ply_dir="$DATASET_PATH/$seq/$PLYS_SUBFOLDER"
    scores_dir="$DATASET_PATH/$seq/$SCORES_SUBFOLDER"

    if [[ ! -d "$ply_dir" ]]; then
        echo "[imp-scores] WARN missing ply dir for $seq: $ply_dir"
        skipped=$((skipped + 1))
        continue
    fi

    # Collect the available frame indices for this sequence.
    indices=()
    while IFS= read -r ply_file; do
        base="$(basename "$ply_file")"
        idx="${base#"$PLY_PREFIX"}"
        idx="${idx%"$PLY_SUFFIX"}"
        indices+=("$idx")
    done < <(find "$ply_dir" -maxdepth 1 -name "$PLY_GLOB" | sort)

    if (( ${#indices[@]} == 0 )); then
        echo "[imp-scores] WARN no ply matching '$PLY_GLOB' in: $ply_dir"
        skipped=$((skipped + 1))
        continue
    fi

    # Apply optional start frame / frame count restrictions.
    selected=()
    for idx in "${indices[@]}"; do
        if [[ -n "$START_FRAME" ]] && (( 10#$idx < START_FRAME )); then
            continue
        fi
        selected+=("$idx")
    done
    if [[ -n "$FRAME_COUNT" ]]; then
        selected=("${selected[@]:0:$FRAME_COUNT}")
    fi

    if (( ${#selected[@]} == 0 )); then
        echo "[imp-scores] WARN no frame selected for $seq"
        skipped=$((skipped + 1))
        continue
    fi

    echo "[imp-scores] sequence: $seq (${#selected[@]} frames)"
    [[ "$DRY_RUN" -eq 1 ]] || mkdir -p "$scores_dir"

    for idx in "${selected[@]}"; do
        score_file="$scores_dir/${PLY_PREFIX}${idx}${SCORES_SUFFIX}"
        if [[ -f "$score_file" && "$OVERWRITE" -eq 0 ]]; then
            echo "[imp-scores]   skip existing: $score_file"
            skipped=$((skipped + 1))
            continue
        fi

        cmd=(
            "$BIN"
            "--source=$ply_dir/$PLY_PATTERN"
            "--startFrame=$((10#$idx))"
            "--frameCount=1"
            "--impScore=$score_file"
        )
        [[ -n "$VIEWPOINT" ]] && cmd+=("--viewpoint=$VIEWPOINT")
        [[ -n "$DOWNSCALE_FACTOR" ]] && cmd+=("--downscaleFactor=$DOWNSCALE_FACTOR")
        [[ -n "$DOWNSCALE_MIN_WIDTH" ]] && cmd+=("--downscaleMinWidth=$DOWNSCALE_MIN_WIDTH")
        [[ -n "$THREADS" ]] && cmd+=("--threads=$THREADS")
        [[ "$VERBOSE" -eq 1 ]] && cmd+=("--verbose=1")

        if [[ "$DRY_RUN" -eq 1 ]]; then
            printf '[imp-scores]   DRY-RUN'
            printf ' %q' "${cmd[@]}"
            printf '\n'
            processed=$((processed + 1))
            continue
        fi

        echo "[imp-scores]   frame $idx -> $score_file"
        if { [[ "$VERBOSE" -eq 1 ]] && "${cmd[@]}"; } || { [[ "$VERBOSE" -eq 0 ]] && "${cmd[@]}" > /dev/null; }; then
            processed=$((processed + 1))
        else
            echo "[imp-scores]   ERROR failed: $seq frame $idx" >&2
            rm -f "$score_file"
            failed=$((failed + 1))
        fi
    done
done

echo "[imp-scores] summary: processed=$processed skipped=$skipped failed=$failed"

if [[ "$failed" -gt 0 ]]; then
    exit 1
fi
