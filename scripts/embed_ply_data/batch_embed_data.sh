#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

usage() {
    cat <<'EOF'
Usage:
    batch_embed_data.sh \
    --dataset-path <path> \
    --config <dataset.config.*.sh> \
    --gstools-bin <path-to-gstools-binary> \
    [--plys-subfolder <name>] \
    [--colmap-subfolder <name>] \
    [--overwrite] \
    [--embed-scores] \
    [--scores-only] \
    [--scores-subfolder <name>] \
    [--scores-suffix <suffix>] \
    [--scores-output-suffix <suffix>] \
    [--camera-output-suffix <suffix>] \
    [--camera-score-output-suffix <suffix>] \
    [--score-name <attribute-name>] \
    [--score-percentile <float>] \
    [--no-score-log-scale] \
    [--frame-count <N>] \
    [--python <python-executable>]

Description:
  Reads sequences from a dataset config file (same format used by gen-script.sh,
    i.e. ADD_SEQ <sequence> <width> <height>) and runs embed_data.py on each
  sequence PLY folder.

Modes:
    1) Cameras only (default)
    2) Scores only (--scores-only)
    3) Scores + cameras (--embed-scores)

Expected per-sequence layout:
  <dataset-path>/<sequence>/<plys-subfolder>/frame_XXXX.ply
  <dataset-path>/<sequence>/<colmap-subfolder>/frame_XXXX/sparse/{cameras.bin,images.bin}

Examples:
    Cameras only:
        ./batch_embed_data.sh \
            --dataset-path /media/gs_storage/DATASET/MPEG/CTC_156 \
                --config ./testset_inputs/dataset.config.ctc.sh \
            --gstools-bin /mnt/data01/mpeg_suite/mpeg-gsc-tools/gsTools/build/msvc/Release/bin/gstools \
            --overwrite

    Scores only:
        ./batch_embed_data.sh \
            --dataset-path <DATASET_ROOT> \
            --config <PATH_TO_DATASET_CONFIG_SH> \
            --gstools-bin <PATH_TO_GSTOOLS_BINARY> \
            --scores-only \
            --scores-subfolder scores \
            --scores-suffix .impScores.txt \
            --scores-output-suffix .score.ply \
            --score-name score \
            --score-percentile 1.0 \
            --overwrite

    Camera-only output suffix:
        ./batch_embed_data.sh \
            --dataset-path <DATASET_ROOT> \
            --config <PATH_TO_DATASET_CONFIG_SH> \
            --gstools-bin <PATH_TO_GSTOOLS_BINARY> \
            --camera-output-suffix .camerapos.ply \
            --overwrite

    Camera+scores output suffix:
        ./batch_embed_data.sh \
            --dataset-path <DATASET_ROOT> \
            --config <PATH_TO_DATASET_CONFIG_SH> \
            --gstools-bin <PATH_TO_GSTOOLS_BINARY> \
            --embed-scores \
            --scores-subfolder scores \
            --scores-suffix .impScores.txt \
            --camera-score-output-suffix .camerapos.ply \
            --score-name score \
            --score-percentile 1.0 \
            --overwrite

    Scores and Cameras:
        ./batch_embed_data.sh \
            --dataset-path /media/gs_storage/DATASET/MPEG/CTC_156 \
            --config ./testset_inputs/dataset.config.ctc.sh \
            --gstools-bin /mnt/data01/mpeg_suite/mpeg-gsc-tools/gsTools/build/msvc/Release/bin/gstools \
            --embed-scores \
            --scores-subfolder scores \
            --scores-suffix .impScores.txt \
            --scores-output-suffix .camerapos.ply \
            --score-name score \
            --score-percentile 1.0 \
            --overwrite

Procedure:
    1) For each sequence, the script will look for PLY files in the specified PLY subfolder.
    2) For each PLY file, the script will look for a corresponding COLMAP sparse reconstruction folder
       (frame_XXXX/sparse) in the specified COLMAP subfolder.
    3) The script will look for a score file in the specified scores subfolder, 
       with the same base name as the PLY file and the specified scores suffix.
    4) The script will run embed_data.py with the specified options, 
       which will embed camera positions and/or scores into the PLY files.

Score files:
    Text files with scores are expected to be in the scores subfolder, with the same base name 
    as the PLY files, e.g.:
        <dataset-path>/<sequence>/<scores-subfolder>/frame_XXXX.impScores.txt
    
    There must be a score file for each PLY file, and the number of scores in the file must match 
    the number of points in the PLY file.

The resulting PLY files with embedded scores and cameras will be written to the same PLY folder, 
    with the specified output suffix, e.g.:
        <dataset-path>/<sequence>/<plys-subfolder>/frame_XXXX.camerapos.ply.

EOF
}

DATASET_PATH=""
CONFIG_FILE=""
GSTOOLS_BIN=""
PLYS_SUBFOLDER="plys"
COLMAP_SUBFOLDER="colmap_data"
OVERWRITE=0
EMBED_SCORES=0
SCORES_ONLY=0
SCORES_SUBFOLDER="scores"
SCORES_SUFFIX=".impScores.txt"
SCORES_OUTPUT_SUFFIX=".with_score.ply"
CAMERA_OUTPUT_SUFFIX=".camerapos.ply"
CAMERA_SCORE_OUTPUT_SUFFIX=".camerapos.ply"
SCORE_NAME="score"
SCORE_PERCENTILE="1.0"
SCORE_LOG_SCALE=1
FRAME_COUNT=""
PYTHON_BIN="python3"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dataset-path)
            DATASET_PATH="$2"
            shift 2
            ;;
        --config)
            CONFIG_FILE="$2"
            shift 2
            ;;
        --gstools-bin|--camerapos-bin)
            GSTOOLS_BIN="$2"
            shift 2
            ;;
        --plys-subfolder)
            PLYS_SUBFOLDER="$2"
            shift 2
            ;;
        --colmap-subfolder)
            COLMAP_SUBFOLDER="$2"
            shift 2
            ;;
        --overwrite)
            OVERWRITE=1
            shift
            ;;
        --embed-scores)
            EMBED_SCORES=1
            shift
            ;;
        --scores-only)
            SCORES_ONLY=1
            shift
            ;;
        --scores-subfolder)
            SCORES_SUBFOLDER="$2"
            shift 2
            ;;
        --scores-suffix)
            SCORES_SUFFIX="$2"
            shift 2
            ;;
        --scores-output-suffix)
            SCORES_OUTPUT_SUFFIX="$2"
            shift 2
            ;;
        --camera-output-suffix)
            CAMERA_OUTPUT_SUFFIX="$2"
            shift 2
            ;;
        --camera-score-output-suffix)
            CAMERA_SCORE_OUTPUT_SUFFIX="$2"
            shift 2
            ;;
        --score-name)
            SCORE_NAME="$2"
            shift 2
            ;;
        --score-percentile)
            SCORE_PERCENTILE="$2"
            shift 2
            ;;
        --no-score-log-scale)
            SCORE_LOG_SCALE=0
            shift
            ;;
        --frame-count)
            FRAME_COUNT="$2"
            shift 2
            ;;
        --python)
            PYTHON_BIN="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: Unknown argument: $1" >&2
            usage
            exit 1
            ;;
    esac
done

if [[ -z "$DATASET_PATH" || -z "$CONFIG_FILE" || -z "$GSTOOLS_BIN" ]]; then
    echo "ERROR: --dataset-path, --config and --gstools-bin are required." >&2
    usage
    exit 1
fi

if [[ ! -d "$DATASET_PATH" ]]; then
    echo "ERROR: dataset path does not exist: $DATASET_PATH" >&2
    exit 1
fi

if [[ ! -f "$CONFIG_FILE" ]]; then
    fallback_config="${CONFIG_FILE/ctyc/ctc}"
    if [[ "$fallback_config" != "$CONFIG_FILE" && -f "$fallback_config" ]]; then
        echo "[embed-cameras] WARN config not found, using fallback: $fallback_config"
        CONFIG_FILE="$fallback_config"
    else
        echo "ERROR: config file does not exist: $CONFIG_FILE" >&2
        exit 1
    fi
fi

if [[ ! -f "$GSTOOLS_BIN" ]]; then
    echo "ERROR: gstools binary does not exist: $GSTOOLS_BIN" >&2
    exit 1
fi

if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
    echo "ERROR: python executable not found in PATH: $PYTHON_BIN" >&2
    exit 1
fi

EMBED_SCRIPT="$SCRIPT_DIR/embed_data.py"
if [[ ! -f "$EMBED_SCRIPT" ]]; then
    echo "ERROR: embed script not found: $EMBED_SCRIPT" >&2
    exit 1
fi

declare -A PATHS
ADD_SEQ() {
    local seq="$1"
    PATHS["$seq"]="$DATASET_PATH/$seq/$PLYS_SUBFOLDER"
}

# shellcheck source=/dev/null
source "$CONFIG_FILE"

if (( ${#PATHS[@]} == 0 )); then
    echo "ERROR: no sequences were registered by config: $CONFIG_FILE" >&2
    exit 1
fi

echo "[embed-cameras] dataset: $DATASET_PATH"
echo "[embed-cameras] config:  $CONFIG_FILE"
echo "[embed-cameras] total sequences: ${#PATHS[@]}"
if [[ "$SCORES_ONLY" -eq 1 ]]; then
    echo "[embed-cameras] mode: scores-only"
fi
if [[ "$EMBED_SCORES" -eq 1 || "$SCORES_ONLY" -eq 1 ]]; then
    echo "[embed-cameras] score embedding: enabled"
    echo "[embed-cameras] scores subfolder: $SCORES_SUBFOLDER"
    echo "[embed-cameras] scores suffix:    $SCORES_SUFFIX"
    if [[ "$SCORES_ONLY" -eq 1 ]]; then
        echo "[embed-cameras] scores output suffix: $SCORES_OUTPUT_SUFFIX"
    fi
    echo "[embed-cameras] score attribute:  $SCORE_NAME"
    echo "[embed-cameras] score percentile: $SCORE_PERCENTILE"
    if [[ "$SCORE_LOG_SCALE" -eq 1 ]]; then
        echo "[embed-cameras] score log-scale:  enabled"
    else
        echo "[embed-cameras] score log-scale:  disabled"
    fi
fi

processed=0
skipped=0
failed=0

for seq in "${!PATHS[@]}"; do
    ply_dir="${PATHS[$seq]%/}"
    colmap_root="$DATASET_PATH/$seq/$COLMAP_SUBFOLDER"
    scores_root="$DATASET_PATH/$seq/$SCORES_SUBFOLDER"

    if [[ ! -d "$ply_dir" ]]; then
        echo "[embed-cameras] WARN missing ply dir for $seq: $ply_dir"
        skipped=$((skipped + 1))
        continue
    fi

    echo "[embed-cameras] sequence: $seq"

    cmd=(
        "$PYTHON_BIN" "$EMBED_SCRIPT"
        --ply-dir "$ply_dir"
        --camerapos-bin "$GSTOOLS_BIN"
        --colmap-root "$colmap_root"
        --camera-output-suffix "$CAMERA_OUTPUT_SUFFIX"
        --camera-score-output-suffix "$CAMERA_SCORE_OUTPUT_SUFFIX"
    )

    if [[ "$EMBED_SCORES" -eq 1 || "$SCORES_ONLY" -eq 1 ]]; then
        cmd+=(
            --embed-scores
            --scores-root "$scores_root"
            --scores-suffix "$SCORES_SUFFIX"
            --scores-output-suffix "$SCORES_OUTPUT_SUFFIX"
            --score-name "$SCORE_NAME"
            --score-percentile "$SCORE_PERCENTILE"
        )
        if [[ "$SCORES_ONLY" -eq 1 ]]; then
            cmd+=(--scores-only)
        fi
        if [[ "$SCORE_LOG_SCALE" -eq 0 ]]; then
            cmd+=(--no-score-log-scale)
        fi
    fi

    if [[ -n "$FRAME_COUNT" ]]; then
        cmd+=(--frame-count "$FRAME_COUNT")
    fi

    if [[ "$OVERWRITE" -eq 1 ]]; then
        cmd+=(--overwrite)
    fi

    if "${cmd[@]}"; then
        processed=$((processed + 1))
    else
        echo "[embed-cameras] ERROR failed sequence: $seq" >&2
        failed=$((failed + 1))
    fi
done

echo "[embed-cameras] summary: processed=$processed skipped=$skipped failed=$failed"

if [[ "$failed" -gt 0 ]]; then
    exit 1
fi
