# The copyright in this software is being made available under the BSD
# Licence, included below.  This software may be subject to other third
# party and contributor rights, including patent rights, and no such
# rights are granted under this licence.
#
# Copyright (c) 2025, ISO/IEC
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
#
# * Neither the name of the ISO/IEC nor the names of its contributors
#   may be used to endorse or promote products derived from this
#   software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

#!/bin/bash
set -euo pipefail

# Script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

## Helper functions for checks
ENSURE_EXEC() { [ -x "$1" ] || { echo "ERROR: Executable not found or not executable: $1"; exit 1; }; }
ENSURE_FILE() { [ -f "$1" ] || { echo "ERROR: File not found: $1"; exit 1; }; }
ENSURE_DIR() { [ -d "$1" ] || { echo "ERROR: Directory not found: $1"; exit 1; }; }
WRITE() { echo "$2" >> "$1"; }
MKDIRCK() { [[ -d "$1" ]] || mkdir -p "$1"; }

## Experiment name
JOBNAME="TM13-GS-v2"

## Encoder flags
enc_flags=""

## Experiment settings
CFG_DIR=cfg
AUTO_RUN=2              # 0:off, 1:run_local, 2:run_slurm 
FLAG_1F=1               # 0:NF, 1:1F
RUN_MODE=0              # 0 = anchor(octree-hybrid), 1 = alternative

# Whether to sample viewpoints for library sequence.
# 0:off, 1:on (set VIEW_PERIOD=40 for library sequence)
# If you use the latest version of the library sequence, set to 0 because the cameras are already sampled with a period of 40.
# If you use the older version of the library sequence, set to 1 to sample viewpoints with a period of 40.
SAMPLE_VIEWS_FOR_LIBRARY=0

## Dataset path
DATASET_PATH=/path/to/CTC/dataset
PLYS_SUBFOLDER="plys"

## Output path
BASE_OUTPUT="/path/to/output"

## Binary and Configuration Paths
TMC13_DIR="/path/to/tmc13-GS"
TMC13_BIN="$TMC13_DIR/build/tmc3/tmc3"
INPUT_SEQUENCES_CFG="${SCRIPT_DIR}/testset_inputs/dataset.config.ctc.sh"
CFG_PATH="$TMC13_DIR/$CFG_DIR" 
METRIC="/path/to/mpeg-gsc-metrics/build/bin/mpeg-gsc-metrics"
GSCTOOL="/path/to/mpeg-gsc-tools/"
PYTHON="python3"

## External Pre-Post processing
FLAG_EXT_PP=0  # 0 = No external pre/post, 1 = External pre-/post-processing
PRE_SCRIPT="/path/to/gs_pre_process.py"
POST_SCRIPT="/path/to/gs_post_process.py"

# Initial checks
ENSURE_DIR "$DATASET_PATH"
ENSURE_DIR "$TMC13_DIR"
ENSURE_DIR "$GSCTOOL"
ENSURE_EXEC "$TMC13_BIN"
ENSURE_EXEC "$METRIC"
ENSURE_FILE "$INPUT_SEQUENCES_CFG"

if [[ $RUN_MODE -eq 0 ]]; then
    MODE_TOOLs=(
        "lossless-geom-lossy-attrs octree-hybrid"
        "lossless-geom-lossless-attrs octree-predlift"
    )
else
    MODE_TOOLs=(
        "lossless-geom-lossy-attrs octree-hybrid"
        "lossless-geom-lossy-attrs octree-predlift"
        "lossless-geom-lossy-attrs octree-raht"
        "lossless-geom-lossy-attrs predgeom-hybrid"
        "lossless-geom-lossy-attrs predgeom-predlift"
        "lossless-geom-lossy-attrs predgeom-raht"
        "lossless-geom-lossless-attrs octree-predlift"
        "lossless-geom-lossless-attrs predgeom-predlift"
    )
fi

###########################################################################################
# Prepare dataset sequences
###########################################################################################
declare -A PATHs       # PLY folder paths
declare -A REND_RESs   # Render resolution: "width height"

# Function that external files will use to register sequences
ADD_SEQ() {
    local seq="$1"
    local width="$2"
    local height="$3"
    
    PATHs["$seq"]="$DATASET_PATH/$seq/$PLYS_SUBFOLDER/"
    REND_RESs["$seq"]="$width $height"

    if [[ "$seq" == "library" ]]; then
        # Check if the camerapos.ply file exists for the library sequence
        if [ -f "$DATASET_PATH/$seq/$PLYS_SUBFOLDER/frame_0000.camerapos.ply" ]; then
            camera_count=$(grep -v "camera_position" "$DATASET_PATH/$seq/$PLYS_SUBFOLDER/frame_0000.camerapos.ply" | wc -l)
            if [ "$camera_count" -gt 100 ]; then
                echo "[WARN] Please update the library dataset from 4,000 views to the 100-sampled COLMAP version on the content server. The viewpoint subsampling feature will be deprecated soon."
                if [ "$SAMPLE_VIEWS_FOR_LIBRARY" -eq 0 ]; then
                    echo "[ERROR] If you are using an older version of the library (with 4000 views), please set SAMPLE_VIEWS_FOR_LIBRARY=1."
                    exit 1
                fi
            fi
        fi
    fi
}

if [ "$FLAG_1F" -eq 1 ]; then
    PREFIX="1F"
else
    PREFIX="NF"
fi

# Load sequences from external file
source "$INPUT_SEQUENCES_CFG"

# Validate dataset paths
for SEQ in "${!PATHs[@]}"; do
    PLY_PATH="${PATHs[$SEQ]}"
    if [[ ! -d "$PLY_PATH" ]]; then
        echo "❌ $PLY_PATH does not exist"
        exit 1
    fi
    shopt -s nullglob
    files=("$PLY_PATH"/*.camerapos.ply)
    if (( ${#files[@]} < 1 )); then
        echo "❌ No *.camerapos.ply in $PLY_PATH"
        exit 1
    fi
done

SCRIPT_PATH="$(readlink -f "$0" 2>/dev/null || realpath "$0")"
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"

CURTIME="$(date "+%y-%m-%d-%H%M%S")"
OUT="$BASE_OUTPUT/$JOBNAME-$CURTIME"
BATCH_GPCC="$OUT/jobs/batch_gpcc"
BATCH_METRIC="$OUT/jobs/batch_metric"

MKDIRCK "$OUT/metrics"
MKDIRCK "$BATCH_GPCC"
MKDIRCK "$BATCH_METRIC"
MKDIRCK "$OUT/bin"
MKDIRCK "$OUT/jobs/logs"
MKDIRCK "$OUT/jobs/done/gpcc"
MKDIRCK "$OUT/jobs/done/metric"

# Validate build path
if [[ ! -d "$TMC13_DIR/build" ]]; then
    echo "❌ $TMC13_DIR/build does not exist"
    exit 1
fi

## automatic cfg generation
mkdir -p "$TMC13_DIR/$CFG_DIR" 
cd "$TMC13_DIR/$CFG_DIR"
if [ $FLAG_1F == 1 ]; then
    if [ $RUN_MODE == 0 ]; then
    bash "$TMC13_DIR/scripts/gen-cfg.sh" --octree --hybrid --gaussiansplat --1F --cfgdir=$TMC13_DIR/cfg
    bash "$TMC13_DIR/scripts/gen-cfg.sh" --octree --pred-lift --gaussiansplat --1F --cfgdir=$TMC13_DIR/cfg
    else 
    bash "$TMC13_DIR/scripts/gen-cfg.sh" --all --gaussiansplat --1F --cfgdir=$TMC13_DIR/cfg
    fi
else
    if [ $RUN_MODE == 0 ]; then
    bash "$TMC13_DIR/scripts/gen-cfg.sh" --octree --hybrid --gaussiansplat --NF --cfgdir=$TMC13_DIR/cfg
    bash "$TMC13_DIR/scripts/gen-cfg.sh" --octree --pred-lift --gaussiansplat --NF --cfgdir=$TMC13_DIR/cfg
    else 
    bash "$TMC13_DIR/script/gen-cfg.sh" --all --gaussiansplat --NF --cfgdir=$TMC13_DIR/cfg
    fi
fi
ENSURE_DIR "$CFG_PATH"

# Build binaries
cd "$TMC13_DIR/build"
make -j"$(nproc)"
cp "$TMC13_BIN" "$OUT/bin/tmc3"
TMC13_BIN="$OUT/bin/tmc3"

cd "$OUT/jobs"
RUN_GPCC="$OUT/jobs/step1_run_gpcc.sh"
RUN_METRIC="$OUT/jobs/step2_run_metric.sh"
RUN_SUMMARY="$OUT/jobs/step3_run_summary.sh"

# Initialize run scripts
: > "$RUN_GPCC"
: > "$RUN_METRIC"
: > "$RUN_SUMMARY"

for MODE_TOOL in "${MODE_TOOLs[@]}"; do
    read -r MODE TOOL <<< "$MODE_TOOL"

    [[ ! -d "$CFG_PATH/$TOOL/$MODE" ]] && continue

    csv_file="$OUT/metrics/${PREFIX}_TMC13-GS-$TOOL-$MODE-$JOBNAME.csv"
    echo "seq, rate, total_Byte, pos_Byte, sh0_Byte, sh1_Byte, sh2_Byte, sh3_Byte, rot_Byte, scale_Byte, opacity_Byte, RGB-PSNR, YUV-PSNR, YUV-SSIM, YUV-IVSSIM, GT-RGB-PSNR, min_RGB-PSNR, min_YUV-PSNR, min_YUV-SSIM, max_RGB-PSNR, max_YUV-PSNR, max_YUV-SSIM, total_encT, pos_encT, att_encT, total_decT, pos_decT, att_decT, enc_MaxRSS, dec_MaxRSS" > "$csv_file"

    for SEQ in "${!PATHs[@]}"; do
        PLY_PATH="${PATHs[$SEQ]}"
        files=("$PLY_PATH"/*.camerapos.ply)
        read -r REN_W REN_H <<< "${REND_RESs[$SEQ]}"
        
        # set the start frame, the number of frames, viewpoints sampling period
        if [ $FLAG_1F == 1 ]; then NUM_FRAMES=1; else NUM_FRAMES=32; fi
        # For only the old version of library sequence (SAMPLE_VIEWS_FOR_LIBRARY=1), we sample viewpoints with a period of 40. 
        if [ $SAMPLE_VIEWS_FOR_LIBRARY == 1 ] && [[ "$SEQ" == "library" ]]; then VIEW_PERIOD=40; else VIEW_PERIOD=1; fi

        SelFILEs=("${files[@]:0:$NUM_FRAMES}")
        if [[ "$MODE" == "lossless-geom-lossless-attrs" ]]; then
            CFG_DIRs=("$CFG_PATH/$TOOL/$MODE/$SEQ")
        else
            CFG_DIRs=($CFG_PATH/$TOOL/$MODE/$SEQ/*)
        fi
        for CFG in "${CFG_DIRs[@]}"; do
            if [[ "$MODE" == "lossless-geom-lossless-attrs" ]]; then
                RATE="lossless"   
                RESULT_DIR="$OUT/results/$TOOL/$MODE/$SEQ"
                MKDIRCK "$RESULT_DIR"
            else
                RATE="$(basename "$CFG")"
                RESULT_DIR="$OUT/results/$TOOL/$MODE/$SEQ/$RATE"
                MKDIRCK "$RESULT_DIR"
            fi

            echo "$TOOL/$MODE/$SEQ/$RATE"

            for FILE in "${SelFILEs[@]}"; do
                FRAME="$(basename "$FILE" .camerapos.ply)"
                JOBFILE="$BATCH_GPCC/$JOBNAME-$CURTIME-$TOOL-$MODE-$SEQ-$RATE-$FRAME.sh"

                {
                    echo "#!/bin/bash"
                    echo "#SBATCH --output=$OUT/jobs/logs/%j.out"
                    echo "#SBATCH --ntasks-per-node=1"
                    echo "#SBATCH --nodes=1"
                    echo "#SBATCH --mem-per-cpu=8G"
                    echo "set -euo pipefail"
                    SRC_PLY="${FILE/.camerapos/}"
                    if [ $FLAG_EXT_PP == 1 ]; then
                        echo "$PYTHON ${PRE_SCRIPT} -i \"$SRC_PLY\" -c \"$RESULT_DIR/$FRAME.json\" -o \"$RESULT_DIR/$FRAME.pre.ply\""
                        echo "make -f $TMC13_DIR/scripts/Makefile.tmc13-step -C \"$RESULT_DIR\" VPATH=\"$CFG\" ENCODER=\"$TMC13_BIN\" DECODER=\"$TMC13_BIN\" SRCSEQ=\"$RESULT_DIR/$FRAME.pre.ply\" ENCFLAGS=\"${enc_flags}\""
                        echo "$PYTHON \"$POST_SCRIPT\" -i \"$RESULT_DIR/$FRAME.pre.ply.bin.decoded.ply\" -c \"$RESULT_DIR/$FRAME.json\" -o \"$RESULT_DIR/$FRAME.ply.bin.decoded.ply\""
                    else
                        echo "make -f $TMC13_DIR/scripts/Makefile.tmc13-step -C \"$RESULT_DIR\" VPATH=\"$CFG\" ENCODER=\"$TMC13_BIN\" DECODER=\"$TMC13_BIN\" SRCSEQ=\"$SRC_PLY\" ENCFLAGS=\"${enc_flags}\""
                    fi
                    echo "touch $OUT/jobs/done/gpcc/\${SLURM_JOB_ID:-UNDEFINED}"
                } > "$JOBFILE"
                chmod +x "$JOBFILE"
                echo "$JOBFILE" >> "$RUN_GPCC"

                numVIEWs=$(grep -a "camera_position" "$FILE" | grep -v "camera_position_header" | wc -l)
                for (( VIEW=0; VIEW<numVIEWs; VIEW+=VIEW_PERIOD )); do
                    JOBFILE_MET="$BATCH_METRIC/$JOBNAME-$CURTIME-$TOOL-$MODE-$SEQ-$RATE-$FRAME-$VIEW.sh"
                    {
                        echo "#!/bin/bash"
                        echo "#SBATCH --output=$OUT/jobs/logs/%j.out"
                        echo "#SBATCH --ntasks-per-node=1"
                        echo "#SBATCH --nodes=1"
                        echo "#SBATCH --mem-per-cpu=8G"
                        echo "set -euo pipefail"
                        echo "$METRIC -v 1 --startFrame=0 --frameCount=1 --width=$REN_W --height=$REN_H --cpu=1 --threads=1 --useCameraPosition=1 --onlyViewpoint=$VIEW --computeSsim=0 --computeOm=0 --source=\"$FILE\" --decode=\"$RESULT_DIR/$FRAME.ply.bin.decoded.ply\" --output=\"$RESULT_DIR/$FRAME.ply.bin.decoded.ply.v$VIEW.metric\""
                        echo "touch $OUT/jobs/done/metric/\${SLURM_JOB_ID:-UNDEFINED}"
                    } > "$JOBFILE_MET"
                    chmod +x "$JOBFILE_MET"
                    echo "$JOBFILE_MET" >> "$RUN_METRIC"
                done
            done
            echo "$PYTHON $SCRIPT_DIR/gsc_parse_gpcc_results.py -d \"$RESULT_DIR\" -s \"$SEQ\" -r \"$RATE\" -f $FLAG_1F >> \"$csv_file\"" >> "$RUN_SUMMARY"
        
        done
    done
done

if [[ $RUN_MODE -eq 0 ]]; then
    read -r MODE TOOL <<< "${MODE_TOOLs[0]}"
    csv_C1="$OUT/metrics/${PREFIX}_TMC13-GS-$TOOL-$MODE-$JOBNAME.csv"

    read -r MODE TOOL <<< "${MODE_TOOLs[1]}"
    csv_CW="$OUT/metrics/${PREFIX}_TMC13-GS-$TOOL-$MODE-$JOBNAME.csv"

    echo "$PYTHON $SCRIPT_DIR/gsc_csv2excel.py -f $FLAG_1F \
        -t \"$SCRIPT_DIR/templates/WG04N0826_template.xlsm\" \
        -x \"$OUT/metrics/${PREFIX}_CTC_$JOBNAME.xlsm\" \
        --aC1 \"$SCRIPT_DIR/anchor_results/${PREFIX}_TMC13-GS-octree-hybrid-lossless-geom-lossy-attrs.csv\" \
        --pC1 \"$csv_C1\" \
        --aCW \"$SCRIPT_DIR/anchor_results/${PREFIX}_TMC13-GS-octree-predlift-lossless-geom-lossless-attrs.csv\" \
        --pCW \"$csv_CW\"" >> "$RUN_SUMMARY"

else
    for MODE_TOOL in "${MODE_TOOLs[@]}"; do
        read -r MODE TOOL <<< "${MODE_TOOL}"
        csv_file="$OUT/metrics/${PREFIX}_TMC13-GS-$TOOL-$MODE-$JOBNAME.csv"
        
        CMD_BASE="$PYTHON $SCRIPT_DIR/gsc_csv2excel.py -f $FLAG_1F \
            -t \"$SCRIPT_DIR/templates/WG04N0826_template.xlsm\" \
            -x \"$OUT/metrics/${PREFIX}_CTC_$JOBNAME.xlsm\""

        if [[ "$MODE" == "lossless-geom-lossless-attrs" ]]; then
            echo "$CMD_BASE \
                --aCW \"$SCRIPT_DIR/anchor_results/${PREFIX}_TMC13-GS-octree-predlift-lossless-geom-lossless-attrs.csv\" \
                --pCW \"$csv_file\"" >> "$RUN_SUMMARY"
        else
            echo "$CMD_BASE \
                --aC1 \"$SCRIPT_DIR/anchor_results/${PREFIX}_TMC13-GS-octree-hybrid-lossless-geom-lossy-attrs.csv\" \
                --pC1 \"$csv_file\"" >> "$RUN_SUMMARY"
        fi
    done
fi

# Function to run jobs in parallel with progress display
run_with_progress() {
    local jobfile="$1"
    local nproc="${2:-$(nproc)}"
    local total
    total=$(wc -l < "$jobfile")
    local count=0

    echo "Total jobs: $total"

    # Run jobs in parallel, emit "DONE" after each job finishes
    cat "$jobfile" | xargs -P "${nproc}" -I {} bash -c '{} >/dev/null 2>&1; echo DONE' |
    while read -r signal; do
        count=$((count+1))
        percent=$((count*100/total))
        echo -ne "Progress: $count/$total (${percent}%)\r"
    done
    echo
}

run_with_progress_sbatch_wait() {
    local jobfile="$1"
    local done_dir="$2"
    mkdir -p "$done_dir"

    local total
    total=$(wc -l < "$jobfile")
    local submitted=0
    local done_count=0
    echo "Submitting $total jobs to SLURM..."
    while read -r job; do
        while (( submitted - done_count > 10000 )); do
            sleep 4
            submit_percent=$((submitted * 100 / total))
            done_count=$(find "$done_dir" -type f -printf '.' | wc -c)
            done_percent=$((done_count * 100 / total))
            echo -ne "Submitted: $submitted/$total (${submit_percent}%). Progress: $done_count/$total (${done_percent}%)\r"
        done

        sbatch "$job" >/dev/null
        submitted=$((submitted + 1))
        submit_percent=$((submitted * 100 / total))

        if (( submitted % 100 == 0 )); then
            done_count=$(find "$done_dir" -type f -printf '.' | wc -c)
        fi

        done_percent=$((done_count * 100 / total))
        echo -ne "Submitted: $submitted/$total (${submit_percent}%). Progress: $done_count/$total (${done_percent}%)\r"
    done < "$jobfile"

    echo "All jobs submitted. Waiting for completion...                              "
    sleep 2

    while (( done_count < total )); do
        done_count=$(find "$done_dir" -type f -printf '.' | wc -c)
        percent=$((done_count*100/total))
        echo -ne "Progress: $done_count/$total (${percent}%)\r"
        sleep 1
    done
    echo
}

## Master run script
MASTER_RUN_SLURM="$OUT/jobs/run_all_slurm.sh"
MASTER_RUN="$OUT/jobs/run_all.sh"

# Export function for master script
declare -f run_with_progress_sbatch_wait > "$MASTER_RUN_SLURM"
{
    echo "#!/bin/bash"
    echo "set -euo pipefail"
    echo
    echo "SCRIPT_PATH=\"$(readlink -f "$0" 2>/dev/null || realpath \"$0\")\""
    echo "SCRIPT_DIR=\"$(dirname "$SCRIPT_PATH")\""
    echo "cd $SCRIPT_DIR"
    echo
    echo "echo '▶ Step 1: GPCC encode/decode in parallel'"
    echo "run_with_progress_sbatch_wait \"$RUN_GPCC\" \"$OUT/jobs/done/gpcc/\" "
    echo
    echo "echo '▶ Step 2: Metric calculation in parallel'"
    echo "run_with_progress_sbatch_wait \"$RUN_METRIC\" \"$OUT/jobs/done/metric/\""
    echo
    echo "echo '▶ Step 3: Summary generation'"
    echo "bash \"$RUN_SUMMARY\""
    echo "echo '✅ All tasks completed'"
} >> "$MASTER_RUN_SLURM"
chmod +x "$MASTER_RUN_SLURM"

# Export function for master script
declare -f run_with_progress > "$MASTER_RUN"
{
    echo "#!/bin/bash"
    echo "set -euo pipefail"
    echo
    echo "SCRIPT_PATH=\"$(readlink -f "$0" 2>/dev/null || realpath \"$0\")\""
    echo "SCRIPT_DIR=\"$(dirname "$SCRIPT_PATH")\""
    echo "cd $SCRIPT_DIR"
    echo
    echo "echo '▶ Step 1: GPCC encode/decode in parallel'"
    echo "run_with_progress \"$RUN_GPCC\""
    echo
    echo "echo '▶ Step 2: Metric calculation in parallel'"
    echo "run_with_progress \"$RUN_METRIC\""
    echo
    echo "echo '▶ Step 3: Summary generation'"
    echo "bash \"$RUN_SUMMARY\""
    echo "echo '✅ All tasks completed'"
} >> "$MASTER_RUN"
chmod +x "$MASTER_RUN"

echo
echo "[Generated scripts]"
echo "$RUN_GPCC"
echo "$RUN_METRIC"
echo "$RUN_SUMMARY"
echo "$MASTER_RUN"
echo "$MASTER_RUN_SLURM"
echo
echo "💡 Run (Local): bash \"$MASTER_RUN\""
echo "💡 Run (Slurm): bash \"$MASTER_RUN_SLURM\""

if [ $AUTO_RUN == 1 ]; then
bash "$MASTER_RUN"
elif [ $AUTO_RUN == 2 ]; then
bash "$MASTER_RUN_SLURM"
fi
