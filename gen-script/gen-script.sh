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

## Experiment name
JOBNAME="G-PCC-for-GSC-mpeg152"

## Experiment settings
CFG_DIR=cfg_CTC/JEE6.6  # cfg_CTC/JEE6.6 or cfg_CTC/JEE6.2 or cfg_3dgs_sample 
AUTO_RUN=1              # 0:off, 1:run_local, 2:run_slurm 
FLAG_1F=1               # 0:NF, 1:1F

## Dataset path
DATASET_PATH=/path/to/CTC/

## Output path
BASE_OUTPUT="/path/to/output"

## Binary and Configuration Paths
TMC13_DIR="/path/to/g-pcc-for-gsc"
TMC13_BIN="$TMC13_DIR/build/tmc3/tmc3"
CFG_PATH="$TMC13_DIR/$CFG_DIR" 
METRIC="/path/to/mpeg-gsc-metrics/build/bin/mpeg-gsc-metrics"
GSCTOOL="/path/to/mpeg-gsc-tools/"
PYTHON="python3"

## Conditions and Tools Candidates
MODEs=("lossless-geom-lossy-attrs" "lossless-geom-lossless-attrs")
TOOLs=("octree-raht" "octree-lift" "octree-pred" "predgeom-lift" "predgeom-pred" "predgeom-raht")

###########################################################################################

## Test sequences
declare -A PATHs
ADD_SEQ() { PATHs[$1]="$DATASET_PATH/$1/plys/"; }
# 1F & NF
ADD_SEQ bartender_semitracked
ADD_SEQ cinema_semitracked
ADD_SEQ breakfast_semitracked
ADD_SEQ breakfast_untracked
ADD_SEQ breakdance_untracked
ADD_SEQ manwithfruit_tracked
# NF optional
ADD_SEQ bartender_tracked 
ADD_SEQ cinema_tracked 
ADD_SEQ breakfast_tracked 
# 1F 
if [ $FLAG_1F == 1 ]; then
ADD_SEQ lego_ferrari 
ADD_SEQ lego_bugatti 
ADD_SEQ cricket_player 
ADD_SEQ plant 
ADD_SEQ solo_tango_female 
ADD_SEQ solo_tango_male 
ADD_SEQ tango_duo 
ADD_SEQ tennis_player 
ADD_SEQ library 
ADD_SEQ flowerdance 
ADD_SEQ gymnast 
fi

# resolution from colmap 
declare -A REND_RESs
# 1F & NF
REND_RESs[bartender_semitracked]="1920 1080"
REND_RESs[cinema_semitracked]="1920 1080"
REND_RESs[breakfast_semitracked]="1920 1080"
REND_RESs[breakfast_untracked]="1920 1080"
REND_RESs[breakdance_untracked]="1920 1080"
REND_RESs[bartender_tracked]="1920 1080"
REND_RESs[cinema_tracked]="1920 1080"
REND_RESs[breakfast_tracked]="1920 1080"
REND_RESs[manwithfruit_tracked]="3840 2160"
REND_RESs[lego_ferrari]="4594 5514"
REND_RESs[lego_bugatti]="3852 2868"
REND_RESs[cricket_player]="4520 2540"
REND_RESs[plant]="2954 3968"
REND_RESs[solo_tango_female]="4534 2542"
REND_RESs[solo_tango_male]="4528 2544"
REND_RESs[tango_duo]="4522 2538"
REND_RESs[tennis_player]="4518 2540"
REND_RESs[library]="3793 2131"
REND_RESs[flowerdance]="2456 2054"
REND_RESs[gymnast]="2456 2054"

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

WRITE() { echo "$2" >> "$1"; }
MKDIRCK() { [[ -d "$1" ]] || mkdir -p "$1"; }

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

for MODE in "${MODEs[@]}"; do
    for TOOL in "${TOOLs[@]}"; do
        
        [[ ! -d "$CFG_PATH/$TOOL/$MODE" ]] && continue

        if [ $FLAG_1F == 1 ]; then
            csv_file="$OUT/metrics/1F_GPCC-for-GSC-$TOOL-$MODE-$JOBNAME.csv"
        else
            csv_file="$OUT/metrics/NF_GPCC-for-GSC-$TOOL-$MODE-$JOBNAME.csv"
        fi
        echo "seq, rate, total_Byte, pos_Byte, sh0_Byte, sh1_Byte, sh2_Byte, sh3_Byte, rot_Byte, scale_Byte, opacity_Byte, RGB-PSNR, YUV-PSNR, YUV-SSIM, YUV-IVSSIM, GT-RGB-PSNR, min_RGB-PSNR, min_YUV-PSNR, min_YUV-SSIM, max_RGB-PSNR, max_YUV-PSNR, max_YUV-SSIM, total_encT, pos_encT, att_encT, total_decT, pos_decT, att_decT, enc_MaxRSS, dec_MaxRSS" > "$csv_file"

        for SEQ in "${!PATHs[@]}"; do
            PLY_PATH="${PATHs[$SEQ]}"
            files=("$PLY_PATH"/*.camerapos.ply)
            read -r REN_W REN_H <<< "${REND_RESs[$SEQ]}"
            
            # set the start frame, the number of frames, viewpoints sampling period
            if [ $FLAG_1F == 1 ]; then NUM_FRAMES=1; else NUM_FRAMES=32; fi
            if [[ "$SEQ" == "library" ]]; then VIEW_PERIOD=40; else VIEW_PERIOD=1; fi  

            SelFILEs=("${files[@]:0:$NUM_FRAMES}")
            CFG_DIRs=($CFG_PATH/$TOOL/$MODE/*)
            for CFG in "${CFG_DIRs[@]}"; do
                RATE="$(basename "$CFG")"
                RESULT_DIR="$OUT/results/$TOOL/$MODE/$SEQ/$RATE"
                MKDIRCK "$RESULT_DIR"

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
                        echo "$PYTHON $GSCTOOL/pre_post_processing/gs_pre_process.py -i \"$FILE\" -c \"$RESULT_DIR/$FRAME.json\" -o \"$RESULT_DIR/$FRAME.pre.ply\""
                        echo "make -f $TMC13_DIR/scripts/Makefile.tmc13-step -C \"$RESULT_DIR\" VPATH=\"$CFG\" ENCODER=\"$TMC13_BIN\" DECODER=\"$TMC13_BIN\" SRCSEQ=\"$RESULT_DIR/$FRAME.pre.ply\""
                        echo "$PYTHON $GSCTOOL/pre_post_processing/gs_post_process.py -i \"$RESULT_DIR/$FRAME.pre.ply.bin.decoded.ply\" -c \"$RESULT_DIR/$FRAME.json\" -o \"$RESULT_DIR/$FRAME.ply.bin.decoded.ply\""
                        echo "touch $OUT/jobs/done/gpcc/\$SLURM_JOB_ID"
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
                            echo "$METRIC -v 1 --startFrame=0 --frameCount=1 --width=$REN_W --height=$REN_H --cpu=1 --useCameraPosition=1 --onlyViewpoint=$VIEW --source=\"$FILE\" --decode=\"$RESULT_DIR/$FRAME.ply.bin.decoded.ply\" --output=\"$RESULT_DIR/$FRAME.ply.bin.decoded.ply.v$VIEW.metric\""
                            echo "touch $OUT/jobs/done/metric/\$SLURM_JOB_ID"
                        } > "$JOBFILE_MET"
                        chmod +x "$JOBFILE_MET"
                        echo "$JOBFILE_MET" >> "$RUN_METRIC"
                    done
                done

                echo "$PYTHON $SCRIPT_DIR/gsc_parse_gpcc_results.py -d \"$RESULT_DIR\" -s \"$SEQ\" -r \"$RATE\" -f $FLAG_1F >> \"$csv_file\"" >> "$RUN_SUMMARY"
            
            done
        done

        if [ $FLAG_1F == 1 ]; then
        echo "$PYTHON $SCRIPT_DIR/gsc_csv2excel.py -f 1 -t \"$SCRIPT_DIR/templates/WG04N0751 template.xlsm\" -a \"$SCRIPT_DIR/anchor_results/1F_GPCC-for-GSC-octree-raht-lossless-geom-lossy-attrs-release-mpeg152.csv\" -p \"$csv_file\" -x \"$OUT/metrics/1F_anchor__vs__GPCC-for-GSC-$TOOL-$MODE-$JOBNAME.xlsm\"" >> "$RUN_SUMMARY"
        else
        echo "$PYTHON $SCRIPT_DIR/gsc_csv2excel.py -f 0 -t \"$SCRIPT_DIR/templates/WG04N0751 template.xlsm\" -a \"$SCRIPT_DIR/anchor_results/NF_GPCC-for-GSC-octree-raht-lossless-geom-lossy-attrs-release-mpeg151.csv\" -p \"$csv_file\" -x \"$OUT/metrics/NF_anchor__vs__GPCC-for-GSC-$TOOL-$MODE-$JOBNAME.xlsm\"" >> "$RUN_SUMMARY"
        fi
    done
done

# Function to run jobs in parallel with progress display
run_with_progress() {
    local jobfile="$1"
    local total
    total=$(wc -l < "$jobfile")
    local count=0

    echo "Total jobs: $total"

    # Run jobs in parallel, emit "DONE" after each job finishes
    cat "$jobfile" | xargs -P "$(nproc)" -I {} bash -c '{} >/dev/null 2>&1; echo DONE' |
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
