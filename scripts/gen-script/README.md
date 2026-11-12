# The gen-script package for TMC13-GS simulations

This guide explains how to prepare data, configure paths, and run scripts for JEE6.2, JEE6.6 and G-PCC CEs experiments using **TMC13-GS**.

### Package Structure
The `gen-script` package is organized as follows:
```
/gen-script
├── gen-script.sh                      # Main script to generate experiment scripts
├── gsc_parse_gpcc_results.py          # Parses GPCC output into CSV
├── gsc_csv2excel.py                   # Export CSV to XLSM
├── README.md                          # Documentation for the package
│
├── anchor_results/                    # Reference (anchor) results for comparison
│   ├── 1F_TMC13-GS-octree-hybrid-lossless-geom-lossy-attrs.csv
│   ├── 1F_TMC13-GS-octree-predlift-lossless-geom-lossless-attrs.csv
│   ├── NF_TMC13-GS-octree-hybrid-lossless-geom-lossy-attrs.csv
│   └── NF_TMC13-GS-octree-predlift-lossless-geom-lossless-attrs.csv
│
├── templates/                         # XLSM Templates from WG05N0420 for Reporting
│   └── WG05N0420_template.xlsm
│
└── testset_inputs/
    └── dataset.config.ctc.sh          # Dataset configuration for evaluation

```
> Note: The provided anchor result files (.csv) were generated in the LGE environment. **After completing the anchor cross‑check, please update these CSV files accordingly**.

> Note: [MDS26500_WG05_N00420	
](https://dms.mpeg.expert/doc_end_user/documents/155_Geneva/wg11/MDS26500_WG05_N00420-v2.zip) - Common test conditions for Gaussian splat coding



## 1. Preparation 

### MPEG-GSC-Metrics & MPEG-GSC-Tools Software
Clone the repositories:
```
git clone https://git.mpeg.expert/MPEG/Explorations/GSC/gsc-software/mpeg-gsc-metrics
git clone https://git.mpeg.expert/MPEG/Explorations/GSC/gsc-software/mpeg-gsc-tools
```
Build each project separately.
### PLY Files with Camera Positions and Importance Scores

#### 1. Embed Camera Positions into PLY Files (`*.ply`)

Embed the camera parameters into the source PLY files (`*.ply`) using using **embed_ply_data** located in:

```text
scripts/embed_ply_data
```
An example script is provided in `run-embed-camera.sh`.

```bash
bash ./batch_embed_data.sh \
    --dataset-path /path/to/dataset \
    --config /path/to/dataset.config.sh \
    --gstools-bin /path/to/gstools \
    --embed-scores \
    --overwrite
```

#### 2. Generate Importance Score Files

Generate importance score files (`*.impScores.txt`) for all PLY files using **batch_importance_scores** located in:

```text
scripts/batch_importance_scores
```

An example script is provided in `run.sh`.

```bash
cd scripts/batch_importance_scores

bash ./batch_importance_scores.sh \
  --dataset-path "/path/to/<dataset-root>/" \
  --config ./testset_inputs/dataset.config.ctc.sh \
  --bin /path/to/mpeg-gsc-importance/build/Release/bin/importance-scores \
  --plys-subfolder plys
```

#### 3. Embed Importance Scores into PLY Files (`*.camerapos.ply`)

Embed the generated importance scores into the camera-position-enabled PLY files (`*.camerapos.ply`) using **embed_ply_data** located in:

```text
scripts/embed_ply_data
```
An example script is provided in `run-embed-score.sh`.

```bash
cd scripts/batch_importance_scores

bash ./batch_embed_data.sh \
  --dataset-path "/path/to/<dataset-root>/" \
  --config ./testset_inputs/dataset.config.ctc.sh \
  --gstools-bin "/path/to/mpeg-gsc-tools/gsTools/build/msvc/Release/bin/gstools" \
  --embed-scores \
  --scores-subfolder scores \
  --scores-suffix .impScores.txt \
  --scores-output-suffix .camerapos.ply \
  --score-name score \
  --score-percentile 1.0 \
  --overwrite
```

### Python Package Dependency
Install the `openpyxl` Python package (version 3.1.5) to handle XLSM files:
```
python3 -m pip install openpyxl==3.1.5
```

## 2. Build TMC13-GS
Clone the repository:
```bash
git clone https://git.mpeg.expert/MPEG/3dgh/g-pcc/software/tm/mpeg-pcc-tmc13-gs
cd mpeg-pcc-tmc13-gs
git fetch origin tag v4.0
git checkout tags/v4.0
```
Build the project:
```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```

## 3. Experiment and Path Settings
Edit `gen-script.sh` to set the sequence, binary, and configuration directories.

### Experiment Name
```bash
JOBNAME="TM13-GS-v4.0"
```
### Experiment Settings
For 1F evaluation:
```
AUTO_RUN=0              # 0:off, 1:run_local, 2:run_slurm
FLAG_1F=1               # 0:NF, 1:1F
FLAG_CONCAT=0           # 0:per-frame (standard), 1:ConcatGPCC (NF frames concatenated)
```
For NF evaluation:
```
AUTO_RUN=0              # 0:off, 1:run_local, 2:run_slurm
FLAG_1F=0               # 0:NF, 1:1F
FLAG_CONCAT=1           # 0:per-frame (standard), 1:ConcatGPCC (NF frames concatenated)
```
> **Note**: If `AUTO_RUN` is set to 1 or 2, the generated scripts will run immediately after creation.

### Test Sequences
Set the `DATASET_PATH` path:
```bash
## Dataset path
DATASET_PATH=/path/to/CTC/
```
> **Note**: `DATASET_PATH` must follow MPEG Datasets CTC (https://content.mpeg.expert/data/Explorations/GSC/CTC).

### Output Path
```bash
BASE_OUTPUT="/path/to/output"
```

### Binary and Configuration Paths
```bash
TMC13_DIR="/path/to/mpeg-pcc-tmc13-gs"     # TMC13 directory
TMC13_BIN="$TMC13_DIR/build/tmc3/tmc3" # TMC13 binary
CFG_PATH="$TMC13_DIR/cfg"              # Configuration directory
METRIC="/path/to/mpeg-gsc-metrics/build/bin/mpeg-gsc-metrics"  # Metric binary
PYTHON="python3"
```

## 4. Running
### Generate All Job Scripts
```bash
bash gen-script.sh
```
### Example Output
```bash
/path/to/output/jobs/step1_run_gpcc.sh
/path/to/output/jobs/step2_run_metric.sh
/path/to/output/jobs/step3_run_summary.sh
/path/to/output/jobs/run_all.sh
/path/to/output/jobs/run_all_slurm.sh

💡 Run in parallel (local): bash "/path/to/output/jobs/run_all.sh"
💡 Run in parallel (Slurm): bash "/path/to/output/jobs/run_all_slurm.sh"
```

When script generation is complete, the following scripts will be created in the `output/jobs` directory:

* `step0_preprocess.sh` – Pre-processing jobs
* `step1_run_gpcc.sh` – GPCC encode/decode jobs
* `step2_run_metric.sh` – Metric calculation jobs
* `step3_run_summary.sh` – Summary results
* `run_all.sh` – Run all jobs in parallel locally
* `run_all_slurm.sh` – Run all jobs in parallel on Slurm

## 5. Execution
Local execution:
```bash
bash /path/to/output/jobs/run_all.sh
```
Slurm execution:
```bash
bash /path/to/output/jobs/run_all_slurm.sh
```

## 6. `run_all.sh` Example Log
```
▶ GPCC encode/decode in parallel
Submitting 100 jobs to SLURM...
All jobs submitted. Waiting for completion...                              
Progress: 100/100 (100%)
▶ Metric calculation in parallel
Submitting 4930 jobs to SLURM...
All jobs submitted. Waiting for completion...                              
Progress: 4930/4930 (100%)
▶ Summary generation
✅ All tasks completed
```

## 7. Simulation Results
All simulation outputs are organized under:
```
/path/to/output/
└── TM13-GS-v4.0/
    ├── bin/                                           # Executables or binaries used for the run
    ├── cfg/                                           # Configurations for each sequences
    ├── jobs/                                          # Job scripts and logs
    ├── results/octree-hybrid/lossless-geom-lossy-attrs/
    │   ├── bartender_semitracked/                       
    │   │   ├── r01/                               
    │   │   │   ├── frame000.pre.ply.bin.log           # Encoder log
    │   │   │   ├── frame000.pre.ply.bin.time          # Encoding time and RSS
    │   │   │   ├── ...                               
    │   │   │   ├── frame000.pre.ply.bin.decoded.log   # Decoder log
    │   │   │   ├── frame000.pre.ply.bin.decoded.time  # Decoding time and RSS
    │   │   │   └── frame000.pre.ply.bin.decoded.ply   # Decoded PLY
    │   │   ├── ...                                
    │   │   └── r04/
    │   ├── ...                                   
    │   └── ManWithFruit
    └── metrics/
        ├── 1F_TMC13-GS-octree-hybrid-lossless-geom-lossy-attrs-TM13-GS-v3.0-rc2.csv
        ├── 1F_TMC13-GS-octree-predlift-lossless-geom-lossless-attrs-TM13-GS-v3.0-rc2.csv
        │   # Raw results data in CSV format
        └── 1F_TMC13-GS-octree-predlift-lossless-geom-lossless-attrs-TM13-GS-v3.0-rc2.xlsm
            # XLSM comparing Anchor vs Proposed results
```

## Contact
If you have any questions or issues, please contact Jongseok Lee (jongseok1.lee@lge.com).