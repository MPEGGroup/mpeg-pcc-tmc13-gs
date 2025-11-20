# The gen-script package for G-PCC-for-GSC simulations

This guide explains how to prepare data, configure paths, and run scripts for JEE6.2, JEE6.6 and G-PCC CEs experiments using **G-PCC-for-GSC**.

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
│   ├── 1F_GPCC-for-GSC-octree-raht-lossless-geom-lossy-attrs-release-mpeg152.csv
│   └── NF_GPCC-for-GSC-octree-raht-lossless-geom-lossy-attrs-release-mpeg151.csv
│
└── templates/                         # XLSM Templates from WG04N0751 for Reporting
    ├── WG04N0751 template.xlsm
```
> Note: The provided anchor result files (.csv) were generated in the LGE environment. **After completing the anchor cross‑check, please update these CSV files accordingly**.

> Note: [MDS25782_WG04_N00751	
](https://dms.mpeg.expert/doc_end_user/documents/152_Geneva/wg11/MDS25782_WG04_N00751.zip) - Common test conditions for Gaussian splat coding



## 1. Preparation 

### MPEG-GSC-Metrics & MPEG-GSC-Tools Software
Clone the repositories:
```
git clone https://git.mpeg.expert/MPEG/Explorations/GSC/gsc-software/mpeg-gsc-metrics
git clone https://git.mpeg.expert/MPEG/Explorations/GSC/gsc-software/mpeg-gsc-tools
```
Build each project separately.
### PLY Files (with Camera Positions)
Prepare all input PLY files with camera positions (`*.camerapos.ply`) using **gsTools** located in `mpeg-gsc-tools/gsTools`

### Python Package Dependency
Install the `openpyxl` Python package (version 3.1.5) to handle XLSM files:
```
python3 -m pip install openpyxl==3.1.5
```

## 2. Build G-PCC-for-GSC
Clone the repository:
```bash
git clone https://git.mpeg.expert/MPEG/Explorations/GSC/gsc-software/g-pcc-for-gsc
cd g-pcc-for-gsc
git fetch origin tag release_mpeg152-rc1
git checkout tags/release_mpeg152-rc1
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
JOBNAME="G-PCC-for-GSC-mpeg152"
```
### Experiment Settings
```
MODEs=("lossless-geom-lossy-attrs")
TOOLs=("octree-raht") # or octree-predlift
AUTO_RUN=1            # 0:off, 1:run_local, 2:run_slurm
FLAG_1F=1             # 0: NF, 1: 1F
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
TMC13_DIR="/path/to/g-pcc-for-gsc"     # TMC13 directory
TMC13_BIN="$TMC13_DIR/build/tmc3/tmc3" # TMC13 binary
CFG_PATH="$TMC13_DIR/$CFG_DIR"         # Configuration directory
METRIC="/path/to/mpeg-gsc-metrics/build/bin/mpeg-gsc-metrics"  # Metric binary
GSCTOOL="/path/to/mpeg-gsc-tools/" # Metric binary
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
▶ Step 1: GPCC encode/decode in parallel
Submitting 100 jobs to SLURM...
All jobs submitted. Waiting for completion...                              
Progress: 100/100 (100%)
▶ Step 2: Metric calculation in parallel
Submitting 4930 jobs to SLURM...
All jobs submitted. Waiting for completion...                              
Progress: 4930/4930 (100%)
▶ Step 3: Summary generation
✅ All tasks completed
```

## 7. Simulation Results
All simulation outputs are organized under:
```
/path/to/output/
└── G-PCC-for-GSC-mpeg152/
    ├── bin/                                           # Executables or binaries used for the run
    ├── jobs/                                          # Job scripts and logs
    ├── results/octree-raht/lossless-geom-lossy-attrs/
    │   ├── bartender_track/                       
    │   │   ├── r01/                               
    │   │   │   ├── frame000.pre.ply.bin.log           # Encoder log
    │   │   │   ├── frame000.pre.ply.bin.time          # Encoding time and RSS
    │   │   │   ├── ...                               
    │   │   │   ├── frame000.pre.ply.bin.decoded.log   # Decoder log
    │   │   │   ├── frame000.pre.ply.bin.decoded.time  # Decoding time and RSS
    │   │   │   └── frame000.pre.ply.bin.decoded.ply   # Decoded PLY
    │   │   ├── ...                                
    │   │   └── r05/
    │   ├── ...                                   
    │   └── ManWithFruit
    └── metrics/
        ├── 1F_anchor__vs__GPCC-for-GSC-octree-raht-lossless-geom-lossy-attrs-G-PCC-for-GSC-mpeg152.csv
        │   # Raw results data in CSV format
        └── 1F-GPCC-for-GSC-octree-raht-lossless-geom-lossy-attrs-G-PCC-for-GSC-mpeg152.xlsm
            # XLSM comparing Anchor vs Proposed results
```

## Contact
If you have any questions or issues, please contact Jongseok Lee (jongseok1.lee@lge.com).