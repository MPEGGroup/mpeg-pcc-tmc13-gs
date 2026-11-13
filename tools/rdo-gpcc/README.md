# RDO-GPCC <!-- omit in toc -->
Codebase of  [m73559 [GSC][JEE 6.2-related] Improving G-PCC Anchor with Gaussian Self-Finetuning](https://git.mpeg.expert/MPEG/Explorations/GSC/contributions/-/issues/142).

## Table of Contents <!-- omit in toc -->
- [Requirements](#requirements)
- [Build a Container](#build-a-container)
- [Install Third-party Libraries and Python Dependencies](#install-third-party-libraries-and-python-dependencies)
- [Prepare Datasets](#prepare-datasets)
- [Generate Scripts](#generate-scripts)
- [Run Scripts](#run-scripts)


## Requirements
- [NVIDIA CUDA GPU](https://developer.nvidia.com/cuda-gpus) 
- [Singularity](https://github.com/sylabs/singularity)
- [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html)

Tested environment is
- NVIDIA A100 40GB
- Singularity-CE 3.11.5
- NVIDIA Container Toolkit 1.17.5

Although not tested, the code would work on other GPUs and newer version of Singularity and NVIDIA Container Toolkit.

For newer GPUs with CUDA capabilities > sm_90 (e.g., GeForce RTX 5090), you need to upgrade PyTorch to the latest version.
We have prepared an example Singularity definition file and configuration file in `examples/pt29cu130`.
Please replace `singularity/rdo-gpcc.def` and `pyproject.toml` with those under `examples/pt29cu130`.

If you cannot install Singularity and NVIDIA Container Toolkit, 
you can manually install dependencies such as CUDA, cmake, uv to the host machine.
See `singularity/rdo-gpcc.def` for details.

## Build a Container
(If you have manually installed dependencies to the host machine, skip this section.)

You can build a Singularity container from the Singularity definition file `singularity/rdo-gpcc.def` by `singularity build` command.
```
$ cd singularity/
$ sudo singularity build rdo-gpcc.sif rdo-gpcc.def  
```

You can launch the container by `singularity shell` command.
```
$ singularity shell --nv --nvccli  -H $(realpath ./singularity/home) ./singularity/rdo-gpcc.sif  
Singularity>
```
We use `--nv --nvccli` options to enable GPU access inside the container and `-H` option to use [sanitized home](https://singularityhub.github.io/singularityware-archive/docs/v2-2/faq#solution-1-specify-the-home-to-mount).

All the following commands should be executed inside the container.

## Install Third-party Libraries and Python Dependencies

This project depends on the following third-party libraries. Clone and build each one, following the build instructions in each repository (newer compatible versions would also work).

| Library | Repository | Tested version |
|---|---|---|
| mpeg-pcc-tmc13-gs | https://git.mpeg.expert/MPEG/3dgh/g-pcc/software/tm/mpeg-pcc-tmc13-gs | `release-v4.0-rc1` |
| mpeg-gsc-metrics | https://git.mpeg.expert/MPEG/Explorations/GSC/gsc-software/mpeg-gsc-metrics | `release_mpeg155` |
| mpeg-gsc-tools | https://git.mpeg.expert/MPEG/Explorations/GSC/gsc-software/mpeg-gsc-tools | `release_mpeg155` |
| mpeg-gsc-importance | https://git.mpeg.expert/MPEG/Explorations/GSC/gsc-software/mpeg-gsc-importance | `release_mpeg155` |

Then install pycolmap and the Python dependencies:
```
Singularity> bash scripts/install/install.sh
```
This clones [pycolmap](https://github.com/rmbrualla/pycolmap) (commit `cc7ea4b7301720ac29287dbe450952511b32125e`), applies a compatibility patch, and installs the Python dependencies from `pyproject.toml`.


## Prepare Datasets
Download the CTC datasets from the MPEG content server and preprocess them so that each frame has a Gaussian-splat PLY with embedded camera positions (`*.camerapos.ply`) and its corresponding COLMAP training views (using the `gstools` tool from `mpeg-gsc-tools`).

## Generate Scripts
`gen-script-with-finetuning.sh` is based on the original `gen-script.sh` of `mpeg-pcc-tmc13-gs`, with the fine-tuning step added.

Copy [`scripts/gen-script/gen-script-with-finetuning.sh`](scripts/gen-script/gen-script-with-finetuning.sh) into `mpeg-pcc-tmc13-gs/scripts/gen-script/` (so it can find its helper files), then edit the settings and placeholder paths at the top. Variables marked with **`*`** are the ones added for fine-tuning (they do not exist in the original `gen-script.sh`) and only need to be set when `FLAG_FINETUNE` ≥ 1:

| Variable | Description |
|---|---|
| `DATASET_PATH` | Processed dataset directory |
| `BASE_OUTPUT` | Output directory for the generated jobs and results |
| `TMC13_DIR` | Path to the built `mpeg-pcc-tmc13-gs` |
| `METRIC` | Path to the `mpeg-gsc-metrics` binary |
| `PYTHON` | Python launcher used to run `preprocess_rdo.py` and the summary scripts (default `python3`; inside the container use `uv run --project ${RDOGPCC_DIR}`) |
| `AUTO_RUN` | Auto-run mode after generation (`0`: off, `1`: local, `2`: Slurm) |
| `FLAG_FINETUNE`* | Fine-tuning mode (see the table below) |
| `RDOGPCC_DIR`* | This repository (defaults to `${TMC13_DIR}/tools/rdo-gpcc`) |
| `GSTOOLS_BIN`* | `gstools` binary from `mpeg-gsc-tools` |
| `IMP_BIN`* | `importance-scores` binary from the pruning build of `mpeg-gsc-metrics` |
| `FLAG_MULTILAYER`* | Multi-layer joint fine-tuning of the SH layers (m77600) after the fine-tuning (`0`: off, `1`: on; see [`README.multilayer.md`](README.multilayer.md)) |
| `MULTILAYER_STEPS`* | Training steps of the multi-layer joint fine-tuning (default `10000`) |

`FLAG_FINETUNE` selects the fine-tuning mode:

| Value | Behavior |
|---|---|
| `0` | Off — skip fine-tuning and encode the input PLY directly (same as the original `gen-script.sh`) |
| `1` | Fine-tune with `cfg/default.cfg` |
| `2` | Fine-tune with per-sequence BD-rate configs (`cfg/per-seq-bdrate/<seq>.cfg`) |
| `3` | Fine-tune with per-sequence BD-area configs (`cfg/per-seq-bdarea/<seq>.cfg`) |

Run the copied script:
```
Singularity> bash mpeg-pcc-tmc13-gs/scripts/gen-script/gen-script-with-finetuning.sh
```
The run scripts are generated under `<BASE_OUTPUT>/<JOBNAME>-<timestamp>/jobs/` as `run_all.sh` (local) and `run_all_slurm.sh` (Slurm).

## Run Scripts
This runs fine-tuning, G-PCC encoding/decoding, and metric computation in sequence.
```
Singularity> bash <BASE_OUTPUT>/<JOBNAME>-<timestamp>/jobs/run_all.sh
```

GPU concurrency during fine-tuning introduces variability in pre-processed GS models across runs.
This results in minor differences in bitsize, PSNR, and BD-rate.
In our experiments, RGB-PSNR BD-rate deviations were typically within 1% (e.g., -27.6% to -27.9% on bartender) under consistent conditions.
