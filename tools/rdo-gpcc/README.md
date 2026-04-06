# RDO-GPCC <!-- omit in toc -->
Codebase of  [m73559 [GSC][JEE 6.2-related] Improving G-PCC Anchor with Gaussian Self-Finetuning](https://git.mpeg.expert/MPEG/Explorations/GSC/contributions/-/issues/142).

## Table of Contents <!-- omit in toc -->
- [Setup](#setup)
  - [Requirements](#requirements)
  - [Build a Container](#build-a-container)
  - [Install Third-party Libraries and Python Dependencies](#install-third-party-libraries-and-python-dependencies)
- [Reproduce Results for Cross-checking](#reproduce-results-for-cross-checking)
  - [Download Datasets](#download-datasets)
  - [Process Datasets](#process-datasets)
  - [Generate Scripts](#generate-scripts)
  - [Run Scripts](#run-scripts)
- [Contact](#contact)


## Setup

### Requirements
- [NVIDIA CUDA GPU](https://developer.nvidia.com/cuda-gpus) 
- [Singularity](https://github.com/sylabs/singularity)
- [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html)

Tested environment is
- NVIDIA A100 40GB
- Singularity-CE 3.11.5
- NVIDIA Container Toolkit 1.17.5

Although not tested, the code would work on other GPUs and newer version of Singularity and NVIDIA Container Toolkit.

For newer GPUs with CUDA capabilities > sm_90 (e.g., GeForce RTX 5090), you need to upgrade PyTorch to latest version.
We have prepared an example Singularity definition file and configuration file in `examples/pt29cu130`.
Please replace `singularity/rdo-gpcc.def` and `pyproject.toml` with those under `examples/pt29cu130`.

If you cannot install Singularity and NVIDIA Container Toolkit, 
you can mannually install dependencies such as CUDA, cmake, uv to the host machine.
See `singularity/rdo-gpcc.def` for details.

### Build a Container
(If you have mannually installed dependencies to the host machine, skip this section.)

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

### Install Third-party Libraries and Python Dependencies
You can build and install third-party libraries and python dependencies by `scripts/install/install.sh`.
```
Singularity> bash scripts/install/install.sh
```

This script clones and builds the following third-party libraries:
- [mpeg-pcc-tmc13-gs](https://git.mpeg.expert/MPEG/3dgh/g-pcc/software/tm/mpeg-pcc-tmc13-gs) (`third_party/mpeg-pcc-tmc13-gs`)
- [mpeg-gsc-metrics](https://git.mpeg.expert/MPEG/Explorations/GSC/gsc-software/mpeg-gsc-metrics.git) (`third_party/mpeg-gsc-metrics`)
- [mpeg-gsc-tools](https://git.mpeg.expert/MPEG/Explorations/GSC/gsc-software/mpeg-gsc-tools.git) (`third_party/mpeg-gsc-tools`)
- [pycolmap](https://github.com/rmbrualla/pycolmap) (`third_party/pycolmap`)

It applys a patch for pycolmap to fix minor compatibility issues.

If you have already cloned mpeg-pcc-tmc13-gs, mpeg-gsc-metrics and mpeg-gsc-tools elsewhere, you can create symbolic links under `third_party`.

It also installs python dependencies in `pyproject.toml`.



## Reproduce Results for Cross-checking
### Download Datasets
Download CTC datasets from MPEG content server and place them under `data/CTC`.

After downloading:
```
data/CTC
 |- bartender_semitracked
 |- bartender_tracked
 |- breakdance_untracked
 ...
```

### Process Datasets
You can process datasets by `scripts/dataset/process_datasets.sh`.

It will add camera parameters to ply files and generate training views by Minimum Spanning Tree.
```
Singularity> bash scripts/dataset/process_datasets.sh
```
Processed datasets will be saved to `data/CTC.processed`.

### Generate Scripts
You can generate scripts by `scripts/gen-script/gen-script-with-finetuning.sh`.
```
Singularity> bash scripts/gen-script/gen-script-with-finetuning.sh
```
The scripts will be generated under `output/gen-script-output/<JOBNAME>/jobs/run_all.sh`

### Run Scripts
You can run experiments by running generated scripts.
```
Singularity> bash output/gen-script-output/<JOBNAME>/jobs/run_all.sh
```

GPU concurrency during fine-tuning introduces variability in pre-processed GS models across runs.
This results in minor differences in bitsize, PSNR, and BD-rate.
In our experiments, RGB-PSNR BD-rate deviations were typically within 1% (e.g., -27.6% to -27.9% on bartender) under consistent conditions.
