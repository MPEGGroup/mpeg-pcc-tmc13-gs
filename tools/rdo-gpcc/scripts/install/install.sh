#! /bin/bash
set -uxeo pipefail

# Clean
# rm -rf third_party/*/

# Clone and build tmc13
(
    cd third_party
    git clone -b release-v2.0-rc2 --depth=1 https://git.mpeg.expert/MPEG/3dgh/g-pcc/software/tm/mpeg-pcc-tmc13-gs.git
    cd mpeg-pcc-tmc13-gs
    mkdir -p build
    cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release
    make -j 16
    # Generate config files.
    cd ../cfg
    ../scripts/gen-cfg.sh --gaussiansplat --all
)

# Clone and build mpeg gsc metrics.
(
    cd third_party
    git clone -b release_mpeg154 --depth=1 https://git.mpeg.expert/MPEG/Explorations/GSC/gsc-software/mpeg-gsc-metrics.git
    cd mpeg-gsc-metrics
    bash ./build.sh
)

# Clone and build mpeg gsc tools.
(
    cd third_party
    git clone -b release_mpeg154 --depth=1 https://git.mpeg.expert/MPEG/Explorations/GSC/gsc-software/mpeg-gsc-tools.git
    cd mpeg-gsc-tools
    cd gsTools
    bash ./build.sh
)

# Clone and build mpeg gsc metrics with pruning.
(
    cd third_party
    git clone -b mpeg153-gs-final --depth=1 https://git.mpeg.expert/MPEG/Explorations/GSC/gsc-software/mpeg-gsc-metrics.git mpeg-gsc-metrics-pruning
    cd mpeg-gsc-metrics-pruning
    bash ./build.sh
)

# Clone  pycolmap
(
    cd third_party
    git clone  https://github.com/rmbrualla/pycolmap.git
    cd pycolmap
    git checkout cc7ea4b7301720ac29287dbe450952511b32125e
    git apply ../pycolmap.patch
)

# Install python packages.
# Install fused-ssim seperately because its dependency list is incomplete.
# See https://docs.astral.sh/uv/concepts/projects/config/#build-isolation
uv sync
uv sync --all-extras
