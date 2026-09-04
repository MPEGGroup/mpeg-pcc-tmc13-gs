#! /bin/bash
set -uxeo pipefail

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
