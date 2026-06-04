#!/bin/bash
set -xeuo pipefail

GSCTOOL="third_party/mpeg-gsc-tools"

input_dir=data/CTC
output_dir=data/CTC.processed


datasets=(
    # Forward-facing scenes
    bartender_tracked,0
    bartender_semitracked,0
    cinema_tracked,0
    cinema_semitracked,0
    breakfast_tracked,0
    breakfast_semitracked,0
    breakfast_untracked,0
    breakdance_untracked,0
    trio,0
    makeup,0
    musical,0
    # Object-centric scenes
    manwithfruit_tracked,81
    lego_ferrari,0
    lego_bugatti,0
    cricket_player,0
    plant,0
    solo_tango_female,0
    solo_tango_male,0
    tango_duo,0
    tennis_player,0
    library,0
    flowerdance,0
    gymnast,200
)

for dataset in "${datasets[@]}"; do
    IFS=',' read -r dataset frame <<< "$dataset"
    printf -v frame "%04d" "$frame"

    # Add camera parameters to ply files.
    mkdir -p ${output_dir}/${dataset}/plys
    ${GSCTOOL}/gsTools/build/msvc/Release/bin/cameraPosition \
        --input=${input_dir}/${dataset}/plys/frame_${frame}.ply \
        --camera=${input_dir}/${dataset}/colmap_data/frame_${frame}/sparse/cameras.bin \
        --image=${input_dir}/${dataset}/colmap_data/frame_${frame}/sparse/images.bin \
        --output=${output_dir}/${dataset}/plys/frame_${frame}.camerapos.ply \
        -v 1

    # Copy original views.
    mkdir -p ${output_dir}/${dataset}/colmap_data/frame_${frame}
    cp -r ${input_dir}/${dataset}/colmap_data/frame_${frame}/sparse ${output_dir}/${dataset}/colmap_data/frame_${frame}/
done

