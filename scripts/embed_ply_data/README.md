# Embed PLY Data

This directory contains utilities for embedding COLMAP camera positions and
per-point importance scores into Gaussian-splat PLY files.

## Requirements

- Bash
- Python 3
- The `gstools` binary that provides the `camera` and `attribute` commands
- A dataset config file containing `ADD_SEQ` entries, for example:

```bash
ADD_SEQ sequence_name 1920 1080
```

The width and height values are accepted by the config format. The batch
wrapper uses the sequence names to locate the data.

## Dataset Layout

By default, each sequence is expected to contain:

```text
<dataset-root>/<sequence>/plys/frame_XXXX.ply
<dataset-root>/<sequence>/colmap_data/frame_XXXX/sparse/cameras.{bin,txt}
<dataset-root>/<sequence>/colmap_data/frame_XXXX/sparse/images.{bin,txt}
```

When both formats are present, the binary model is preferred.

For score embedding, score files are expected at:

```text
<dataset-root>/<sequence>/scores/frame_XXXX.impScores.txt
```

Each score file must contain one score for every point in its corresponding
PLY file.

## Usage

Run the batch wrapper from this directory or provide its path directly:

```bash
./batch_embed_data.sh \
    --dataset-path /path/to/dataset \
    --config /path/to/dataset.config.sh \
    --gstools-bin /path/to/gstools \
    --overwrite
```

Camera-only mode is the default. Outputs are written beside the input PLYs as
`frame_XXXX.camerapos.ply`.

Embed scores and camera positions together:

```bash
./batch_embed_data.sh \
    --dataset-path /path/to/dataset \
    --config /path/to/dataset.config.sh \
    --gstools-bin /path/to/gstools \
    --embed-scores \
    --overwrite
```

Embed scores without camera positions:

```bash
./batch_embed_data.sh \
    --dataset-path /path/to/dataset \
    --config /path/to/dataset.config.sh \
    --gstools-bin /path/to/gstools \
    --scores-only \
    --overwrite
```

The default score-only output suffix is `.with_score.ply`. Use
`--scores-subfolder`, `--scores-suffix`, `--score-name`, and
`--score-percentile` when the score files or attribute settings differ from
the defaults. Use `--no-score-log-scale` to disable the default log transform.

## Options

The complete option list and examples are available from the wrapper itself:

```bash
./batch_embed_data.sh --help
```

Useful overrides include:

- `--plys-subfolder <name>` and `--colmap-subfolder <name>` for non-default
  dataset layouts
- `--frame-count <N>` to process only the first N frames per sequence
- `--camera-output-suffix <suffix>` and
  `--camera-score-output-suffix <suffix>` to choose output names
- `--python <executable>` to select a Python interpreter

The wrapper skips existing outputs unless `--overwrite` is supplied. It exits
with a non-zero status if any sequence fails or required COLMAP/score data is
missing.
