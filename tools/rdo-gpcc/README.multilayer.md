# Multi-layer joint fine-tuning of SH layers (m77600)

Encoder-side tool for CE13.87. It adds a second fine-tuning stage that jointly
optimises all spherical-harmonic layers, so a scalable SH bitstream truncated to
a lower degree still renders well, instead of only the full degree-3 model being
optimal.

Everything below is on top of `tools/rdo-gpcc/README.md`; read that first for the
container-based setup. This file covers only what differs for this contribution
and the pitfalls found while running the full CTC set.

---

## 1. Requirements

* NVIDIA GPU. **Check the compute capability first:**

  ```
  nvidia-smi --query-gpu=name,compute_cap --format=csv
  ```

  `pyproject.toml` pins torch 2.4.1 + cu124, which has no kernels for compute
  capability 12.0. On sm_80/sm_90 (A100 and similar) it works unchanged. On
  compute capability 12.0 the run aborts at the first CUDA kernel; use the
  torch 2.9.1 + cu130 `pyproject.toml` (and Singularity definition) in
  `examples/pt29cu130/`, as the main README describes.

  The results reported for this contribution were produced on compute
  capability 12.0 with torch 2.9.1 + cu130. Fine-tuning is not bit-reproducible
  across different GPU, torch or CUDA versions -- see section 9.
* CUDA toolkit — see below. The container route in the main README avoids this.
* **`mpeg-gsc-metrics` built from branch `m76524_support_flexiable_sh`.** The
  per-SH-degree PSNR in this contribution needs `--shDegreeSource` /
  `--shDegreeDecode`, which only that branch provides. Check with:

  ```
  mpeg-gsc-metrics --help | grep shDegree
  ```

  If nothing prints, the binary is from the wrong branch and
  `run_m77600_ctc.sh` will refuse to start.
* ~60 GB free disk for a full 23-sequence run; ~35 GB for the uv package cache.

## 2. CUDA toolkit for building `fused-ssim` (host install only)

`fused-ssim` and `gsplat` compile CUDA at install time. If you are not using the
Singularity container you need a local toolkit. **The main README does not
mention this, and two of the five components are easy to miss:**

| component | why |
|---|---|
| `cuda_nvcc` | the compiler |
| `cuda_cudart` | runtime headers |
| `cuda_cccl` | thrust/cub headers |
| **`cuda_crt`** | CUDA 13 moved `host_config.h` out of `include/crt/`; without it: `fatal error: crt/host_config.h: No such file or directory` |
| **`libnvvm`** | `cicc` is no longer bundled with nvcc; without it: `nvvm/bin/cicc: not found` |

Fetch them from the CUDA redistributable archive and unpack all five into one
prefix (`--strip-components=1` so they merge), then:

```
ln -sfn lib .cuda13/lib64      # torch's cpp_extension looks for lib64
export CUDA_HOME=$PWD/.cuda13
export PATH="$CUDA_HOME/bin:$PATH"
export TORCH_CUDA_ARCH_LIST="12.0"   # your compute_cap, avoids building every arch
```

## 3. Install

```
cd tools/rdo-gpcc
bash scripts/install/install.sh     # third-party libs + uv sync
```

Verify before going further:

```
.venv/bin/python -c "import torch, gsplat, fused_ssim; \
  print(torch.__version__, torch.cuda.get_arch_list())"
```

Your compute capability must appear in the arch list.

## 4. Data

Datasets go under `data/CTC`, processed views under `data/CTC.processed`, exactly
as the main README describes. The processed tree must be regenerated with the
current `cameraPosition`: the camera header gained `width`/`height` columns in
gsTools 4.0, and older processed data lacks them.

> **`library`**: if your copy has 4000 views, the pipeline aborts with
> `[ERROR] ... please set SAMPLE_VIEWS_FOR_LIBRARY=1`. Either download the
> 100-sampled COLMAP version from the content server, or set that flag, which
> subsamples every 40th view. Use the *same* period for fine-tuning
> (`--dataset.subsample-every 40`) and for the metrics, or the two disagree.

## 5. Running the whole experiment (one command)

`scripts/run_m77600_ctc.sh` runs steps 5 to 7 for every sequence in the CTC test
set: fine-tuning (all three arms), encode/decode at each rate, PSNR at SH decode
degrees 0-3, then the CSVs and the XLSM workbooks.

```
cd tools/rdo-gpcc
TMC13_DIR=/path/to/mpeg-pcc-tmc13-gs \
METRIC_BIN=/path/to/mpeg-gsc-metrics/build/Release/bin/mpeg-gsc-metrics \
DATASET_PATH=$PWD/data/CTC.processed \
bash scripts/run_m77600_ctc.sh
```

`METRIC_BIN` must be the `m76524_support_flexiable_sh` build (see section 1);
the script checks this and stops immediately if it is not.

Check what it will do before committing hours to it:

```
DRYRUN=1 STAGES=A bash scripts/run_m77600_ctc.sh
```

Useful settings (all overridable from the environment):

| variable | default | meaning |
|---|---|---|
| `STAGES` | `ABCD` | A fine-tuning, B codec, C metrics, D results. Any subset. |
| `P_CODEC` / `P_METRIC` | 4 / 8 | parallel jobs for stages B and C |
| `STEPS` | 10000 | steps per fine-tuning pass |
| `DEGREE_WEIGHTS` | `0.05 0.05 0.05 1.0` | SH degree weights in the joint loss |
| `SAMPLE_VIEWS_FOR_LIBRARY` | 0 | set to 1 for the older 4000-view library |
| `FT_DIR` / `OUT_DIR` | `data/CTC.finetuned.m77600`, `output/m77600` | where results go |

It is **resumable**: anything whose output already exists is skipped, so a run
that is interrupted can simply be restarted. Stage A is serial because it needs
the GPU; B and C are parallel. Expect roughly a day for the full set, dominated
by the metrics.

Outputs land in `$OUT_DIR/results`:

* `1F_TMC13-GS-octree-hybrid-lossless-geom-lossy-attrs-<arm>.csv` — one per arm
* `1F_CTC_m77600.xlsm` — Anchor.C1 = no fine-tuning, Proposal.C1 = multi-layer
* `1F_CTC_m77600_vs_finetune.xlsm` — Anchor.C1 = single-stage fine-tuning
* `per_sh_degree.csv` — PSNR at each SH decode degree (the scalability result)
* `rd-curves/rd_<seq>.png` — RD curves per sequence: four panels, one per
  decoded SH degree, each holding the three arms
* `rd-curves/rd_all_sequences_sh3.png` — contact sheet of every sequence at SH3

The rate plotted for degree *d* is the truncated stream
`total_Byte - sum(sh_k_Byte for k > d)`, since a decoder that stops at degree
*d* never receives the higher SH layers.  `plot_rd_curves.py --full-rate`
redraws them against the untruncated `total_Byte` instead.

The rest of this section explains what the script does, for anyone who needs to
run the steps individually.

## 5b. Running the fine-tuning by hand

The contribution is a **two-pass** run. The first pass is the ordinary
single-stage fine-tuning; the second starts from its output and jointly
optimises the SH layers.

```
PY=.venv/bin/python
SEQ=bartender_tracked; FRAME=frame_0000
DS=data/CTC.processed/$SEQ
COMMON="--input-ply-path $DS/plys/$FRAME.camerapos.ply
        --dataset-dir    $DS/colmap_data/$FRAME
        --dataset.subsample-every 1
        --quantizer.position-bitdepth 18
        --trainer.neighbor-loss-weight 5e-10 --trainer.lr-scale 0.1
        --pruner.prune-ratio 0.0 --use-rdo True"

# pass 1 - the shared first stage
$PY scripts/preprocess_rdo.py $COMMON --trainer.max-steps 10000 \
    --output-dir  out/$FRAME.stage0.work \
    --output-ply-path out/$FRAME.stage0.ply

# pass 2 - multi-layer joint fine-tuning, starting from pass 1
$PY scripts/preprocess_rdo.py $COMMON --trainer.max-steps 0 \
    --init-ply-path out/$FRAME.stage0.ply \
    --multilayer.enabled True --multilayer.max-steps 10000 \
    --multilayer.degree-weights 0.05 0.05 0.05 1.0 \
    --output-dir  out/$FRAME.finetune.work \
    --output-ply-path out/$FRAME.ply
```

`--input-ply-path` stays the **original** ply in both passes: it defines the
pseudo-GT views. `--init-ply-path` only changes which model is optimised.

### Why two passes rather than one

Running the multi-layer arm with its own first stage also works
(`--trainer.max-steps 10000` and no `--init-ply-path`), but then each arm draws
an independent sample of a stage that is not bit-reproducible on GPU, and that
noise lands in the comparison between arms. Sharing the first pass reduces it:
in CE13.87 Test 1 (five runs of five sequences on two GPUs), the spread of the
difference between the arms was 1.1-1.4x lower with a shared first pass (1.19x
for the deg3 BD-rate). The gain is modest because the multi-layer arm mostly
varies with its first pass, amplified by the second. Sharing is also cheaper:
the first pass runs once instead of twice.

### Via the per-sequence config files

The options work unchanged in the `key: value` config files under `cfg/`:

```
use-rdo: True
trainer.max-steps: 0
pruner.prune-ratio: 0.0
multilayer.enabled: True
multilayer.max-steps: 10_000
multilayer.degree-weights: 0.05 0.05 0.05 1.0
```

Booleans need an explicit `True`/`False` (the CLI uses `FlagConversionOff`);
a bare `--multilayer.enabled` is rejected.

### Via the gen-script

`scripts/gen-script/gen-script-with-finetuning.sh` runs the same two passes
when `FLAG_MULTILAYER=1` (with `FLAG_FINETUNE` >= 1). Each fine-tuning job then
runs pass 1 with the selected config, keeps its output as `<frame>.pass1.ply`,
and runs pass 2 from it for `MULTILAYER_STEPS` steps without pruning. The
encode and metric jobs use the pass-2 model, so the generated results are those
of the multi-layer arm against the anchor results of the gen-script. For the
paired comparison with single-stage fine-tuning and the per-SH-degree PSNR, use
`run_m77600_ctc.sh` (section 5).

## 6. Metrics per SH layer (by hand)

Per-degree PSNR needs the `m76524_support_flexiable_sh` branch of
`mpeg-gsc-metrics`, which adds `--shDegreeSource` / `--shDegreeDecode`. Confirm
with `mpeg-gsc-metrics --help | grep shDegree`; if nothing prints, your binary is
stale — rebuild it.

```
for DEG in 0 1 2 3; do
  mpeg-gsc-metrics -v 0 --startFrame=0 --frameCount=1 \
    --width=$W --height=$H --cpu=1 --useCameraPosition=1 \
    --shDegreeSource=3 --shDegreeDecode=$DEG \
    --source="$ORIGINAL_PLY" --decode="$DECODED_PLY" --output=deg$DEG.metric
done
```

`--shDegreeSource=3` keeps the reference at full quality while the decoded model
is truncated to `$DEG`.

## 7. Results and Excel (by hand)

```
python scripts/gen-script/gsc_parse_gpcc_results.py -d <rate_dir> -s <seq> -r <rate> -f 1
python scripts/gen-script/gsc_csv2excel.py -f 1 \
    -t scripts/gen-script/templates/WG05N0420_template.xlsm \
    -x results.xlsm --aC1 anchor.csv --pC1 proposal.csv
```

`gsc_parse_gpcc_results.py` globs `*.metric` in the directory and **averages
everything it finds**. If you keep all four per-degree files together it will
silently blend the SH degrees. Point it at a directory holding only the
full-degree metric.

## 8. Known issues

* **`CAMPOS_BIN` and `IMP_BIN` must be set even with pruning disabled.** The
  generator `ENSURE_EXEC`s both whenever fine-tuning is on. Their defaults are a
  Windows path (`gsTools/build/msvc/Release/bin/cameraPosition`, which a Linux
  `build.sh` never produces) and the literal placeholder
  `/path/to/mpeg-gsc-metrics/...`. Override both.
* **Pruning goes in pass 1 only.** Pass 1 may prune; pass 2 then starts from the
  pruned model. Keep `--pruner.prune-ratio 0.0` in pass 2, or it prunes that
  model again.
* **Do not launch training with `uv run` if your `pyproject.toml` differs from
  the installed venv** — it re-resolves and can downgrade torch mid-experiment.
  Call `.venv/bin/python` directly, and put `.venv/bin` on `PATH` so gsplat's JIT
  finds `ninja`.
* **gsplat JIT-compiles its kernels on first use** (~3 min) before training
  starts. This is normal, not a hang.

## 9. What agreement to expect

Fine-tuning is **not** bit-reproducible on GPU: gsplat accumulates gradients with
`gpuAtomicAdd`, so summation order varies with block scheduling and diverges over
10 000 steps. Everything downstream is exact — the anchor arm is deterministic
(`--use-rdo False` never enters the trainer; two independent runs produce the
same md5), tmc3 re-encodes to the same md5, and the metrics are bit-identical
across `--threads` values. So all run-to-run movement comes from stage A.

In CE13.87 Test 1, 5 sequences spanning the range of the effect (`musical`
−16.9 % to `flowerdance` −0.6 %) were run 5 times on each of two GPUs (RTX 5090
and RTX PRO 5000). The worst-case std of the deg3 BD-rate was 0.39 pp and
0.60 pp, and that of ΔPSNR at any SH degree 0.085 dB and 0.078 dB. No
systematic difference between the GPUs was found. The CE report proposes these
acceptable ranges for the difference between two independent results:

| quantity | per sequence | test-set average |
|---|---|---|
| deg3 BD-rate, any pair of arms | ±1.5 pp | ±0.3 pp |
| ΔPSNR, SH degrees 0–2 | ±0.2 dB | ±0.05 dB |
| ΔPSNR, SH degree 3 | ±0.3 dB | ±0.05 dB |

The sign of the gain never changed — even `flowerdance`, whose true effect is
only −0.6 %, stayed negative in all ten runs.

**Report BD-rate at degree 3 only.** Below degree 3 the multi-layer curve sits
entirely above the single-stage one, so the two curves barely share a PSNR
interval and the 4-point cubic Bjøntegaard fit is ill-conditioned — at degree 0
it swung between +37.9 % and +61.2 % across identical runs, and for most
sequences it does not evaluate at all. That is a property of the metric, not
noise in the codec, and repeating the run does not fix it. For degrees 0–2 quote
dB deltas, which reproduce to better than 0.09 dB.

Results from a different GPU, torch, or CUDA version will not match exactly.

### Reproducing the deviation measurement

```
bash tools/rdo-gpcc/scripts/run_m77600_deviation.sh
```

It builds the anchor once, repeats the two fine-tuned arms `REPS` times, and
writes `deviation_raw.csv`, `deviation_report.txt` and `deviation.png` under
`output/m77600-deviation/`. Environment overrides:

| variable | default | meaning |
|---|---|---|
| `REPS` | `4` | independent repeats |
| `SEQS` | 5 cheap sequences | sequences to repeat |
| `DEV_FT` / `DEV_OUT` | `data/CTC.deviation`, `output/m77600-deviation` | where to work |

Expect roughly 14–16 h with the defaults on a single GPU; lower `REPS` or `SEQS`
to fit the time available. The run is resumable, so it can be interrupted and
restarted. To fold a full CTC run in as one more sample, pass its results
directory to the analyser:

```
.venv/bin/python tools/rdo-gpcc/scripts/analyze_deviation.py \
    -d output/m77600-deviation -a output/m77600-deviation/anchor/results \
    -e output/m77600/results
```

`run_m77600_ctc.sh` also accepts `SEQS` (subset of the test set) and `ARMS`
(subset of `anchor finetune finetune-multilayer`), which is what the deviation
driver uses to skip rebuilding the deterministic anchor.
