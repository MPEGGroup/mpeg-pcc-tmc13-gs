import time
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
import torch
import tyro
import yaml
from tqdm import tqdm
from typing_extensions import Literal

from rdo_gpcc.dataset import Dataset, DatasetConfig
from rdo_gpcc.mpeg.gs_read_write import read3DG_ply, write3DG_ply
from rdo_gpcc.pruner import Pruner, PrunerConfig
from rdo_gpcc.quantizer import Quantizer, QuantizerConfig
from rdo_gpcc.renderer import Renderer, RendererConfig
from rdo_gpcc.splat_model import SplatModel
from rdo_gpcc.trainer import Trainer, TrainerConfig
from rdo_gpcc.util import CONSOLE


@dataclass
class Config:
    input_ply_path: Path
    output_dir: Path
    dataset_dir: Path
    use_rdo: bool
    output_ply_path: Path | None = None
    device: Literal["cuda", "cpu"] = "cuda"
    max_memory: int | None = None
    quantizer: QuantizerConfig = field(default_factory=lambda: QuantizerConfig())  # noqa: PLW0108
    pruner: PrunerConfig = field(default_factory=lambda: PrunerConfig())  # noqa: PLW0108
    trainer: TrainerConfig = field(default_factory=lambda: TrainerConfig(result_dir=Path()))
    dataset: DatasetConfig = field(default_factory=lambda: DatasetConfig(dataset_dir=Path()))
    renderer: RendererConfig = field(default_factory=lambda: RendererConfig())  # noqa: PLW0108

    def __post_init__(self):
        self.trainer.result_dir = self.output_dir / "rdo"
        self.dataset.dataset_dir = self.dataset_dir
        if self.output_ply_path is None:
            self.output_ply_path = self.output_dir / "point_cloud.ply"


def finetune_splats(
    pos: np.ndarray,
    sh: np.ndarray,
    opacity: np.ndarray,
    scale: np.ndarray,
    rot: np.ndarray,
    config: Config,
):
    start_time = time.time()

    # Create output directory.
    config.output_dir.mkdir(parents=True, exist_ok=True)

    # Dump the whole config.
    with (config.output_dir / "cfg.yml").open("w") as f:
        yaml.dump(asdict(config), f)

    # Set memory limit.
    if config.max_memory is not None and config.device != "cpu":
        total_memory = torch.cuda.get_device_properties(config.device).total_memory / (1024**3)
        CONSOLE.log(f"Setting memory limit to {config.max_memory} / {total_memory} GB")
        torch.cuda.set_per_process_memory_fraction(config.max_memory / total_memory)

    CONSOLE.log("Removing invalid gaussians.")
    n_orig_splats = len(pos)
    mask = (
        np.isfinite(pos).all(1)
        & np.isfinite(sh).all((1, 2))
        & np.isfinite(opacity).all(1)
        & np.isfinite(scale).all(1)
        & np.isfinite(rot).all(1)
    )
    pos = pos[mask]
    sh = sh[mask]
    opacity = opacity[mask]
    scale = scale[mask]
    rot = rot[mask]
    n_new_splats = len(pos)
    CONSOLE.log(
        f"Removed {n_orig_splats - n_new_splats} gaussians ({(n_orig_splats - n_new_splats) / n_orig_splats * 100:.2f} %)."
    )

    CONSOLE.log("Initializing splats.")

    def make_tensor(arr, make_parameter=False):
        tensor = torch.tensor(arr, device=config.device, dtype=torch.float32)
        if make_parameter:
            tensor = torch.nn.Parameter(tensor)
        return tensor

    # Splats to generate pseudo-gt images.
    target_splats = SplatModel(
        sh_degree=3,
        means=make_tensor(pos),
        scales=make_tensor(scale),
        quats=make_tensor(rot),
        opacities=make_tensor(opacity.squeeze(-1)),
        sh0=make_tensor(sh[:, :1, :]),
        shN=make_tensor(sh[:, 1:, :]),
    )

    # The model being optimised may differ from the pseudo-GT reference above.
    if init is None:
        i_pos, i_sh, i_opacity, i_scale, i_rot = pos, sh, opacity, scale, rot
    else:
        i_pos, i_sh, i_opacity, i_scale, i_rot = init
        if len(i_pos) == n_orig_splats:
            # Same point set as the reference: drop the invalid gaussians too.
            i_pos, i_sh = i_pos[mask], i_sh[mask]
            i_opacity, i_scale, i_rot = i_opacity[mask], i_scale[mask], i_rot[mask]
        # Otherwise the init model comes from an earlier stage, which has already
        # removed the invalid gaussians and may also have pruned some.
        CONSOLE.log(f"Starting optimisation from {config.init_ply_path}.")

    # Splats to finetune.
    # We fix positions and optimize the other attributes.
    splats = SplatModel(
        sh_degree=3,
        means=make_tensor(pos),
        scales=make_tensor(scale, True),
        quats=make_tensor(rot, True),
        opacities=make_tensor(opacity.squeeze(-1), True),
        sh0=make_tensor(sh[:, :1, :], True),
        shN=make_tensor(sh[:, 1:, :], True),
    )

    CONSOLE.log("Initializing quantizer.")
    quantizer = Quantizer(config.quantizer)

    CONSOLE.log("Initializing renderer.")
    renderer = Renderer(config.renderer)

    CONSOLE.log("Initializing dataset.")
    dataset = Dataset(config.dataset, target_splats, renderer, config.device)

    # Run pruning.
    if config.pruner.prune_ratio > 0:
        CONSOLE.log("Initializing pruner.")
        pruner = Pruner(config.pruner, dataset, renderer)
        CONSOLE.log("Pruning gaussians.")
        splats = pruner.prune(splats)

    # Run fine-tuning.
    if config.trainer.max_steps > 0:
        CONSOLE.log("Initializing trainer.")
        trainer = Trainer(splats, dataset, renderer, quantizer, config.device, config.trainer)
        CONSOLE.log("Training gaussians.")
        splats = trainer.train()

    # Convert torch tensors to numpy arrays.
    f_pos = splats.means.data.cpu().numpy()
    f_sh = torch.cat((splats.sh0, splats.shN), dim=1).data.cpu().numpy()
    f_opacity = splats.opacities.unsqueeze(-1).data.cpu().numpy()
    f_scale = splats.scales.data.cpu().numpy()
    f_rot = splats.quats.data.cpu().numpy()

    CONSOLE.log("Printing memory stats.")
    CONSOLE.log(torch.cuda.memory.memory_summary())
    CONSOLE.log(f"Elapsed time (s): {time.time() - start_time}")

    return f_pos, f_sh, f_opacity, f_scale, f_rot


def preprocess(config: Config):

    CONSOLE.log("Loading gaussians.")
    pos, sh, opacity, scale, rot, comments = read3DG_ply(config.input_ply_path, tqdm=tqdm)

    # Run fine-tuning.
    if config.use_rdo:
        CONSOLE.log("Finetuning gaussians.")
        pos, sh, opacity, scale, rot = finetune_splats(pos, sh, opacity, scale, rot, config)

    CONSOLE.log("Saving gaussians.")
    write3DG_ply(pos, sh, opacity, scale, rot, True, config.output_ply_path, tqdm, comments=comments)


if __name__ == "__main__":
    # FlagConversionOff makes boolean fields expect an explicit `True`/`False`
    # value (e.g. `--use-rdo True`) instead of `--flag`/`--no-flag` pairs, so
    # external scripts can assemble `--key value` arguments uniformly.
    config: Config = tyro.cli(tyro.conf.FlagConversionOff[Config])
    CONSOLE.log(config)

    CONSOLE.log("Running preprocess.")
    preprocess(config)
