import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd
import torch

from rdo_gpcc.dataset import Dataset
from rdo_gpcc.mpeg.gs_read_write import write3DG_ply
from rdo_gpcc.renderer import Renderer
from rdo_gpcc.splat_model import SplatModel
from rdo_gpcc.util import CONSOLE, run_cmd


def _subsample_camerapos_views(src_ply_path: Path, dst_ply_path: Path, subsample_every: int) -> int:
    """Keep every ``subsample_every``-th camera_position line in a camerapos PLY."""
    sep = b"end_header\n"
    header, _, body = src_ply_path.read_bytes().partition(sep)
    out, n_kept, cam_idx = [], 0, 0
    for line in (header + sep).splitlines(keepends=True):
        # 'comment camera_position ' (trailing space) are view lines; the
        # 'camera_position_header' line has no space and is always kept.
        if line.startswith(b"comment camera_position "):
            if cam_idx % subsample_every == 0:
                out.append(line)
                n_kept += 1
            cam_idx += 1
        else:
            out.append(line)
    dst_ply_path.write_bytes(b"".join(out) + body)
    return n_kept


@dataclass
class PrunerConfig:
    # Pruning ratio (0.0 = no pruning, 1.0 = prune all).
    prune_ratio: float = 0.0
    # Importance score computation mode.
    # "gscm" = use mpeg-gsc-metrics to compute importance scores. CPU-based.
    # "gsplat" = use gsplat rasterizer and autograd to compute importance scores. GPU-based. Faster but slightly worse performance than gscm.
    mode: Literal["gscm", "gsplat"] = "gscm"
    # Path to importance score binary of MPEG GSCM pruning method. Only needed if mode is 'gscm'.
    imp_bin: Path = (
        Path(__file__).parents[2] / "third_party/mpeg-gsc-metrics-pruning/build/Release/bin/importance-scores"
    ).resolve()
    # Path to camera position binary of MPEG GSC tools. Only needed if mode is 'gscm'.
    campos_bin: Path = (
        Path(__file__).parents[2] / "third_party/mpeg-gsc-tools/gsTools/build/msvc/Release/bin/cameraPosition"
    ).resolve()
    # Number of workers for GSCM importance score computation. Only used if mode is 'gscm'.
    # Useful for limiting number of threads to save memory.
    n_gscm_workers: int = 0

    def __post_init__(self):
        # Check that prune_ratio is between 0 and 1.
        if not (0.0 <= self.prune_ratio < 1.0):
            raise ValueError(f"Prune ratio must be between 0 and 1. Got {self.prune_ratio}.")
        # Check that imp_bin and campos_bin exist.
        if self.prune_ratio > 0 and self.mode == "gscm":
            if not self.campos_bin.exists():
                raise FileNotFoundError(f"Camera position binary not found at {self.campos_bin}.")
            if not self.imp_bin.exists():
                raise FileNotFoundError(f"Importance score binary not found at {self.imp_bin}.")


class Pruner:
    def __init__(self, config: PrunerConfig, dataset: Dataset, renderer: Renderer) -> None:
        self.config = config
        self.dataset = dataset
        self.renderer = renderer

    @torch.no_grad()
    def prune(self, splats: SplatModel) -> SplatModel:
        n_before = len(splats.means)

        scores = self.compute_importance_score(splats)
        n_prune = int(n_before * self.config.prune_ratio)
        mask = torch.zeros(n_before, dtype=torch.bool, device=scores.device)
        if n_prune > 0:
            mask[torch.argsort(scores)[:n_prune]] = True  # True = to be removed
        new_splats = self._prune(mask, splats)

        n_after = len(new_splats.means)
        CONSOLE.log(
            f"Pruned {n_before - n_after} / {n_before} Gaussians ({(n_before - n_after) / n_before * 100:.1f} %). Remaining: {n_after}."
        )
        CONSOLE.log(f"Importance scores (sum):\n{pd.Series(scores.cpu().numpy()).describe()}")
        return new_splats

    @torch.no_grad()
    def _prune(self, mask: torch.Tensor, splats: SplatModel) -> SplatModel:
        sel = torch.where(~mask)[0]

        new_splats = SplatModel(
            sh_degree=splats.sh_degree,
            means=splats.means[sel],
            scales=torch.nn.Parameter(splats.scales[sel]),
            quats=torch.nn.Parameter(splats.quats[sel]),
            opacities=torch.nn.Parameter(splats.opacities[sel]),
            sh0=torch.nn.Parameter(splats.sh0[sel]),
            shN=torch.nn.Parameter(splats.shN[sel]),
        )

        return new_splats

    def compute_importance_score(self, splats: SplatModel) -> torch.Tensor:
        if self.config.mode == "gsplat":
            return self._compute_importance_score_gsplat(splats)
        elif self.config.mode == "gscm":
            return self._compute_importance_score_gscm(splats)
        else:
            raise ValueError(f"Unknown prune mode: {self.config.mode}")

    @torch.enable_grad()
    def _compute_importance_score_gsplat(self, splats: SplatModel) -> torch.Tensor:
        """Compute per-Gaussian importance score. Higher = more important."""
        start_time = time.time()
        n_splats = len(splats.means)
        device = splats.means.device
        scores_sum = torch.zeros(n_splats, device=device)

        trainloader = torch.utils.data.DataLoader(
            self.dataset.get_trainset(skip_image=True), batch_size=1, shuffle=False, num_workers=1
        )

        white_colors = torch.nn.parameter.Parameter(torch.ones((n_splats, 1), dtype=torch.float32, device=device))

        for i, data in enumerate(trainloader):
            camtoworlds = data["camtoworld"].to(device)  # [1, 4, 4]
            Ks = data["K"].to(device)  # [1, 3, 3]
            height, width = data["height"], data["width"]

            colors, alphas, info = self.renderer.rasterize_splats(
                splats=splats,
                camtoworlds=camtoworlds,
                Ks=Ks,
                width=width,
                height=height,
                override_colors=white_colors,
            )
            sum_pixel_value = colors.sum()
            scores_per_view = torch.autograd.grad(sum_pixel_value, white_colors)[0]
            assert scores_per_view.shape == (n_splats, 1)
            assert torch.isclose(scores_per_view.sum(), sum_pixel_value, rtol=1e-3)
            with torch.no_grad():
                scores_sum += scores_per_view.squeeze(-1)

        torch.cuda.synchronize()
        CONSOLE.log(f"Computed importance scores using gsplat in {time.time() - start_time:.1f} seconds.")

        return scores_sum

    def _compute_importance_score_gscm(self, splats: SplatModel) -> torch.Tensor:
        """Compute per-Gaussian importance score using GSCM."""

        # Temporary directory for storing ply and importance score file.
        temp_dir = tempfile.TemporaryDirectory()
        input_ply_path = Path(temp_dir.name) / "point_cloud.ply"
        camerapos_ply_path = Path(temp_dir.name) / "point_cloud.camerapos.ply"
        importance_score_path = Path(temp_dir.name) / "importance_scores.txt"
        CONSOLE.log(f"Using temporary directory {temp_dir.name} for GSCM importance score computation.")

        CONSOLE.log("Saving gaussians.")
        pos = splats.means.data.cpu().numpy()
        sh = torch.cat((splats.sh0, splats.shN), dim=1).data.cpu().numpy()
        opacity = splats.opacities.unsqueeze(-1).data.cpu().numpy()
        scale = splats.scales.data.cpu().numpy()
        rot = splats.quats.data.cpu().numpy()
        write3DG_ply(pos, sh, opacity, scale, rot, True, input_ply_path)

        CONSOLE.log("Adding camera parameters to ply file.")
        run_cmd(
            f"{self.config.campos_bin}"
            f" --input={input_ply_path}"
            f" --output={camerapos_ply_path}"
            f" --image={self.dataset.config.dataset_dir / 'sparse' / 'images.bin'}"
            f" --camera={self.dataset.config.dataset_dir / 'sparse' / 'cameras.bin'}"
            " -v 1"
        )

        # cameraPosition embeds ALL cameras, ignoring dataset.subsample_every; match it here.
        subsample_every = self.dataset.config.subsample_every
        if subsample_every > 1:
            n_views = _subsample_camerapos_views(camerapos_ply_path, camerapos_ply_path, subsample_every)
            CONSOLE.log(f"Subsampled pruning views (every {subsample_every}): kept {n_views}.")

        CONSOLE.log("Computing importance scores.")
        camera = self.dataset.get_trainset(skip_image=True)[0]
        width, height = camera["width"], camera["height"]
        assert all(c["width"] == width and c["height"] == height for c in self.dataset.get_trainset(skip_image=True))

        start_time = time.time()
        run_cmd(
            f"{self.config.imp_bin}"
            f" -a {camerapos_ply_path}"
            f" --impScore={importance_score_path}"
            f" --width={width} --height={height}"
            f" -v 1"
            f" --threads={self.config.n_gscm_workers}"
        )
        CONSOLE.log(f"Computed importance scores using GSCM in {time.time() - start_time:.1f} seconds.")

        CONSOLE.log("Reading importance scores.")
        importance_scores = np.loadtxt(importance_score_path)
        importance_scores = torch.from_numpy(importance_scores).to(splats.means)
        assert importance_scores.shape == (len(splats.means),)
        return importance_scores
