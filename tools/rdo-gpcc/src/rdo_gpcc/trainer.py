import json
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import List

import imageio
import numpy as np
import torch
import torch.nn.functional as F
import yaml
from fused_ssim import fused_ssim
from rich.progress import Progress
from sklearn.neighbors import NearestNeighbors
from torch.utils.tensorboard import SummaryWriter
from torchmetrics.image import PeakSignalNoiseRatio, StructuralSimilarityIndexMeasure

from rdo_gpcc.dataset import Dataset
from rdo_gpcc.gsplat.utils import (
    set_random_seed,
)
from rdo_gpcc.quantizer import Quantizer
from rdo_gpcc.renderer import Renderer
from rdo_gpcc.splat_model import SplatModel
from rdo_gpcc.util import CONSOLE


@dataclass
class TrainerConfig:
    # Directory to save results
    result_dir: Path
    # Quantization-aware training
    use_qat: bool = True
    # Number of training steps
    max_steps: int = 10_000
    # Steps to evaluate the model
    eval_steps: List[int] = field(default_factory=lambda: [1, 5_000, 10_000])
    # Weight for SSIM loss
    ssim_lambda: float = 0.2
    # Weight for neighbor loss
    neighbor_loss_weight: float = 2e-9
    # Number of neighbors.
    n_neighbors: int = 16
    # Learning rate scale.
    lr_scale: float = 1.0
    # Dump information to tensorboard every this steps
    tb_every: int = 100
    # Save training images to tensorboard
    tb_save_image: bool = False
    # Number of workers for training dataloaders.
    n_train_workers: int = 8
    # Number of workers for validation dataloaders.
    n_val_workers: int = 1


class Trainer:
    """Engine for training and testing."""

    def __init__(
        self,
        splats: SplatModel,
        dataset: Dataset,
        renderer: Renderer,
        quantizer: Quantizer,
        device: torch.device,
        cfg: TrainerConfig,
    ) -> None:
        set_random_seed(42)

        self.cfg = cfg
        self.device = device

        self.splats = splats
        self.dataset = dataset
        self.renderer = renderer
        self.quantizer = quantizer

        # Setup output directories.
        cfg.result_dir.mkdir(exist_ok=True)
        self.stats_dir = cfg.result_dir / "stats"
        self.stats_dir.mkdir(exist_ok=True)
        self.render_dir = cfg.result_dir / "render"
        self.render_dir.mkdir(exist_ok=True)

        # Tensorboard
        self.writer = SummaryWriter(log_dir=cfg.result_dir / "tb")

        # Optimizer.
        CONSOLE.log("Initializing optimizer.")
        self.optimizer = torch.optim.Adam(
            [
                {"params": splats.scales, "lr": 5e-3 * cfg.lr_scale, "name": "scales"},
                {"params": splats.quats, "lr": 1e-3 * cfg.lr_scale, "name": "quats"},
                {"params": splats.opacities, "lr": 5e-2 * cfg.lr_scale, "name": "opacities"},
                {"params": splats.sh0, "lr": 2.5e-3 * cfg.lr_scale, "name": "sh0"},
                {"params": splats.shN, "lr": 2.5e-3 / 20 * cfg.lr_scale, "name": "shN"},
            ],
            eps=1e-15,
        )

        # Losses & Metrics.
        self.ssim = StructuralSimilarityIndexMeasure(data_range=1.0).to(self.device)
        self.psnr = PeakSignalNoiseRatio(data_range=1.0).to(self.device)

        # Build KNN.
        self.knn_graph_dists = None
        self.knn_graph_indices = None
        if cfg.neighbor_loss_weight > 0:
            CONSOLE.log("Computing KNN graph.")
            self.knn_graph_indices, self.knn_graph_dists = self.compute_knn_graph(
                self.splats.means, self.cfg.n_neighbors
            )

    def compute_knn_graph(self, positions, k):
        n_splats = len(positions)

        x = positions.cpu().numpy()
        model = NearestNeighbors(n_neighbors=k + 1, metric="euclidean", n_jobs=32).fit(x)
        dists, indices = model.kneighbors(x)

        dists = torch.tensor(dists[:, 1:], dtype=torch.float32, device=self.device)
        indices = torch.tensor(indices[:, 1:], dtype=torch.long, device=self.device)
        assert dists.shape == (n_splats, k)
        assert indices.shape == (n_splats, k)

        return indices, dists

    def train(self):
        cfg = self.cfg
        device = self.device

        # Dump cfg.
        with (cfg.result_dir / "cfg.yml").open("w") as f:
            yaml.dump(vars(cfg), f)

        trainloader = torch.utils.data.DataLoader(
            self.dataset.get_trainset(),
            batch_size=1,
            shuffle=True,
            num_workers=cfg.n_train_workers,
            persistent_workers=True,
            pin_memory=True,
        )
        trainloader_iter = iter(trainloader)

        # Training loop.
        with Progress(console=CONSOLE) as progress:
            for step in progress.track(range(cfg.max_steps)):
                try:
                    data = next(trainloader_iter)
                except StopIteration:
                    trainloader_iter = iter(trainloader)
                    data = next(trainloader_iter)

                camtoworlds = data["camtoworld"].to(device)  # [1, 4, 4]
                Ks = data["K"].to(device)  # [1, 3, 3]
                pixels = data["image"].to(device).float() / 255  # [1, H, W, 3]

                height, width = pixels.shape[1:3]

                # Forward
                if cfg.use_qat:
                    # Quantize and dequantize splats.
                    splats, t_splats = self.quantize_dequantize_splats(self.splats, train=True)

                else:
                    splats = t_splats = self.splats

                colors, alphas, info = self.renderer.rasterize_splats(
                    splats=splats,
                    camtoworlds=camtoworlds,
                    Ks=Ks,
                    width=width,
                    height=height,
                )

                # loss
                l1loss = F.l1_loss(colors, pixels)
                ssimloss = 1.0 - fused_ssim(
                    colors.permute(0, 3, 1, 2),
                    pixels.permute(0, 3, 1, 2),
                    padding="valid",
                )
                # Neighbor loss
                neighbor_loss = torch.zeros(1).to(l1loss)
                if cfg.neighbor_loss_weight > 0:

                    def neighbor_diff(attrs):
                        targets = attrs.detach()[self.knn_graph_indices].mean(1)
                        l1_loss = F.l1_loss(attrs, targets, reduction="none")
                        return l1_loss.sum()

                    neighbor_loss += neighbor_diff(t_splats.sh0)
                    for shN in t_splats.shN.split(split_size=5, dim=1):
                        neighbor_loss += neighbor_diff(shN)
                    neighbor_loss += neighbor_diff(t_splats.scales)
                    neighbor_loss += neighbor_diff(t_splats.quats)
                    neighbor_loss += neighbor_diff(t_splats.opacities)

                loss = (
                    l1loss * (1.0 - cfg.ssim_lambda)
                    + ssimloss * cfg.ssim_lambda
                    + neighbor_loss * cfg.neighbor_loss_weight
                )
                loss.backward()
                if step % 1000 == 0:
                    CONSOLE.log(f"[Step {step}] loss={loss.item():.3g}")

                # optimize
                self.optimizer.step()
                self.optimizer.zero_grad(set_to_none=True)

                # Logging.
                if cfg.tb_every > 0 and step % cfg.tb_every == 0:
                    mem = torch.cuda.max_memory_allocated() / 1024**3
                    self.writer.add_scalar("train/loss", loss.item(), step)
                    self.writer.add_scalar("train/l1loss", l1loss.item(), step)
                    self.writer.add_scalar("train/ssimloss", ssimloss.item(), step)
                    self.writer.add_scalar("train/neighbor_loss", neighbor_loss.item() / len(t_splats.means), step)
                    self.writer.add_scalar("train/num_GS", len(self.splats.means), step)
                    self.writer.add_scalar("train/mem", mem, step)
                    if cfg.tb_save_image:
                        canvas = torch.cat([pixels, colors], dim=2).detach().cpu().numpy()
                        canvas = canvas.reshape(-1, *canvas.shape[2:])
                        self.writer.add_image("train/render", canvas, step)
                    self.writer.flush()

                # eval the full set
                if step in [i - 1 for i in cfg.eval_steps]:
                    self.eval(step)
                    torch.cuda.empty_cache()

        return self.splats

    @torch.no_grad()
    def eval(self, step: int, stage: str = "val", save_images: bool = False):
        """Entry for evaluation."""
        CONSOLE.log("Running evaluation...")
        device = self.device

        valloader = torch.utils.data.DataLoader(
            self.dataset.get_valset(),
            batch_size=1,
            shuffle=False,
            num_workers=self.cfg.n_val_workers,
            pin_memory=True,
        )
        ellipse_time = 0
        metrics_list = []
        for i, data in enumerate(valloader):
            camtoworlds = data["camtoworld"].to(device)
            Ks = data["K"].to(device)
            pixels = data["image"].to(device).float() / 255  # [1, H, W, 3]
            height, width = pixels.shape[1:3]

            torch.cuda.synchronize()
            tic = time.time()
            # Forward
            if self.cfg.use_qat:
                # Quantize and dequantize splats.
                splats, _ = self.quantize_dequantize_splats(self.splats, train=False)
            else:
                splats = self.splats
            colors, _, _ = self.renderer.rasterize_splats(
                splats=splats,
                camtoworlds=camtoworlds,
                Ks=Ks,
                width=width,
                height=height,
            )  # [1, H, W, 3]
            torch.cuda.synchronize()
            ellipse_time += time.time() - tic

            colors = torch.clamp(colors, 0.0, 1.0)

            # write images
            if save_images:
                for name, img in {"render": colors, "target": pixels}.items():
                    img = (img.squeeze(0).cpu().numpy() * 255).astype(np.uint8)
                    imageio.imwrite(f"{self.render_dir}/{stage}_step{step}_{i:04d}_{name}.png", img)

            pixels_p = pixels.permute(0, 3, 1, 2)  # [1, 3, H, W]
            colors_p = colors.permute(0, 3, 1, 2)  # [1, 3, H, W]
            metrics = {}
            metrics["mse"] = torch.nn.functional.mse_loss(colors_p, pixels_p).item()
            metrics["psnr"] = self.psnr(colors_p, pixels_p).item()
            metrics["ssim"] = self.ssim(colors_p, pixels_p).item()
            metrics["name"] = data["name"][0]
            metrics_list.append(metrics)

        ellipse_time /= len(valloader)

        stats = {}
        for k, v in metrics_list[0].items():
            if isinstance(v, (int, float)):
                stats[k] = sum(metrics[k] for metrics in metrics_list) / len(metrics_list)
        stats.update(
            {
                "ellipse_time": ellipse_time,
                "num_GS": len(self.splats.means),
            }
        )
        stats["metrics"] = metrics_list
        CONSOLE.log(
            f"PSNR: {stats['psnr']:.3f}, SSIM: {stats['ssim']:.4f} "
            f"Time: {stats['ellipse_time']:.3f}s/image "
            f"Number of GS: {stats['num_GS']}"
        )
        # save stats as json
        with (self.stats_dir / f"{stage}_step{step:04d}.json").open("w") as f:
            json.dump(stats, f, indent=4)
        # save stats to tensorboard
        for k, v in stats.items():
            if isinstance(v, (int, float)):
                self.writer.add_scalar(f"{stage}/{k}", v, step)
        self.writer.flush()

    def quantize_dequantize_splats(self, in_splats: SplatModel, train=False):
        """Quantize and dequantize splats."""
        t_splats = self.quantizer._transform(in_splats)
        with torch.no_grad():
            q_splats, metadata = self.quantizer._quantize(t_splats)
            dq_splats = self.quantizer._dequantize(q_splats, metadata)
        # STE
        if train:
            dq_splats = SplatModel(
                sh_degree=self.splats.sh_degree,
                means=dq_splats.means,
                scales=dq_splats.scales.detach() + (t_splats.scales - t_splats.scales.detach()),
                quats=dq_splats.quats.detach() + (t_splats.quats - t_splats.quats.detach()),
                opacities=dq_splats.opacities.detach() + (t_splats.opacities - t_splats.opacities.detach()),
                sh0=dq_splats.sh0.detach() + (t_splats.sh0 - t_splats.sh0.detach()),
                shN=dq_splats.shN.detach() + (t_splats.shN - t_splats.shN.detach()),
            )
        out_splats = self.quantizer._inverse_transform(dq_splats)

        return out_splats, t_splats
