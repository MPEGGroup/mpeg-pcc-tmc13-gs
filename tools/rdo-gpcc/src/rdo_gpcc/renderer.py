from dataclasses import dataclass
from typing import Dict, Tuple

import torch
from gsplat.rendering import rasterization
from torch import Tensor

from rdo_gpcc.splat_model import SplatModel


@dataclass
class RendererConfig:
    pass

class Renderer:
    def __init__(self, config: RendererConfig):
        self.config = config

    def rasterize_splats(
        self,
        splats: SplatModel,
        camtoworlds: Tensor,
        Ks: Tensor,
        width: int,
        height: int,
        **kwargs,
    ) -> Tuple[Tensor, Tensor, Dict]:
        means = splats.means  # [N, 3]
        # quats = F.normalize(self.splats["quats"], dim=-1)  # [N, 4]
        # rasterization does normalization internally
        quats = splats.quats  # [N, 4]
        scales = torch.exp(splats.scales)  # [N, 3]
        opacities = torch.sigmoid(splats.opacities)  # [N,]

        colors = torch.cat([splats.sh0, splats.shN], 1)  # [N, K, 3]

        render_colors, render_alphas, info = rasterization(
            means=means,
            quats=quats,
            scales=scales,
            opacities=opacities,
            colors=colors,
            viewmats=torch.linalg.inv(camtoworlds),  # [C, 4, 4]
            Ks=Ks,  # [C, 3, 3]
            width=width,
            height=height,
            packed=False,
            sh_degree=splats.sh_degree,
            **kwargs,
        )
        return render_colors, render_alphas, info
    