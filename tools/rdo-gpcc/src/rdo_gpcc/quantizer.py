import math
from dataclasses import dataclass
from typing import Literal

import numpy as np
import torch

from rdo_gpcc.splat_model import SplatModel


@dataclass
class QuantizerConfig:
    position_bitdepth: int = 18
    sh_dc_bitdepth: int = 12
    sh_ac_bitdepth: int = 12
    opacity_bitdepth: int = 12
    scale_bitdepth: int = 12
    rotation_bitdepth: int = 12
    colorspace: Literal["bt470", "bt709"] = "bt709"


class Quantizer:
    def __init__(self, config: QuantizerConfig):
        self.config = config

    def quantize(self, in_splats: SplatModel):
        t_splats = self._transform(in_splats)
        q_splats, metadata = self._quantize(t_splats)
        return q_splats, metadata

    def dequantize(self, quant_splats: SplatModel, metadata: dict):
        dq_splats = self._dequantize(quant_splats, metadata)
        it_splats = self._inverse_transform(dq_splats)
        return it_splats

    def _transform(self, in_splats: SplatModel):
        n_splats = in_splats.means.shape[0]

        # Position.
        t_means = in_splats.means
        # SH.
        t_sh0 = self.rgb2yuv(in_splats.sh0, colorspace=self.config.colorspace)
        t_shN = self.rgb2yuv(in_splats.shN, colorspace=self.config.colorspace)
        # Opacity and scale.
        t_opacities = in_splats.opacities
        t_scales = in_splats.scales
        # Rotation.
        norms = torch.linalg.vector_norm(in_splats.quats, dim=-1, keepdim=True)
        t_quats = torch.where(
            norms > 1e-3,
            torch.where((in_splats.quats[:, 0:1] < 0), -in_splats.quats, in_splats.quats) / norms,
            torch.tensor([[1, 0, 0, 0]]).to(in_splats.quats),
        )

        t_splats = SplatModel(
            sh_degree=in_splats.sh_degree,
            means=t_means,
            scales=t_scales,
            quats=t_quats,
            opacities=t_opacities,
            sh0=t_sh0,
            shN=t_shN,
        )

        return t_splats

    def _quantize(self, t_splats: SplatModel):
        cfg = self.config

        q_means, pos_min, pos_max = self.quantize_scalar(t_splats.means, cfg.position_bitdepth, use_floor_ceil=True)
        q_sh0, sh_dc_min, sh_dc_max = self.quantize_scalar(t_splats.sh0, cfg.sh_dc_bitdepth)
        q_shN, sh_ac_min, sh_ac_max = self.quantize_scalar(t_splats.shN, cfg.sh_ac_bitdepth)
        q_opacities, op_min, op_max = self.quantize_scalar(t_splats.opacities.unsqueeze(-1), cfg.opacity_bitdepth)
        q_opacities = q_opacities.squeeze(-1)
        q_scales, scale_min, scale_max = self.quantize_scalar(t_splats.scales, cfg.scale_bitdepth)
        q_quats, rot_min, rot_max = self.quantize_scalar(t_splats.quats, cfg.rotation_bitdepth)

        metadata = {
            "pos_min": pos_min,
            "pos_max": pos_max,
            "sh_dc_min": sh_dc_min,
            "sh_dc_max": sh_dc_max,
            "sh_ac_min": sh_ac_min,
            "sh_ac_max": sh_ac_max,
            "op_min": op_min,
            "op_max": op_max,
            "scale_min": scale_min,
            "scale_max": scale_max,
            "rot_min": rot_min,
            "rot_max": rot_max,
        }

        quant_splats = SplatModel(
            sh_degree=t_splats.sh_degree,
            means=q_means,
            scales=q_scales,
            quats=q_quats,
            opacities=q_opacities,
            sh0=q_sh0,
            shN=q_shN,
        )

        return quant_splats, metadata

    def _dequantize(self, q_splats: SplatModel, metadata: dict):
        cfg = self.config

        dq_means = self.dequantize_scalar(
            q_splats.means, cfg.position_bitdepth, metadata["pos_min"], metadata["pos_max"]
        )
        dq_sh0 = self.dequantize_scalar(q_splats.sh0, cfg.sh_dc_bitdepth, metadata["sh_dc_min"], metadata["sh_dc_max"])
        dq_shN = self.dequantize_scalar(q_splats.shN, cfg.sh_ac_bitdepth, metadata["sh_ac_min"], metadata["sh_ac_max"])
        dq_opacities = self.dequantize_scalar(
            q_splats.opacities.unsqueeze(-1), cfg.opacity_bitdepth, metadata["op_min"], metadata["op_max"]
        ).squeeze(-1)
        dq_scales = self.dequantize_scalar(
            q_splats.scales, cfg.scale_bitdepth, metadata["scale_min"], metadata["scale_max"]
        )
        dq_quats = self.dequantize_scalar(
            q_splats.quats, cfg.rotation_bitdepth, metadata["rot_min"], metadata["rot_max"]
        )

        dq_splats = SplatModel(
            sh_degree=q_splats.sh_degree,
            means=dq_means,
            scales=dq_scales,
            quats=dq_quats,
            opacities=dq_opacities,
            sh0=dq_sh0,
            shN=dq_shN,
        )

        return dq_splats

    def _inverse_transform(self, dq_splats: SplatModel):
        n_splats = dq_splats.means.shape[0]

        # Position.
        it_means = dq_splats.means
        # SH
        it_sh0 = self.yuv2rgb(dq_splats.sh0, colorspace=self.config.colorspace)
        it_shN = self.yuv2rgb(dq_splats.shN, colorspace=self.config.colorspace)
        # Opacity and scale.
        it_opacities = dq_splats.opacities
        it_scales = dq_splats.scales
        # Rotation.
        norms = torch.linalg.vector_norm(dq_splats.quats, dim=-1, keepdim=True)
        it_quats = torch.where(
            norms > 1e-3,
            torch.where((dq_splats.quats[:, 0:1] < 0), -dq_splats.quats, dq_splats.quats) / norms,
            torch.tensor([[1, 0, 0, 0]]).to(dq_splats.quats),
        )

        it_splats = SplatModel(
            sh_degree=dq_splats.sh_degree,
            means=it_means,
            scales=it_scales,
            quats=it_quats,
            opacities=it_opacities,
            sh0=it_sh0,
            shN=it_shN,
        )

        return it_splats

    @staticmethod
    def rgb2yuv(rgb: torch.Tensor, colorspace: Literal["bt470", "bt709"] = "bt709") -> torch.Tensor:
        """RGB → YUV/YCbCr."""
        if colorspace == "bt709":
            mat = torch.tensor(
                [
                    [+0.2126, +0.7152, +0.0722],
                    [-0.1146, -0.3854, +0.5],
                    [+0.5, -0.4542, -0.0458],
                ]
            ).to(rgb)
        else:  # bt470
            mat = torch.tensor(
                [
                    [+0.29900, +0.58700, +0.11400],
                    [-0.14713, -0.28886, +0.43600],
                    [+0.61500, -0.51498, -0.10001],
                ]
            ).to(rgb)
        return rgb @ mat.T

    @staticmethod
    def yuv2rgb(yuv: torch.Tensor, colorspace: Literal["bt470", "bt709"] = "bt709") -> torch.Tensor:
        """YUV/YCbCr → RGB."""
        if colorspace == "bt709":
            mat = torch.tensor(
                [
                    [1, 0, +1.5748],
                    [1, -0.1873, -0.4681],
                    [1, +1.8556, 0],
                ]
            ).to(yuv)
        else:  # bt470
            mat = torch.tensor(
                [
                    [1, 0, +1.13983],
                    [1, -0.39465, -0.58060],
                    [1, +2.03211, 0],
                ]
            ).to(yuv)
        return yuv @ mat.T

    @staticmethod
    def quantize_scalar(x: torch.Tensor, bits: int, use_floor_ceil: bool = False):
        assert torch.isfinite(x).all()

        scale = (1 << bits) - 1

        # Compute min/max.
        min_val = float(x.min().item())
        max_val = float(x.max().item())
        if use_floor_ceil:
            min_val = math.floor(min_val)
            max_val = math.ceil(max_val)
        else:
            min_val = round(min_val)
            max_val = round(max_val)
        # Edge-case guards (mirrors quantize_scalar).
        if min_val == 0 or max_val == 0:
            if min_val == 0:
                min_val = -1
            if max_val == 0:
                max_val = 1
        elif min_val == max_val:
            min_val, max_val = min_val - 1, max_val + 1
        min_val = int(min_val)
        max_val = int(max_val)

        x_norm = (x - min_val) / (max_val - min_val)
        q = torch.round(x_norm * scale)
        q = torch.clamp(q, 0, scale)

        return q, min_val, max_val

    @staticmethod
    def dequantize_scalar(q: torch.Tensor, bitdepth: int, min_val: float, max_val: float):
        scale = (1 << bitdepth) - 1
        x = q / scale * (max_val - min_val) + min_val
        return x
