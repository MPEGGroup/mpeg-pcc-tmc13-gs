from dataclasses import dataclass
from typing import Union

from torch import Tensor
from torch.nn import Parameter


@dataclass
class SplatModel:
    sh_degree: int
    means: Union[Tensor, Parameter]
    scales: Union[Tensor, Parameter]
    quats: Union[Tensor, Parameter]
    opacities: Union[Tensor, Parameter]
    sh0: Union[Tensor, Parameter]
    shN: Union[Tensor, Parameter]
