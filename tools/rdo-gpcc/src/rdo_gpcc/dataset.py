import concurrent
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import imageio
import numpy as np
import torch
from rich.progress import track
from typing_extensions import Literal

from rdo_gpcc.gsplat.datasets.colmap import Parser
from rdo_gpcc.renderer import Renderer
from rdo_gpcc.splat_model import SplatModel
from rdo_gpcc.util import CONSOLE


@dataclass
class DatasetConfig:
    # Path to the Mip-NeRF 360 dataset
    dataset_dir: Path
    # Target image.
    target: Literal["gt", "rendered"] = "rendered"
    # Downsample factor for the dataset
    data_factor: int = 1
    # Subsample images for training and testing.
    subsample_every: int = 1
    # Every N images there is a test image
    test_every: int = -1
    # Test image names.
    test_image_names: Optional[List[str]] = None
    # A global scaler that applies to the scene size related parameters
    global_scale: float = 1.0
    # Cashe rendered images.
    cache_images: bool = True


class Dataset:
    def __init__(
        self,
        config: DatasetConfig,
        target_splats: Optional[SplatModel] = None,
        renderer: Optional[Renderer] = None,
        device: Optional[torch.device] = None,
    ):
        self.config = config
        self.target_splats = target_splats
        self.renderer = renderer
        self.device = device

        self.parser = Parser(
            data_dir=config.dataset_dir,
            factor=config.data_factor,
            test_every=config.test_every,
            skip_images=config.target == "rendered",
            skip_points3D=True,
        )

        # Subsample images.
        image_indices = None
        if config.subsample_every > 1:
            n_orig_images = len(self.parser.image_names)
            image_indices = [i for i in range(n_orig_images) if i % config.subsample_every == 0]
            CONSOLE.log(f"Subsampled cameras {n_orig_images} -> {len(image_indices)}.")

        # Create cache directory.
        self.cache_dir = None
        if config.cache_images:
            self.cache_dir = tempfile.TemporaryDirectory()
            CONSOLE.log(f"Created image cache directory {self.cache_dir.name} .")

        cameras = self.generate_cameras(image_indices=image_indices)

        assert not (config.test_every > 0 and config.test_image_names is not None)
        if config.test_every > 0:
            train_cameras = [c for i, c in enumerate(cameras) if i % self.parser.test_every != 0]
            test_cameras = [c for i, c in enumerate(cameras) if i % self.parser.test_every == 0]
        elif config.test_image_names is not None:
            CONSOLE.log(f"Test images: {config.test_image_names}")
            train_cameras = []
            test_cameras = []
            for c in cameras:
                if c["name"] in config.test_image_names:
                    test_cameras.append(c)
                else:
                    train_cameras.append(c)
            assert set(c["name"] for c in test_cameras) == set(config.test_image_names)
        else:
            # NOTE: Shallow copy.
            train_cameras = test_cameras = cameras
        self.train_cameras = train_cameras
        self.test_cameras = test_cameras

    @property
    def trainset(self) -> List[Dict]:
        return DatasetWrapper(self.train_cameras)

    @property
    def valset(self) -> List[Dict]:
        return DatasetWrapper(self.test_cameras)

    def generate_cameras(self, image_indices=None, n_workers=8):
        if image_indices is None:
            image_indices = range(len(self.parser.image_names))

        executor = None
        futures = None
        if self.config.cache_images:
            executor = concurrent.futures.ProcessPoolExecutor(n_workers)
            futures = set()

        cameras = []
        for index in track(image_indices):
            camera_id = self.parser.camera_ids[index]
            name = Path(self.parser.image_names[index]).stem
            K = torch.from_numpy(self.parser.Ks_dict[camera_id].copy()).float()  # undistorted K
            camtoworlds = torch.from_numpy(self.parser.camtoworlds[index]).float()
            width, height = self.parser.imsize_dict[camera_id]
            params = self.parser.params_dict[camera_id]
            assert len(params) == 0  # Should be undistorted.

            image = None
            if self.config.target == "gt":
                image = imageio.imread(self.parser.image_paths[index])[..., :3]
                image = torch.from_numpy(image).float() / 255
            elif self.config.target == "rendered":
                device = self.device
                colors, _, _ = self.renderer.rasterize_splats(
                    splats=self.target_splats,
                    camtoworlds=camtoworlds.unsqueeze(0).to(device),
                    Ks=K.unsqueeze(0).to(device),
                    width=width,
                    height=height,
                )
                image = torch.clamp(colors, 0.0, 1.0).squeeze(0).cpu()
            else:
                raise NotImplementedError
            assert image.shape == (height, width, 3)

            # Cache image.
            if self.config.cache_images:
                save_path = Path(self.cache_dir.name) / f"{name}.png"
                assert not save_path.exists()
                # Wait until worker gets available.
                while len(futures) >= n_workers:
                    done, futures = concurrent.futures.wait(futures, return_when=concurrent.futures.FIRST_COMPLETED)
                    # Check results.
                    for future in done:
                        future.result()
                assert len(futures) < n_workers
                # Submit job.
                executor.submit(
                    imageio.imwrite,
                    save_path,
                    (image * 255).round().clamp(0, 255).numpy().astype(np.uint8),
                )
                image = save_path

            cameras.append({"K": K, "camtoworld": camtoworlds, "image": image, "name": name})

        if self.config.cache_images:
            for future in concurrent.futures.as_completed(futures):
                future.result()
            executor.shutdown()

        return cameras


class DatasetWrapper:
    def __init__(self, dataset: List[Dict]):
        self.dataset = dataset

    def __len__(self):
        return len(self.dataset)

    def __getitem__(self, idx: int):
        data = self.dataset[idx].copy()

        # Load cached image.
        image = data["image"]
        if isinstance(image, Path):
            image_path = image
            image = imageio.imread(image_path)
            image = torch.from_numpy(image)
            data["image"] = image

        return data
