import random
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import networkx as nx
import numpy as np
import scipy
import tyro
from scipy.spatial.transform import Rotation, Slerp

from rdo_gpcc.mpeg.read_write_model import Camera, Image, qvec2rotmat, read_model, rotmat2qvec, write_model
from rdo_gpcc.util import CONSOLE


@dataclass
class Config:
    dataset_dir: Path
    output_dir: Path
    method: Literal["sequential", "mst"] = "mst"
    save_merged_model: bool = False


def main(config: Config):
    random.seed(0)
    np.random.seed(0)

    CONSOLE.log(config)

    CONSOLE.log("Reading COLMAP data.")
    cameras, images, points3D = read_model(path=config.dataset_dir / "sparse")
    CONSOLE.log(f"# of cameras: {len(cameras)}")
    CONSOLE.log(f"# of images: {len(images)}")
    CONSOLE.log(f"# of points3D: {len(points3D)}")

    # Generate validation views.
    new_images = {}
    new_cameras = {}
    new_points3D = {}
    image_ids = list(images.keys())
    if config.method == "sequential":
        for i in image_ids:
            new_image_id = i + 10000
            new_image, new_camera = interpolate_views(images, cameras, i, i + 1, new_image_id)
            new_images[new_image_id] = new_image
            new_cameras[new_image_id] = new_camera
    elif config.method == "mst":
        CONSOLE.log("Computing distance matrix.")
        camera_centers = []
        n_views = len(images)
        for i in image_ids:
            image = images[i]
            R = qvec2rotmat(image.qvec)  # W2C
            camera_center = -R.T @ image.tvec
            camera_centers.append(camera_center)
        dist_mat = scipy.spatial.distance_matrix(camera_centers, camera_centers)
        assert dist_mat.shape == (n_views, n_views)

        CONSOLE.log("Computing minimum spanning tree.")
        dist_graph = nx.from_numpy_array(dist_mat, nodelist=image_ids)
        assert len(dist_graph.edges) == n_views * (n_views - 1) / 2
        mst = nx.minimum_spanning_tree(dist_graph)
        assert len(mst.edges) == n_views - 1
        for i, edge in enumerate(mst.edges):
            new_image_id = i + 10000
            new_image, new_camera = interpolate_views(images, cameras, edge[0], edge[1], new_image_id)
            new_images[new_image_id] = new_image
            new_cameras[new_image_id] = new_camera

    CONSOLE.log("Saving generated COLMAP data.")

    def save_model(*args, path: Path, **kwds):
        path.mkdir(parents=True, exist_ok=True)
        write_model(*args, path=path, **kwds)

    save_model(new_cameras, new_images, new_points3D, path=config.output_dir / "sparse")

    CONSOLE.log("Saving merged COLMAP data for visualization.")
    if config.save_merged_model:
        save_model(
            cameras | new_cameras,
            images | new_images,
            points3D | new_points3D,
            path=config.output_dir / "sparse.merged.d",
            ext=".txt",
        )


def interpolate_views(images, cameras, id1, id2, new_id):
    image1 = images[id1]
    image2 = images[id2]
    camera1 = cameras[image1.camera_id]
    camera2 = cameras[image2.camera_id]
    assert camera1.model == camera2.model and camera1.model in ("SIMPLE_PINHOLE", "PINHOLE")
    assert camera1.width == camera2.width
    assert camera1.height == camera2.height

    # Interpolate pose in world coordinate.
    R1 = qvec2rotmat(image1.qvec)  # W2C
    R2 = qvec2rotmat(image2.qvec)  # W2C
    slerp = Slerp(
        [0, 1],
        Rotation.concatenate([Rotation.from_matrix(R1.T), Rotation.from_matrix(R2.T)]),
    )
    new_R = slerp(0.5).as_matrix().T  # W2C
    new_qvec = rotmat2qvec(new_R)
    # Interpolate camera center.
    camera_center1 = -R1.T @ image1.tvec
    camera_center2 = -R2.T @ image2.tvec
    new_camera_center = (camera_center1 + camera_center2) / 2
    new_tvec = -new_R @ new_camera_center
    # Interpolate camera intrinsics.
    new_camera_params = (camera1.params + camera2.params) / 2
    # Update image name.
    image_name1 = Path(image1.name)
    image_name2 = Path(image2.name)
    new_image_name = f"{image_name1.stem}-{image_name2.stem}{image_name1.suffix}"

    new_image = Image(
        id=new_id,
        qvec=new_qvec,
        tvec=new_tvec,
        camera_id=new_id,
        name=new_image_name,
        xys=[],
        point3D_ids=[],
    )
    new_camera = Camera(
        id=new_id,
        model=camera1.model,
        width=camera1.width,
        height=camera1.height,
        params=new_camera_params,
    )

    return new_image, new_camera


if __name__ == "__main__":
    main(tyro.cli(Config))
