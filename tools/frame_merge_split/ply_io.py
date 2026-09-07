from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


PLY_TO_DTYPE = {
    "char": "i1",
    "uchar": "u1",
    "int8": "i1",
    "uint8": "u1",
    "short": "<i2",
    "ushort": "<u2",
    "int16": "<i2",
    "uint16": "<u2",
    "int": "<i4",
    "uint": "<u4",
    "int32": "<i4",
    "uint32": "<u4",
    "float": "<f4",
    "float32": "<f4",
    "double": "<f8",
    "float64": "<f8",
}

DTYPE_TO_PLY = {
    "i1": "char",
    "u1": "uchar",
    "<i2": "short",
    "<u2": "ushort",
    "<i4": "int",
    "<u4": "uint",
    "<f4": "float",
    "<f8": "double",
}


@dataclass
class PlyData:
    path: Path
    data: np.ndarray
    comments: list[str]
    prefix_lines: list[str]
    vertex_count: int


def _parse_header(path: Path) -> tuple[int, int, list[tuple[str, str]], list[str]]:
    offset = 0
    lines: list[str] = []
    with path.open("rb") as fp:
        while True:
            line = fp.readline()
            if not line:
                raise ValueError(f"PLY header has no end_header: {path}")
            offset += len(line)
            text = line.decode("ascii", errors="replace").rstrip("\n")
            lines.append(text)
            if text.strip() == "end_header":
                break

    if not lines or lines[0].strip() != "ply":
        raise ValueError(f"Not a PLY file: {path}")

    fmt = None
    vertex_count = None
    properties: list[tuple[str, str]] = []
    prefix_lines: list[str] = []
    in_vertex = False

    for line in lines[1:]:
        tokens = line.split()
        if not tokens:
            continue
        key = tokens[0]
        if key == "format":
            fmt = tokens[1]
            if fmt != "binary_little_endian":
                raise ValueError(f"Only binary_little_endian PLY is supported: {path} ({fmt})")
        elif key == "element":
            in_vertex = tokens[1] == "vertex"
            if in_vertex:
                vertex_count = int(tokens[2])
            else:
                prefix_lines.append(line)
        elif key == "property" and in_vertex:
            if tokens[1] == "list":
                raise ValueError(f"List properties are not supported: {path}")
            properties.append((tokens[1], tokens[2]))
        elif key == "end_header":
            pass
        else:
            prefix_lines.append(line)

    if fmt is None or vertex_count is None:
        raise ValueError(f"Invalid PLY header: {path}")
    return offset, vertex_count, properties, prefix_lines


def read_ply(path: str | Path, mmap: bool = True) -> PlyData:
    path = Path(path)
    offset, vertex_count, properties, prefix_lines = _parse_header(path)
    dtype = np.dtype([(name, PLY_TO_DTYPE[ptype]) for ptype, name in properties])
    if mmap:
        data = np.memmap(path, mode="r", dtype=dtype, offset=offset, shape=(vertex_count,))
    else:
        with path.open("rb") as fp:
            fp.seek(offset)
            data = np.fromfile(fp, dtype=dtype, count=vertex_count)
    comments = [line for line in prefix_lines if line.startswith("comment")]
    return PlyData(path=path, data=data, comments=comments, prefix_lines=prefix_lines, vertex_count=vertex_count)


def _ply_type(dtype: np.dtype) -> str:
    key = dtype.str
    if key.startswith("|"):
        key = key[1:]
    if key not in DTYPE_TO_PLY:
        raise ValueError(f"Unsupported dtype for PLY output: {dtype}")
    return DTYPE_TO_PLY[key]


def write_ply(
    path: str | Path,
    data: np.ndarray,
    comments: Iterable[str] | None = None,
    extra_lines: Iterable[str] | None = None,
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if data.dtype.names is None:
        raise ValueError("Structured NumPy array required for PLY output")

    comments = list(comments or [])
    extra_lines = list(extra_lines or [])
    with path.open("wb") as fp:
        fp.write(b"ply\n")
        fp.write(b"format binary_little_endian 1.0\n")
        for line in extra_lines:
            if line.startswith("comment"):
                fp.write((line + "\n").encode("ascii", errors="replace"))
        for line in comments:
            if not line.startswith("comment"):
                line = "comment " + line
            fp.write((line + "\n").encode("ascii", errors="replace"))
        fp.write(f"element vertex {len(data)}\n".encode("ascii"))
        for name in data.dtype.names:
            fp.write(f"property {_ply_type(data.dtype.fields[name][0])} {name}\n".encode("ascii"))
        fp.write(b"end_header\n")
        np.asarray(data).tofile(fp)


def drop_field(data: np.ndarray, field: str) -> np.ndarray:
    """Return `data` without column `field` (unchanged if the field is absent)."""
    names = data.dtype.names or ()
    if field not in names:
        return data
    keep = [n for n in names if n != field]
    out = np.empty(len(data), dtype=[(n, data.dtype[n]) for n in keep])
    for n in keep:
        out[n] = data[n]
    return out


def prune_by_score(data: np.ndarray, frac: float, field: str = "score") -> np.ndarray:
    """Drop the least-important `frac` fraction of points (lowest `field` value).

    No-op when frac<=0 or `field` is absent.  Deterministic (stable argsort), and
    the surviving points keep their original order."""
    if frac <= 0.0 or field not in (data.dtype.names or ()):
        return data
    n = len(data)
    drop = int(round(n * float(frac)))
    if drop <= 0:
        return data
    drop = min(drop, n - 1)
    order = np.argsort(np.asarray(data[field]), kind="stable")  # ascending: low = unimportant
    keep = np.sort(order[drop:])                                # keep the highest (1-frac)
    return data[keep]


def concat_same_schema(parts: list[np.ndarray], dtype: np.dtype | None = None) -> np.ndarray:
    parts = [np.asarray(part) for part in parts if len(part) > 0]
    if not parts:
        if dtype is None:
            raise ValueError("Cannot concatenate empty parts without dtype")
        return np.empty(0, dtype=dtype)
    first_dtype = parts[0].dtype
    for part in parts[1:]:
        if part.dtype != first_dtype:
            raise ValueError("All PLY parts must have identical property schema")
    return np.concatenate(parts)

