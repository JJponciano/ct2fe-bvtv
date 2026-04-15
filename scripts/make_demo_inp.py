#!/usr/bin/env python3
"""Create a simple Abaqus C3D8 demo mesh inside a NIfTI CT volume."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from map_bvtv_to_fe import read_nifti  # noqa: E402


def world_from_voxel(affine: np.ndarray, voxel: np.ndarray) -> np.ndarray:
    points = np.vstack([voxel.T, np.ones(voxel.shape[0])])
    return (affine @ points)[:3].T


def write_c3d8_grid(
    path: Path,
    affine: np.ndarray,
    min_ijk: np.ndarray,
    max_ijk: np.ndarray,
    counts: tuple[int, int, int],
) -> None:
    nx, ny, nz = counts
    xs = np.linspace(min_ijk[0], max_ijk[0], nx + 1)
    ys = np.linspace(min_ijk[1], max_ijk[1], ny + 1)
    zs = np.linspace(min_ijk[2], max_ijk[2], nz + 1)

    node_ids: dict[tuple[int, int, int], int] = {}
    nodes: list[tuple[int, np.ndarray]] = []
    next_node_id = 1
    for k, z in enumerate(zs):
        for j, y in enumerate(ys):
            for i, x in enumerate(xs):
                node_ids[(i, j, k)] = next_node_id
                coord = world_from_voxel(affine, np.array([[x, y, z]], dtype=float))[0]
                nodes.append((next_node_id, coord))
                next_node_id += 1

    elements: list[tuple[int, tuple[int, ...]]] = []
    next_element_id = 1
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                elements.append(
                    (
                        next_element_id,
                        (
                            node_ids[(i, j, k)],
                            node_ids[(i + 1, j, k)],
                            node_ids[(i + 1, j + 1, k)],
                            node_ids[(i, j + 1, k)],
                            node_ids[(i, j, k + 1)],
                            node_ids[(i + 1, j, k + 1)],
                            node_ids[(i + 1, j + 1, k + 1)],
                            node_ids[(i, j + 1, k + 1)],
                        ),
                    )
                )
                next_element_id += 1

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("*Heading\n")
        handle.write("** Demo regular C3D8 mesh generated from CT volume bounds\n")
        handle.write("*Node\n")
        for node_id, coord in nodes:
            handle.write(f"{node_id}, {coord[0]:.9g}, {coord[1]:.9g}, {coord[2]:.9g}\n")
        handle.write("*Element, type=C3D8\n")
        for element_id, conn in elements:
            handle.write(f"{element_id}, " + ", ".join(str(node_id) for node_id in conn) + "\n")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Create a small demo Abaqus mesh in the high-density CT region.")
    parser.add_argument("--ct", required=True, help="Input CT volume, .nii or .nii.gz")
    parser.add_argument("--out", required=True, help="Output Abaqus .inp path")
    parser.add_argument("--elements", default="8,8,6", help="Element counts as nx,ny,nz")
    parser.add_argument("--anchor-threshold", type=float, default=300.0, help="Voxel threshold used to find dense tissue")
    parser.add_argument("--padding-voxels", type=float, default=4.0, help="Padding around detected dense-tissue bbox")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    counts = tuple(int(part) for part in args.elements.split(","))
    if len(counts) != 3 or any(value <= 0 for value in counts):
        raise ValueError("--elements must be three positive integers, e.g. 8,8,6")

    volume = read_nifti(args.ct)
    mask = volume.data >= args.anchor_threshold
    if not np.any(mask):
        fallback = float(np.percentile(volume.data[np.isfinite(volume.data)], 99.0))
        mask = volume.data >= fallback
        print(f"No voxels above {args.anchor_threshold}; using 99th percentile {fallback:.6g}")

    ijk = np.argwhere(mask)
    min_ijk = ijk.min(axis=0).astype(float) - args.padding_voxels
    max_ijk = ijk.max(axis=0).astype(float) + args.padding_voxels
    min_ijk = np.maximum(min_ijk, 0)
    max_ijk = np.minimum(max_ijk, np.array(volume.data.shape, dtype=float) - 1)
    write_c3d8_grid(Path(args.out), volume.affine, min_ijk, max_ijk, counts)

    print(f"CT shape: {volume.data.shape}")
    print(f"Dense-tissue bbox voxels: min={min_ijk}, max={max_ijk}")
    print(f"Elements: {counts[0] * counts[1] * counts[2]}")
    print(f"Wrote: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
