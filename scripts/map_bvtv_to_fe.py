#!/usr/bin/env python3
"""Map local CT-derived BV/TV values to Abaqus finite elements.

The workflow follows the method described in the paper in this workspace:
segment CT voxels, sample a spherical region around each FE element centroid,
and export element-wise BV/TV plus density-dependent elastic properties.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import struct
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
from skimage.filters import threshold_otsu


NIFTI_DTYPES = {
    2: np.uint8,
    4: np.int16,
    8: np.int32,
    16: np.float32,
    64: np.float64,
    256: np.int8,
    512: np.uint16,
    768: np.uint32,
    1024: np.int64,
    1280: np.uint64,
}


@dataclass(frozen=True)
class Volume:
    data: np.ndarray
    affine: np.ndarray
    spacing: np.ndarray


@dataclass(frozen=True)
class Mesh:
    nodes: dict[int, np.ndarray]
    elements: dict[int, list[int]]


@dataclass(frozen=True)
class MappingResult:
    element_id: int
    centroid_x: float
    centroid_y: float
    centroid_z: float
    sample_voxels: int
    bvtv: float
    youngs_modulus: float


@dataclass(frozen=True)
class PeriImplantResult:
    centerline_start_x: float
    centerline_start_y: float
    centerline_start_z: float
    centerline_end_x: float
    centerline_end_y: float
    centerline_end_z: float
    implant_diameter_mm: float
    shell_thickness_mm: float
    height_mm: float
    sample_voxels: int
    bvtv: float
    threshold: float


def _open_binary(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rb")
    return path.open("rb")


def _read_header_value(fmt: str, data: bytes, offset: int, endian: str):
    size = struct.calcsize(fmt)
    return struct.unpack(endian + fmt, data[offset : offset + size])


def _qform_affine(header: bytes, endian: str) -> np.ndarray:
    pixdim = np.array(_read_header_value("8f", header, 76, endian), dtype=float)
    b, c, d = _read_header_value("3f", header, 256, endian)
    xoff, yoff, zoff = _read_header_value("3f", header, 268, endian)
    qfac = -1.0 if pixdim[0] < 0 else 1.0
    dx, dy, dz = pixdim[1], pixdim[2], pixdim[3] * qfac

    a_sq = 1.0 - (b * b + c * c + d * d)
    a = math.sqrt(max(a_sq, 0.0))
    rotation = np.array(
        [
            [a * a + b * b - c * c - d * d, 2 * b * c - 2 * a * d, 2 * b * d + 2 * a * c],
            [2 * b * c + 2 * a * d, a * a + c * c - b * b - d * d, 2 * c * d - 2 * a * b],
            [2 * b * d - 2 * a * c, 2 * c * d + 2 * a * b, a * a + d * d - c * c - b * b],
        ],
        dtype=float,
    )
    affine = np.eye(4, dtype=float)
    affine[:3, :3] = rotation @ np.diag([dx, dy, dz])
    affine[:3, 3] = [xoff, yoff, zoff]
    return affine


def read_nifti(path: str | Path) -> Volume:
    """Read a NIfTI-1 .nii or .nii.gz volume without nibabel."""
    path = Path(path)
    with _open_binary(path) as handle:
        header = handle.read(348)
        if len(header) != 348:
            raise ValueError(f"{path} is too small to be a NIfTI-1 image")

        sizeof_hdr_le = struct.unpack("<i", header[:4])[0]
        sizeof_hdr_be = struct.unpack(">i", header[:4])[0]
        if sizeof_hdr_le == 348:
            endian = "<"
        elif sizeof_hdr_be == 348:
            endian = ">"
        else:
            raise ValueError(f"{path} does not look like a NIfTI-1 image")

        dims = _read_header_value("8h", header, 40, endian)
        ndim = int(dims[0])
        if ndim < 3:
            raise ValueError(f"{path} must contain at least a 3D volume")
        shape = tuple(int(v) for v in dims[1 : ndim + 1])
        if len(shape) > 3:
            shape = shape[:3] + (int(np.prod(shape[3:])),)

        datatype = _read_header_value("h", header, 70, endian)[0]
        if datatype not in NIFTI_DTYPES:
            raise ValueError(f"Unsupported NIfTI datatype code: {datatype}")
        dtype = np.dtype(NIFTI_DTYPES[datatype]).newbyteorder(endian)

        pixdim = np.array(_read_header_value("8f", header, 76, endian), dtype=float)
        vox_offset = int(round(_read_header_value("f", header, 108, endian)[0]))
        slope = float(_read_header_value("f", header, 112, endian)[0])
        intercept = float(_read_header_value("f", header, 116, endian)[0])
        qform_code = _read_header_value("h", header, 252, endian)[0]
        sform_code = _read_header_value("h", header, 254, endian)[0]

        if vox_offset < 348:
            vox_offset = 348
        handle.seek(vox_offset)
        count = int(np.prod(shape))
        data = np.frombuffer(handle.read(count * dtype.itemsize), dtype=dtype, count=count)
        if data.size != count:
            raise ValueError(f"{path} ended before all image voxels were read")

    data = data.reshape(shape, order="F")
    if data.ndim == 4:
        data = data[..., 0]
    data = data.astype(np.float32, copy=False)
    if slope not in (0.0, 1.0) and not math.isnan(slope):
        data = data * slope
    if intercept != 0.0 and not math.isnan(intercept):
        data = data + intercept

    if sform_code > 0:
        affine = np.eye(4, dtype=float)
        affine[0, :] = _read_header_value("4f", header, 280, endian)
        affine[1, :] = _read_header_value("4f", header, 296, endian)
        affine[2, :] = _read_header_value("4f", header, 312, endian)
    elif qform_code > 0:
        affine = _qform_affine(header, endian)
    else:
        affine = np.eye(4, dtype=float)
        affine[0, 0] = pixdim[1] or 1.0
        affine[1, 1] = pixdim[2] or 1.0
        affine[2, 2] = pixdim[3] or 1.0

    spacing = np.linalg.norm(affine[:3, :3], axis=0)
    return Volume(data=data, affine=affine, spacing=spacing)


def parse_abaqus_inp(path: str | Path) -> Mesh:
    """Parse nodes and elements from a basic Abaqus .inp file."""
    nodes: dict[int, np.ndarray] = {}
    elements: dict[int, list[int]] = {}
    mode: str | None = None

    with Path(path).open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("**"):
                continue
            if line.startswith("*"):
                keyword = line.lower().split(",", 1)[0]
                if keyword == "*node":
                    mode = "node"
                elif keyword == "*element":
                    mode = "element"
                else:
                    mode = None
                continue

            parts = [part.strip() for part in line.split(",") if part.strip()]
            if mode == "node":
                if len(parts) < 4:
                    raise ValueError(f"Malformed node line: {line}")
                nodes[int(parts[0])] = np.array([float(parts[1]), float(parts[2]), float(parts[3])])
            elif mode == "element":
                if len(parts) < 2:
                    raise ValueError(f"Malformed element line: {line}")
                elements[int(parts[0])] = [int(part) for part in parts[1:]]

    if not nodes:
        raise ValueError(f"No *Node block found in {path}")
    if not elements:
        raise ValueError(f"No *Element block found in {path}")
    return Mesh(nodes=nodes, elements=elements)


def rigid_transform_from_landmarks(
    source_points: np.ndarray,
    target_points: np.ndarray,
    allow_scaling: bool = False,
) -> np.ndarray:
    """Return a 4x4 affine that maps source landmarks onto target landmarks."""
    source = np.asarray(source_points, dtype=float)
    target = np.asarray(target_points, dtype=float)
    if source.shape != target.shape:
        raise ValueError("Source and target landmarks must have the same shape")
    if source.ndim != 2 or source.shape[1] != 3:
        raise ValueError("Landmarks must be an Nx3 array")
    if source.shape[0] < 3:
        raise ValueError("At least three non-collinear landmarks are required")

    source_centroid = source.mean(axis=0)
    target_centroid = target.mean(axis=0)
    source_centered = source - source_centroid
    target_centered = target - target_centroid

    if np.linalg.matrix_rank(source_centered) < 2 or np.linalg.matrix_rank(target_centered) < 2:
        raise ValueError("Landmarks are degenerate; provide non-collinear points")

    covariance = source_centered.T @ target_centered
    u, singular_values, vt = np.linalg.svd(covariance)
    rotation = vt.T @ u.T
    if np.linalg.det(rotation) < 0:
        vt[-1, :] *= -1
        rotation = vt.T @ u.T

    scale = 1.0
    if allow_scaling:
        denominator = float(np.sum(source_centered**2))
        if denominator <= 0:
            raise ValueError("Cannot estimate scale from coincident source landmarks")
        scale = float(np.sum(singular_values) / denominator)

    transform = np.eye(4, dtype=float)
    transform[:3, :3] = scale * rotation
    transform[:3, 3] = target_centroid - scale * rotation @ source_centroid
    return transform


def load_registration_transform(
    landmarks_json: str | Path | None = None,
    transform_json: str | Path | None = None,
    allow_scaling: bool = False,
) -> np.ndarray:
    """Load a mesh-to-CT/world transform from landmarks or a 4x4 matrix."""
    if landmarks_json and transform_json:
        raise ValueError("Use either landmarks_json or transform_json, not both")
    if transform_json:
        payload = json.loads(Path(transform_json).read_text(encoding="utf-8"))
        matrix = payload.get("mesh_to_ct_transform", payload)
        transform = np.asarray(matrix, dtype=float)
        if transform.shape != (4, 4):
            raise ValueError("Transform JSON must contain a 4x4 matrix")
        return transform
    if landmarks_json:
        payload = json.loads(Path(landmarks_json).read_text(encoding="utf-8"))
        source = payload.get("mesh") or payload.get("source") or payload.get("source_landmarks")
        target = payload.get("ct") or payload.get("target") or payload.get("target_landmarks")
        if source is None or target is None:
            raise ValueError(
                "Landmark JSON must contain mesh/source landmarks and ct/target landmarks"
            )
        return rigid_transform_from_landmarks(np.asarray(source), np.asarray(target), allow_scaling=allow_scaling)
    return np.eye(4, dtype=float)


def transform_point(point: np.ndarray, affine: np.ndarray) -> np.ndarray:
    homogeneous = np.array([point[0], point[1], point[2], 1.0], dtype=float)
    return (affine @ homogeneous)[:3]


def segment_volume(data: np.ndarray, threshold: str | float, bone_is_high: bool = True) -> tuple[np.ndarray, float]:
    finite = data[np.isfinite(data)]
    if finite.size == 0:
        raise ValueError("Volume does not contain finite voxel values")

    if isinstance(threshold, str) and threshold.lower() == "otsu":
        threshold_value = float(threshold_otsu(finite))
    else:
        threshold_value = float(threshold)

    mask = data >= threshold_value if bone_is_high else data <= threshold_value
    return mask.astype(np.uint8), threshold_value


def element_centroid(mesh: Mesh, node_ids: Iterable[int]) -> np.ndarray:
    coords = []
    for node_id in node_ids:
        try:
            coords.append(mesh.nodes[node_id])
        except KeyError as exc:
            raise ValueError(f"Element references missing node {node_id}") from exc
    return np.mean(np.vstack(coords), axis=0)


def spherical_sample_indices(
    center_world: np.ndarray,
    radius_mm: float,
    affine: np.ndarray,
    inverse_affine: np.ndarray,
    shape: tuple[int, int, int],
    spacing: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    center_voxel = inverse_affine @ np.array([center_world[0], center_world[1], center_world[2], 1.0])
    radius_voxels = int(math.ceil(radius_mm / float(np.min(spacing)))) + 1

    starts = np.maximum(np.floor(center_voxel[:3]).astype(int) - radius_voxels, 0)
    stops = np.minimum(np.floor(center_voxel[:3]).astype(int) + radius_voxels + 1, np.array(shape))
    if np.any(starts >= stops):
        empty = np.array([], dtype=int)
        return empty, empty, empty

    ii, jj, kk = np.meshgrid(
        np.arange(starts[0], stops[0]),
        np.arange(starts[1], stops[1]),
        np.arange(starts[2], stops[2]),
        indexing="ij",
    )
    flat = np.vstack([ii.ravel(), jj.ravel(), kk.ravel(), np.ones(ii.size)])
    world = affine @ flat
    dist2 = np.sum((world[:3].T - center_world) ** 2, axis=1)
    keep = dist2 <= radius_mm * radius_mm
    return ii.ravel()[keep], jj.ravel()[keep], kk.ravel()[keep]


def hollow_cylinder_sample_indices(
    start_world: np.ndarray,
    end_world: np.ndarray,
    inner_radius_mm: float,
    shell_thickness_mm: float,
    affine: np.ndarray,
    inverse_affine: np.ndarray,
    shape: tuple[int, int, int],
    spacing: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    axis = np.asarray(end_world, dtype=float) - np.asarray(start_world, dtype=float)
    height_mm = float(np.linalg.norm(axis))
    if height_mm <= 0:
        raise ValueError("Implant centerline start and end points must be distinct")
    direction = axis / height_mm
    outer_radius_mm = inner_radius_mm + shell_thickness_mm

    endpoints = np.vstack([start_world, end_world])
    endpoint_voxels = (inverse_affine @ np.vstack([endpoints.T, np.ones(2)]))[:3].T
    padding_voxels = int(math.ceil(outer_radius_mm / float(np.min(spacing)))) + 2
    starts = np.maximum(np.floor(endpoint_voxels.min(axis=0)).astype(int) - padding_voxels, 0)
    stops = np.minimum(np.ceil(endpoint_voxels.max(axis=0)).astype(int) + padding_voxels + 1, np.array(shape))

    if np.any(starts >= stops):
        empty = np.array([], dtype=int)
        return empty, empty, empty, height_mm

    ii, jj, kk = np.meshgrid(
        np.arange(starts[0], stops[0]),
        np.arange(starts[1], stops[1]),
        np.arange(starts[2], stops[2]),
        indexing="ij",
    )
    flat = np.vstack([ii.ravel(), jj.ravel(), kk.ravel(), np.ones(ii.size)])
    world = (affine @ flat)[:3].T
    relative = world - start_world
    axial = relative @ direction
    radial_vectors = relative - np.outer(axial, direction)
    radial = np.linalg.norm(radial_vectors, axis=1)
    keep = (
        (axial >= 0.0)
        & (axial <= height_mm)
        & (radial >= inner_radius_mm)
        & (radial <= outer_radius_mm)
    )
    return ii.ravel()[keep], jj.ravel()[keep], kk.ravel()[keep], height_mm


def parse_point(text: str) -> np.ndarray:
    values = [float(part.strip()) for part in text.split(",")]
    if len(values) != 3:
        raise ValueError(f"Expected point as x,y,z, got: {text}")
    return np.array(values, dtype=float)


def compute_peri_implant_bvtv(
    volume: Volume,
    threshold: str | float,
    centerline_start: np.ndarray,
    centerline_end: np.ndarray,
    implant_diameter_mm: float,
    shell_thickness_mm: float = 1.0,
    bone_is_high: bool = True,
) -> PeriImplantResult:
    mask, threshold_value = segment_volume(volume.data, threshold=threshold, bone_is_high=bone_is_high)
    inverse_affine = np.linalg.inv(volume.affine)
    ii, jj, kk, height_mm = hollow_cylinder_sample_indices(
        centerline_start,
        centerline_end,
        inner_radius_mm=implant_diameter_mm / 2.0,
        shell_thickness_mm=shell_thickness_mm,
        affine=volume.affine,
        inverse_affine=inverse_affine,
        shape=mask.shape,
        spacing=volume.spacing,
    )
    bvtv = float(mask[ii, jj, kk].mean()) if ii.size else float("nan")
    return PeriImplantResult(
        centerline_start_x=float(centerline_start[0]),
        centerline_start_y=float(centerline_start[1]),
        centerline_start_z=float(centerline_start[2]),
        centerline_end_x=float(centerline_end[0]),
        centerline_end_y=float(centerline_end[1]),
        centerline_end_z=float(centerline_end[2]),
        implant_diameter_mm=float(implant_diameter_mm),
        shell_thickness_mm=float(shell_thickness_mm),
        height_mm=float(height_mm),
        sample_voxels=int(ii.size),
        bvtv=bvtv,
        threshold=float(threshold_value),
    )


def map_bvtv(
    volume: Volume,
    mesh: Mesh,
    threshold: str | float = "otsu",
    sphere_diameter_mm: float = 1.25,
    e0: float = 8534.64,
    k: float = 1.63,
    nu: float = 0.246,
    bone_is_high: bool = True,
    mesh_to_ct_transform: np.ndarray | None = None,
) -> tuple[list[MappingResult], float]:
    mask, threshold_value = segment_volume(volume.data, threshold=threshold, bone_is_high=bone_is_high)
    inverse_affine = np.linalg.inv(volume.affine)
    radius_mm = sphere_diameter_mm / 2.0
    transform = np.eye(4, dtype=float) if mesh_to_ct_transform is None else mesh_to_ct_transform
    results: list[MappingResult] = []

    for element_id in sorted(mesh.elements):
        mesh_centroid = element_centroid(mesh, mesh.elements[element_id])
        centroid = transform_point(mesh_centroid, transform)
        ii, jj, kk = spherical_sample_indices(
            centroid,
            radius_mm,
            volume.affine,
            inverse_affine,
            mask.shape,
            volume.spacing,
        )
        if ii.size == 0:
            bvtv = float("nan")
        else:
            bvtv = float(mask[ii, jj, kk].mean())
        youngs_modulus = float(e0 * (bvtv**k)) if np.isfinite(bvtv) else float("nan")
        results.append(
            MappingResult(
                element_id=element_id,
                centroid_x=float(centroid[0]),
                centroid_y=float(centroid[1]),
                centroid_z=float(centroid[2]),
                sample_voxels=int(ii.size),
                bvtv=bvtv,
                youngs_modulus=youngs_modulus,
            )
        )
    return results, threshold_value


def write_csv(path: str | Path, results: list[MappingResult], threshold_value: float, nu: float) -> None:
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "element_id",
                "centroid_x_mm",
                "centroid_y_mm",
                "centroid_z_mm",
                "sample_voxels",
                "bvtv",
                "youngs_modulus_mpa",
                "poissons_ratio",
                "threshold",
            ]
        )
        for row in results:
            writer.writerow(
                [
                    row.element_id,
                    f"{row.centroid_x:.9g}",
                    f"{row.centroid_y:.9g}",
                    f"{row.centroid_z:.9g}",
                    row.sample_voxels,
                    f"{row.bvtv:.9g}",
                    f"{row.youngs_modulus:.9g}",
                    f"{nu:.9g}",
                    f"{threshold_value:.9g}",
                ]
            )


def write_abaqus_material_include(
    path: str | Path,
    results: list[MappingResult],
    nu: float,
    min_modulus: float,
    name_prefix: str = "BONE",
) -> None:
    """Write per-element material and section definitions for Abaqus."""
    with Path(path).open("w", encoding="utf-8") as handle:
        handle.write("** Element-wise BV/TV-derived elastic material assignments\n")
        handle.write("** Include after nodes/elements in an Abaqus input file.\n")
        for row in results:
            modulus = row.youngs_modulus
            if not np.isfinite(modulus):
                continue
            modulus = max(float(modulus), min_modulus)
            label = f"{name_prefix}_{row.element_id}"
            handle.write(f"*Elset, elset={label}\n")
            handle.write(f"{row.element_id}\n")
            handle.write(f"*Material, name={label}\n")
            handle.write("*Elastic\n")
            handle.write(f"{modulus:.9g}, {nu:.9g}\n")
            handle.write(f"*Solid Section, elset={label}, material={label}\n")
            handle.write(",\n")


def write_peri_implant_csv(path: str | Path, result: PeriImplantResult) -> None:
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "centerline_start_x_mm",
                "centerline_start_y_mm",
                "centerline_start_z_mm",
                "centerline_end_x_mm",
                "centerline_end_y_mm",
                "centerline_end_z_mm",
                "implant_diameter_mm",
                "shell_thickness_mm",
                "height_mm",
                "sample_voxels",
                "bvtv",
                "threshold",
            ]
        )
        writer.writerow(
            [
                f"{result.centerline_start_x:.9g}",
                f"{result.centerline_start_y:.9g}",
                f"{result.centerline_start_z:.9g}",
                f"{result.centerline_end_x:.9g}",
                f"{result.centerline_end_y:.9g}",
                f"{result.centerline_end_z:.9g}",
                f"{result.implant_diameter_mm:.9g}",
                f"{result.shell_thickness_mm:.9g}",
                f"{result.height_mm:.9g}",
                result.sample_voxels,
                f"{result.bvtv:.9g}",
                f"{result.threshold:.9g}",
            ]
        )


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Map local BV/TV from a CT NIfTI volume to Abaqus finite elements."
    )
    parser.add_argument("--ct", required=True, help="Input CT volume, .nii or .nii.gz")
    parser.add_argument("--mesh", required=True, help="Input Abaqus .inp mesh")
    parser.add_argument("--out-prefix", required=True, help="Output file prefix")
    parser.add_argument("--threshold", default="otsu", help="Segmentation threshold or 'otsu' (default)")
    parser.add_argument("--sphere-diameter-mm", type=float, default=1.25, help="Sampling sphere diameter")
    parser.add_argument("--e0", type=float, default=8534.64, help="E0 in MPa for E = E0 * BVTV^k")
    parser.add_argument("--k", type=float, default=1.63, help="Power-law exponent for E = E0 * BVTV^k")
    parser.add_argument("--nu", type=float, default=0.246, help="Poisson's ratio")
    parser.add_argument("--min-modulus", type=float, default=1.0, help="Floor for Abaqus elastic modulus output")
    parser.add_argument(
        "--landmarks-json",
        help=(
            "Optional registration landmarks JSON. It must contain mesh/source and ct/target Nx3 landmarks. "
            "The estimated transform maps mesh coordinates into CT/world coordinates."
        ),
    )
    parser.add_argument(
        "--transform-json",
        help="Optional JSON 4x4 mesh-to-CT transform matrix, or {'mesh_to_ct_transform': [[...], ...]}.",
    )
    parser.add_argument(
        "--allow-landmark-scaling",
        action="store_true",
        help="Allow uniform scaling when estimating the landmark transform.",
    )
    parser.add_argument(
        "--peri-implant-out",
        help="Optional CSV path for article-style peri-implant hollow-cylinder BV/TV.",
    )
    parser.add_argument(
        "--implant-centerline-start",
        help="Implant/pilot-hole centerline start in CT/world coordinates as x,y,z.",
    )
    parser.add_argument(
        "--implant-centerline-end",
        help="Implant/pilot-hole centerline end in CT/world coordinates as x,y,z.",
    )
    parser.add_argument("--implant-diameter-mm", type=float, help="Implant diameter for hollow-cylinder BV/TV.")
    parser.add_argument(
        "--peri-implant-shell-thickness-mm",
        type=float,
        default=1.0,
        help="Hollow-cylinder shell thickness around implant diameter (article used 1 mm).",
    )
    parser.add_argument(
        "--bone-is-low",
        action="store_true",
        help="Use voxels <= threshold as bone instead of voxels >= threshold",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    threshold: str | float
    threshold = "otsu" if args.threshold.lower() == "otsu" else float(args.threshold)

    volume = read_nifti(args.ct)
    mesh = parse_abaqus_inp(args.mesh)
    mesh_to_ct_transform = load_registration_transform(
        landmarks_json=args.landmarks_json,
        transform_json=args.transform_json,
        allow_scaling=args.allow_landmark_scaling,
    )
    results, threshold_value = map_bvtv(
        volume,
        mesh,
        threshold=threshold,
        sphere_diameter_mm=args.sphere_diameter_mm,
        e0=args.e0,
        k=args.k,
        nu=args.nu,
        bone_is_high=not args.bone_is_low,
        mesh_to_ct_transform=mesh_to_ct_transform,
    )

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    csv_path = out_prefix.with_suffix(".csv")
    inp_path = out_prefix.with_suffix(".materials.inp")
    write_csv(csv_path, results, threshold_value=threshold_value, nu=args.nu)
    write_abaqus_material_include(inp_path, results, nu=args.nu, min_modulus=args.min_modulus)

    if args.peri_implant_out:
        missing = [
            name
            for name, value in [
                ("--implant-centerline-start", args.implant_centerline_start),
                ("--implant-centerline-end", args.implant_centerline_end),
                ("--implant-diameter-mm", args.implant_diameter_mm),
            ]
            if value is None
        ]
        if missing:
            parser.error("--peri-implant-out requires " + ", ".join(missing))
        peri_result = compute_peri_implant_bvtv(
            volume,
            threshold=threshold,
            centerline_start=parse_point(args.implant_centerline_start),
            centerline_end=parse_point(args.implant_centerline_end),
            implant_diameter_mm=args.implant_diameter_mm,
            shell_thickness_mm=args.peri_implant_shell_thickness_mm,
            bone_is_high=not args.bone_is_low,
        )
        write_peri_implant_csv(args.peri_implant_out, peri_result)

    valid = np.array([row.bvtv for row in results if np.isfinite(row.bvtv)], dtype=float)
    print(f"CT shape: {volume.data.shape}")
    print(f"Voxel spacing (mm): {volume.spacing}")
    print(f"Elements mapped: {len(results)}")
    print(f"Threshold used: {threshold_value:.6g}")
    if valid.size:
        print(f"BV/TV min/mean/max: {valid.min():.6g} / {valid.mean():.6g} / {valid.max():.6g}")
    print(f"Wrote: {csv_path}")
    print(f"Wrote: {inp_path}")
    if args.peri_implant_out:
        print(f"Wrote: {args.peri_implant_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
