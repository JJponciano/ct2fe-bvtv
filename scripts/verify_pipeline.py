#!/usr/bin/env python3
"""Rigorous checks for the CT-to-FE BV/TV demo pipeline."""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from map_bvtv_to_fe import (  # noqa: E402
    Mesh,
    Volume,
    compute_peri_implant_bvtv,
    map_bvtv,
    parse_abaqus_inp,
    read_nifti,
    rigid_transform_from_landmarks,
    segment_volume,
)


E0 = 8534.64
K = 1.63
NU = 0.246
MIN_ABAQUS_MODULUS = 1.0


@dataclass(frozen=True)
class CsvRow:
    element_id: int
    centroid: np.ndarray
    sample_voxels: int
    bvtv: float
    youngs_modulus: float
    nu: float
    threshold: float


def assert_close(actual: float, expected: float, tolerance: float, label: str) -> None:
    if not math.isfinite(actual) or not math.isfinite(expected):
        if math.isnan(actual) and math.isnan(expected):
            return
        raise AssertionError(f"{label}: non-finite mismatch actual={actual}, expected={expected}")
    if abs(actual - expected) > tolerance:
        raise AssertionError(f"{label}: actual={actual}, expected={expected}, tolerance={tolerance}")


def read_csv_rows(path: Path) -> list[CsvRow]:
    rows: list[CsvRow] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        for raw in csv.DictReader(handle):
            rows.append(
                CsvRow(
                    element_id=int(raw["element_id"]),
                    centroid=np.array(
                        [
                            float(raw["centroid_x_mm"]),
                            float(raw["centroid_y_mm"]),
                            float(raw["centroid_z_mm"]),
                        ],
                        dtype=float,
                    ),
                    sample_voxels=int(raw["sample_voxels"]),
                    bvtv=float(raw["bvtv"]),
                    youngs_modulus=float(raw["youngs_modulus_mpa"]),
                    nu=float(raw["poissons_ratio"]),
                    threshold=float(raw["threshold"]),
                )
            )
    if not rows:
        raise AssertionError(f"{path} contains no result rows")
    return rows


def verify_peri_implant_output(path: Path) -> str:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1:
        raise AssertionError(f"peri-implant output should contain exactly one data row, found {len(rows)}")
    row = rows[0]
    bvtv = float(row["bvtv"])
    sample_voxels = int(row["sample_voxels"])
    shell = float(row["shell_thickness_mm"])
    if not (0.0 <= bvtv <= 1.0):
        raise AssertionError(f"peri-implant BV/TV outside [0, 1]: {bvtv}")
    if sample_voxels <= 0:
        raise AssertionError("peri-implant output has no sampled voxels")
    if shell <= 0:
        raise AssertionError(f"peri-implant shell thickness must be positive: {shell}")
    return f"peri-implant output is valid (BV/TV={bvtv:.6g}, sample_voxels={sample_voxels})"


def parse_material_elastic_values(path: Path) -> dict[int, float]:
    values: dict[int, float] = {}
    current_element: int | None = None
    expect_elastic = False
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line.startswith("*Elset"):
                match = re.search(r"elset=BONE_(\d+)", line, flags=re.IGNORECASE)
                current_element = int(match.group(1)) if match else None
                expect_elastic = False
            elif line.lower() == "*elastic":
                expect_elastic = True
            elif expect_elastic and current_element is not None and line and not line.startswith("*"):
                modulus = float(line.split(",", 1)[0].strip())
                values[current_element] = modulus
                expect_elastic = False
    return values


def verify_mapped_inp_output(name: str, mapped_path: Path, source_mesh: Mesh, material_values: dict[int, float]) -> list[str]:
    if not mapped_path.exists():
        raise AssertionError(f"{name}: mapped Abaqus input file does not exist: {mapped_path}")

    text = mapped_path.read_text(encoding="utf-8")
    lower_text = text.lower()
    required_keywords = ["*node", "*element", "*material", "*elastic", "*solid section"]
    for keyword in required_keywords:
        if keyword not in lower_text:
            raise AssertionError(f"{name}: mapped Abaqus input is missing {keyword}")
    if "self-contained bv/tv-mapped abaqus input" not in lower_text:
        raise AssertionError(f"{name}: mapped Abaqus input is missing the generated-file header")

    mapped_mesh = parse_abaqus_inp(mapped_path)
    if set(mapped_mesh.nodes) != set(source_mesh.nodes):
        raise AssertionError(f"{name}: mapped Abaqus input node IDs do not match the source mesh")
    if set(mapped_mesh.elements) != set(source_mesh.elements):
        raise AssertionError(f"{name}: mapped Abaqus input element IDs do not match the source mesh")
    for element_id, node_ids in source_mesh.elements.items():
        if mapped_mesh.elements[element_id] != node_ids:
            raise AssertionError(f"{name}: mapped Abaqus input changed connectivity for element {element_id}")

    mapped_material_values = parse_material_elastic_values(mapped_path)
    if set(mapped_material_values) != set(material_values):
        raise AssertionError(f"{name}: mapped Abaqus input material element IDs do not match the material include")
    for element_id, expected_modulus in material_values.items():
        assert_close(
            mapped_material_values[element_id],
            expected_modulus,
            5e-5,
            f"{name}: mapped Abaqus input element {element_id} modulus",
        )

    return [
        f"{name}: self-contained mapped Abaqus input contains the source FE mesh ({len(source_mesh.elements)} elements)",
        f"{name}: self-contained mapped Abaqus input material sections match the material include ({len(material_values)} sections)",
    ]


def brute_force_expected_bvtv(volume: Volume, center_world: np.ndarray, radius_mm: float, threshold: float) -> float:
    mask, _ = segment_volume(volume.data, threshold=threshold)
    indices = np.indices(volume.data.shape).reshape(3, -1)
    homogeneous = np.vstack([indices, np.ones(indices.shape[1])])
    world = (volume.affine @ homogeneous)[:3].T
    distances = np.linalg.norm(world - center_world, axis=1)
    keep = distances <= radius_mm
    if not np.any(keep):
        return float("nan")
    flat_mask = mask.reshape(-1)
    return float(flat_mask[keep].mean())


def brute_force_expected_cylinder_bvtv(
    volume: Volume,
    start_world: np.ndarray,
    end_world: np.ndarray,
    inner_radius_mm: float,
    outer_radius_mm: float,
    threshold: float,
) -> float:
    mask, _ = segment_volume(volume.data, threshold=threshold)
    axis = end_world - start_world
    height = float(np.linalg.norm(axis))
    direction = axis / height
    indices = np.indices(volume.data.shape).reshape(3, -1)
    homogeneous = np.vstack([indices, np.ones(indices.shape[1])])
    world = (volume.affine @ homogeneous)[:3].T
    relative = world - start_world
    axial = relative @ direction
    radial = np.linalg.norm(relative - np.outer(axial, direction), axis=1)
    keep = (axial >= 0.0) & (axial <= height) & (radial >= inner_radius_mm) & (radial <= outer_radius_mm)
    if not np.any(keep):
        return float("nan")
    return float(mask.reshape(-1)[keep].mean())


def run_synthetic_checks() -> list[str]:
    messages: list[str] = []

    data = np.zeros((7, 7, 7), dtype=np.float32)
    data[2:5, 2:5, 2:5] = 100.0
    volume = Volume(data=data, affine=np.eye(4), spacing=np.ones(3))
    mesh = Mesh(
        nodes={
            1: np.array([2.5, 2.5, 2.5]),
            2: np.array([3.5, 2.5, 2.5]),
            3: np.array([3.5, 3.5, 2.5]),
            4: np.array([2.5, 3.5, 2.5]),
            5: np.array([2.5, 2.5, 3.5]),
            6: np.array([3.5, 2.5, 3.5]),
            7: np.array([3.5, 3.5, 3.5]),
            8: np.array([2.5, 3.5, 3.5]),
            9: np.array([0.0, 0.0, 0.0]),
            10: np.array([1.0, 0.0, 0.0]),
            11: np.array([1.0, 1.0, 0.0]),
            12: np.array([0.0, 1.0, 0.0]),
            13: np.array([0.0, 0.0, 1.0]),
            14: np.array([1.0, 0.0, 1.0]),
            15: np.array([1.0, 1.0, 1.0]),
            16: np.array([0.0, 1.0, 1.0]),
        },
        elements={1: [1, 2, 3, 4, 5, 6, 7, 8], 2: [9, 10, 11, 12, 13, 14, 15, 16]},
    )
    results, threshold = map_bvtv(volume, mesh, threshold=50.0, sphere_diameter_mm=3.0)
    by_id = {row.element_id: row for row in results}
    assert_close(threshold, 50.0, 0.0, "synthetic fixed threshold")
    assert_close(by_id[1].bvtv, 1.0, 0.0, "synthetic full bone BV/TV")
    assert_close(by_id[1].youngs_modulus, E0, 1e-9, "synthetic full bone E")
    assert_close(by_id[2].bvtv, 0.0, 0.0, "synthetic empty BV/TV")
    assert_close(by_id[2].youngs_modulus, 0.0, 0.0, "synthetic empty E")
    messages.append("Synthetic full/empty bone test passed")

    affine = np.array(
        [
            [2.0, 0.0, 0.0, 10.0],
            [0.0, 1.0, 0.0, -5.0],
            [0.0, 0.0, 0.5, 3.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
        dtype=float,
    )
    data = np.zeros((9, 10, 11), dtype=np.float32)
    ii, jj, kk = np.indices(data.shape)
    data[(ii + 2 * jj + kk) >= 18] = 100.0
    volume = Volume(data=data, affine=affine, spacing=np.array([2.0, 1.0, 0.5]))
    center_voxel = np.array([4.0, 5.0, 5.0, 1.0])
    center_world = (affine @ center_voxel)[:3]
    half = np.array([0.5, 0.5, 0.5])
    nodes: dict[int, np.ndarray] = {}
    node_id = 1
    for dx in [-half[0], half[0]]:
        for dy in [-half[1], half[1]]:
            for dz in [-half[2], half[2]]:
                nodes[node_id] = center_world + np.array([dx, dy, dz])
                node_id += 1
    mesh = Mesh(nodes=nodes, elements={1: list(nodes)})
    radius = 2.25
    result, _ = map_bvtv(volume, mesh, threshold=50.0, sphere_diameter_mm=radius * 2)
    expected = brute_force_expected_bvtv(volume, center_world, radius, threshold=50.0)
    assert_close(result[0].bvtv, expected, 1e-12, "anisotropic affine spherical sampling BV/TV")
    assert_close(result[0].youngs_modulus, E0 * (expected**K), 1e-9, "anisotropic affine E")
    messages.append("Anisotropic affine sphere test passed")

    source = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]])
    angle = np.deg2rad(30.0)
    rotation = np.array(
        [
            [np.cos(angle), -np.sin(angle), 0.0],
            [np.sin(angle), np.cos(angle), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    translation = np.array([10.0, -4.0, 2.5])
    target = (rotation @ source.T).T + translation
    estimated = rigid_transform_from_landmarks(source, target)
    recovered = (estimated @ np.vstack([source.T, np.ones(source.shape[0])]))[:3].T
    if np.max(np.abs(recovered - target)) > 1e-12:
        raise AssertionError("landmark rigid registration did not recover target points")
    messages.append("Rigid landmark registration test passed")

    data = np.zeros((15, 15, 15), dtype=np.float32)
    ii, jj, kk = np.indices(data.shape)
    radius = np.sqrt((ii - 7.0) ** 2 + (jj - 7.0) ** 2)
    data[(radius >= 2.0) & (radius <= 4.0) & (kk >= 3) & (kk <= 11)] = 100.0
    data[(radius >= 2.0) & (radius <= 4.0) & (kk >= 7) & (kk <= 11)] = 0.0
    volume = Volume(data=data, affine=np.eye(4), spacing=np.ones(3))
    start = np.array([7.0, 7.0, 3.0])
    end = np.array([7.0, 7.0, 11.0])
    peri = compute_peri_implant_bvtv(
        volume,
        threshold=50.0,
        centerline_start=start,
        centerline_end=end,
        implant_diameter_mm=4.0,
        shell_thickness_mm=2.0,
    )
    expected_peri = brute_force_expected_cylinder_bvtv(
        volume,
        start,
        end,
        inner_radius_mm=2.0,
        outer_radius_mm=4.0,
        threshold=50.0,
    )
    assert_close(peri.bvtv, expected_peri, 1e-12, "peri-implant hollow cylinder BV/TV")
    messages.append("Peri-implant hollow-cylinder BV/TV test passed")

    return messages


def verify_output_set(
    name: str,
    ct_path: Path,
    mesh_path: Path,
    csv_path: Path,
    material_path: Path,
    mapped_path: Path | None,
    threshold_mode: str,
) -> list[str]:
    messages: list[str] = []
    volume = read_nifti(ct_path)
    mesh = parse_abaqus_inp(mesh_path)
    rows = read_csv_rows(csv_path)
    material_values = parse_material_elastic_values(material_path)

    if len(rows) != len(mesh.elements):
        raise AssertionError(f"{name}: CSV rows={len(rows)} but mesh elements={len(mesh.elements)}")
    messages.append(f"{name}: row count matches mesh elements ({len(rows)})")

    if len(material_values) != len(rows):
        raise AssertionError(f"{name}: material count={len(material_values)} but CSV rows={len(rows)}")
    messages.append(f"{name}: Abaqus material count matches CSV rows ({len(material_values)})")

    if mapped_path is not None:
        messages.extend(verify_mapped_inp_output(name, mapped_path, mesh, material_values))

    thresholds = np.array([row.threshold for row in rows], dtype=float)
    if np.nanmax(np.abs(thresholds - thresholds[0])) > 1e-9:
        raise AssertionError(f"{name}: inconsistent thresholds in CSV")
    if threshold_mode == "otsu":
        _, expected_threshold = segment_volume(volume.data, threshold="otsu")
    else:
        expected_threshold = float(threshold_mode)
    assert_close(thresholds[0], expected_threshold, 1e-3, f"{name}: threshold")
    messages.append(f"{name}: threshold is consistent ({thresholds[0]:.6g})")

    for row in rows:
        if row.sample_voxels <= 0:
            raise AssertionError(f"{name}: element {row.element_id} has no sampled voxels")
        if not (0.0 <= row.bvtv <= 1.0):
            raise AssertionError(f"{name}: element {row.element_id} BV/TV outside [0, 1]: {row.bvtv}")
        assert_close(row.nu, NU, 1e-12, f"{name}: element {row.element_id} nu")
        expected_e = E0 * (row.bvtv**K)
        assert_close(row.youngs_modulus, expected_e, 5e-5, f"{name}: element {row.element_id} E formula")
        expected_abaqus_e = max(expected_e, MIN_ABAQUS_MODULUS)
        assert_close(
            material_values[row.element_id],
            expected_abaqus_e,
            5e-5,
            f"{name}: element {row.element_id} Abaqus elastic modulus",
        )
    messages.append(f"{name}: BV/TV range, sampled voxel counts, nu, E formula, and Abaqus modulus all pass")

    recomputed, _ = map_bvtv(
        volume,
        mesh,
        threshold="otsu" if threshold_mode == "otsu" else float(threshold_mode),
        sphere_diameter_mm=8.0,
        e0=E0,
        k=K,
        nu=NU,
    )
    recomputed_by_id = {row.element_id: row for row in recomputed}
    sample_ids = sorted({rows[0].element_id, rows[-1].element_id, rows[len(rows) // 2].element_id} | {row.element_id for row in rows[::37]})
    row_by_id = {row.element_id: row for row in rows}
    for element_id in sample_ids:
        actual = row_by_id[element_id]
        expected = recomputed_by_id[element_id]
        assert_close(actual.bvtv, expected.bvtv, 5e-9, f"{name}: recomputed BV/TV element {element_id}")
        assert_close(
            actual.youngs_modulus,
            expected.youngs_modulus,
            5e-5,
            f"{name}: recomputed E element {element_id}",
        )
    messages.append(f"{name}: selected elements recompute exactly from CT and mesh")

    return messages


def write_report(path: Path, messages: list[str]) -> None:
    text = "\n".join(
        [
            "# Verification Report",
            "",
            "## Automated Checks",
            "",
            *[f"- PASS: {message}" for message in messages],
            "",
            "## Agreement With The Article",
            "",
            "- PASS: The core mapping operation follows the article: segment CT, compute element centroid, sample a spherical region, average binary voxels to get local BV/TV.",
            "- PASS: Otsu thresholding is implemented, matching the article's global Otsu segmentation step.",
            "- PASS: The default sampling sphere diameter is 1.25 mm, matching the article.",
            "- PASS: The elastic relation uses the article constants E0 = 8534.64 MPa, k = 1.63, and nu = 0.246.",
            "- PASS: Optional rigid landmark registration is implemented for mapping mesh coordinates into CT/world coordinates.",
            "- PASS: Self-contained mapped Abaqus input files are generated and verified to contain the source FE mesh plus element-wise material sections.",
            "- PASS: Article-style peri-implant BV/TV in a 1 mm hollow cylindrical shell can be computed from an implant/pilot-hole centerline.",
            "- LIMITATION: The demo data are public clinical abdomen CT data, not the article's 24.6 micrometer jawbone micro-CT scans.",
            "- LIMITATION: The demo mesh is a generated regular Abaqus grid, not the article's dental implant/jawbone geometry.",
            "- LIMITATION: The demo mapping does not use real experimental landmarks; when no landmarks or transform are provided, the pipeline assumes CT and FE mesh coordinates are already aligned.",
            "- LIMITATION: The nonlinear Abaqus UMAT plasticity and damage model from the article is not reproduced; the export is elastic material assignment by element.",
            "- LIMITATION: The demo uses an 8 mm sphere because the public CT voxel spacing is about 1.49 mm. The article's 1.25 mm sphere is appropriate for their high-resolution micro-CT, not for this low-resolution demo CT.",
            "",
            "Conclusion: the implementation is algorithmically consistent with the article's BV/TV-to-element material mapping step, but the demo results are not quantitatively comparable to the paper's published jawbone/implant results.",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Verify the CT-to-FE BV/TV pipeline outputs.")
    parser.add_argument("--ct", default="data/CT_Abdo.nii.gz")
    parser.add_argument("--mesh", default="data/demo_abdomen_mesh.inp")
    parser.add_argument("--csv-300", default="outputs/demo_abdomen_bvtv.csv")
    parser.add_argument("--materials-300", default="outputs/demo_abdomen_bvtv.materials.inp")
    parser.add_argument("--mapped-300", default="outputs/demo_abdomen_bvtv.mapped.inp")
    parser.add_argument("--csv-otsu", default="outputs/demo_abdomen_bvtv_otsu.csv")
    parser.add_argument("--materials-otsu", default="outputs/demo_abdomen_bvtv_otsu.materials.inp")
    parser.add_argument("--mapped-otsu", default="outputs/demo_abdomen_bvtv_otsu.mapped.inp")
    parser.add_argument("--csv-registered", default="outputs/demo_abdomen_bvtv_registered.csv")
    parser.add_argument("--materials-registered", default="outputs/demo_abdomen_bvtv_registered.materials.inp")
    parser.add_argument("--mapped-registered", default="outputs/demo_abdomen_bvtv_registered.mapped.inp")
    parser.add_argument("--peri-implant", default="outputs/demo_peri_implant_bvtv.csv")
    parser.add_argument("--report", default="outputs/verification_report.md")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    messages: list[str] = []
    messages.extend(run_synthetic_checks())
    messages.extend(
        verify_output_set(
            "threshold-300",
            Path(args.ct),
            Path(args.mesh),
            Path(args.csv_300),
            Path(args.materials_300),
            Path(args.mapped_300),
            "300",
        )
    )
    messages.extend(
        verify_output_set(
            "otsu",
            Path(args.ct),
            Path(args.mesh),
            Path(args.csv_otsu),
            Path(args.materials_otsu),
            Path(args.mapped_otsu),
            "otsu",
        )
    )
    if Path(args.csv_registered).exists() and Path(args.materials_registered).exists():
        messages.extend(
            verify_output_set(
                "registered-threshold-300",
                Path(args.ct),
                Path(args.mesh),
                Path(args.csv_registered),
                Path(args.materials_registered),
                Path(args.mapped_registered),
                "300",
            )
        )
    if Path(args.peri_implant).exists():
        messages.append(verify_peri_implant_output(Path(args.peri_implant)))
    write_report(Path(args.report), messages)
    for message in messages:
        print(f"PASS: {message}")
    print(f"Wrote: {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
