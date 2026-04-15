#!/usr/bin/env python3
"""Generate an HTML visual report for CT-to-FE BV/TV mapping results."""

from __future__ import annotations

import argparse
import csv
import html
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from map_bvtv_to_fe import parse_abaqus_inp, read_nifti, segment_volume  # noqa: E402


@dataclass(frozen=True)
class ResultRow:
    element_id: int
    centroid: np.ndarray
    sample_voxels: int
    bvtv: float
    youngs_modulus: float
    threshold: float


def read_results(path: str | Path) -> list[ResultRow]:
    rows: list[ResultRow] = []
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rows.append(
                ResultRow(
                    element_id=int(row["element_id"]),
                    centroid=np.array(
                        [
                            float(row["centroid_x_mm"]),
                            float(row["centroid_y_mm"]),
                            float(row["centroid_z_mm"]),
                        ],
                        dtype=float,
                    ),
                    sample_voxels=int(row["sample_voxels"]),
                    bvtv=float(row["bvtv"]),
                    youngs_modulus=float(row["youngs_modulus_mpa"]),
                    threshold=float(row["threshold"]),
                )
            )
    if not rows:
        raise ValueError(f"No rows found in {path}")
    return rows


def normalize_image(data: np.ndarray, vmin: float, vmax: float) -> np.ndarray:
    if vmax <= vmin:
        vmax = vmin + 1.0
    return np.clip((data - vmin) / (vmax - vmin), 0.0, 1.0)


def slice_for_axis(data: np.ndarray, axis: int, index: int) -> np.ndarray:
    if axis == 0:
        return data[index, :, :].T
    if axis == 1:
        return data[:, index, :].T
    if axis == 2:
        return data[:, :, index].T
    raise ValueError(f"Unsupported axis: {axis}")


def projection_for_axis(points_ijk: np.ndarray, axis: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if axis == 0:
        return points_ijk[:, 1], points_ijk[:, 2], points_ijk[:, 0]
    if axis == 1:
        return points_ijk[:, 0], points_ijk[:, 2], points_ijk[:, 1]
    if axis == 2:
        return points_ijk[:, 0], points_ijk[:, 1], points_ijk[:, 2]
    raise ValueError(f"Unsupported axis: {axis}")


def plane_label(axis: int) -> str:
    return {0: "sagittal", 1: "coronal", 2: "axial"}[axis]


def representative_slice(
    points_ijk: np.ndarray,
    axis: int,
    shape: tuple[int, int, int],
    bvtv: np.ndarray | None = None,
) -> int:
    values = points_ijk[:, axis]
    finite = np.isfinite(values)
    if not np.any(finite):
        return shape[axis] // 2

    if bvtv is not None:
        positive = finite & np.isfinite(bvtv) & (bvtv > 0.0)
        if np.any(positive):
            rounded = np.clip(np.round(values[positive]).astype(int), 0, shape[axis] - 1)
            scores = np.bincount(rounded, weights=bvtv[positive], minlength=shape[axis])
            return int(np.argmax(scores))

    values = values[finite]
    median = np.nanmedian(values)
    nearest = values[np.argmin(np.abs(values - median))]
    return int(np.clip(round(float(nearest)), 0, shape[axis] - 1))


def save_slice(
    output_path: Path,
    data: np.ndarray,
    axis: int,
    index: int,
    title: str,
    vmin: float,
    vmax: float,
) -> None:
    image = normalize_image(slice_for_axis(data, axis, index), vmin=vmin, vmax=vmax)
    fig, ax = plt.subplots(figsize=(6, 6), dpi=150)
    ax.imshow(image, cmap="gray", origin="lower")
    ax.set_title(title)
    ax.set_xlabel("voxel index")
    ax.set_ylabel("voxel index")
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def save_mask_slice(output_path: Path, mask: np.ndarray, axis: int, index: int, title: str) -> None:
    fig, ax = plt.subplots(figsize=(6, 6), dpi=150)
    ax.imshow(slice_for_axis(mask, axis, index), cmap="gray", origin="lower", vmin=0, vmax=1)
    ax.set_title(title)
    ax.set_xlabel("voxel index")
    ax.set_ylabel("voxel index")
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def save_overlay(
    output_path: Path,
    data: np.ndarray,
    points_ijk: np.ndarray,
    bvtv: np.ndarray,
    axis: int,
    index: int,
    title: str,
    vmin: float,
    vmax: float,
    slice_tolerance: float,
) -> None:
    image = normalize_image(slice_for_axis(data, axis, index), vmin=vmin, vmax=vmax)
    x, y, z = projection_for_axis(points_ijk, axis)
    keep = np.abs(z - index) <= slice_tolerance
    zero = keep & np.isfinite(bvtv) & (bvtv <= 0.0)
    positive = keep & np.isfinite(bvtv) & (bvtv > 0.0)

    fig, ax = plt.subplots(figsize=(6, 6), dpi=150)
    ax.imshow(image, cmap="gray", origin="lower")
    if np.any(zero):
        ax.scatter(
            x[zero],
            y[zero],
            c="#d0d0d0",
            s=16,
            edgecolors="#4d4d4d",
            linewidths=0.25,
            alpha=0.55,
            label="BV/TV = 0",
        )
    if np.any(positive):
        scatter = ax.scatter(
            x[positive],
            y[positive],
            c=bvtv[positive],
            cmap="viridis",
            vmin=0.0,
            vmax=1.0,
            s=46,
            edgecolors="white",
            linewidths=0.45,
        )
        fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04, label="BV/TV")
    if np.any(zero):
        ax.legend(loc="lower right", fontsize=8, frameon=True)
    ax.set_title(title)
    ax.set_xlabel("voxel index")
    ax.set_ylabel("voxel index")
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def save_histogram(output_path: Path, values: np.ndarray, title: str, xlabel: str) -> None:
    fig, ax = plt.subplots(figsize=(7, 4), dpi=150)
    ax.hist(values[np.isfinite(values)], bins=30, color="#2374ab", edgecolor="white")
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("element count")
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def save_3d_scatter(output_path: Path, points: np.ndarray, bvtv: np.ndarray) -> None:
    fig = plt.figure(figsize=(7, 6), dpi=150)
    ax = fig.add_subplot(111, projection="3d")
    scatter = ax.scatter(
        points[:, 0],
        points[:, 1],
        points[:, 2],
        c=bvtv,
        cmap="viridis",
        vmin=0.0,
        vmax=1.0,
        s=22,
        depthshade=False,
    )
    fig.colorbar(scatter, ax=ax, shrink=0.65, pad=0.08, label="BV/TV")
    ax.set_title("FE element centroids colored by BV/TV")
    ax.set_xlabel("x voxel")
    ax.set_ylabel("y voxel")
    ax.set_zlabel("z voxel")
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def write_html_report(
    path: Path,
    title: str,
    ct_path: Path,
    mesh_path: Path,
    results_path: Path,
    images: dict[str, str],
    rows: list[ResultRow],
    volume_shape: tuple[int, int, int],
    spacing: np.ndarray,
    threshold: float,
) -> None:
    bvtv = np.array([row.bvtv for row in rows], dtype=float)
    modulus = np.array([row.youngs_modulus for row in rows], dtype=float)
    sample_voxels = np.array([row.sample_voxels for row in rows], dtype=int)
    nonzero = int(np.count_nonzero(bvtv > 0))
    full = int(np.count_nonzero(bvtv >= 1.0))
    finite_modulus = modulus[np.isfinite(modulus)]
    max_rows = sorted(rows, key=lambda row: row.bvtv, reverse=True)[:12]

    def img_tag(key: str, caption: str) -> str:
        src = html.escape(images[key])
        return f'<figure><img src="{src}" alt="{html.escape(caption)}"><figcaption>{html.escape(caption)}</figcaption></figure>'

    table_rows = "\n".join(
        "<tr>"
        f"<td>{row.element_id}</td>"
        f"<td>{row.bvtv:.4f}</td>"
        f"<td>{row.youngs_modulus:.2f}</td>"
        f"<td>{row.sample_voxels}</td>"
        "</tr>"
        for row in max_rows
    )

    html_text = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{html.escape(title)}</title>
  <style>
    body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; margin: 0; color: #172026; background: #f7f8fa; }}
    header {{ background: #172026; color: white; padding: 24px 32px; }}
    main {{ padding: 24px 32px 40px; max-width: 1280px; margin: 0 auto; }}
    h1 {{ margin: 0 0 8px; font-size: 28px; }}
    h2 {{ margin: 28px 0 12px; font-size: 20px; }}
    .meta, .stats {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap: 12px; }}
    .box {{ background: white; border: 1px solid #d8dee4; border-radius: 8px; padding: 14px 16px; }}
    .box strong {{ display: block; color: #5b6770; font-size: 12px; text-transform: uppercase; margin-bottom: 6px; }}
    .grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(290px, 1fr)); gap: 16px; }}
    figure {{ background: white; border: 1px solid #d8dee4; border-radius: 8px; margin: 0; padding: 10px; }}
    img {{ width: 100%; height: auto; display: block; }}
    figcaption {{ font-size: 13px; color: #5b6770; margin-top: 8px; }}
    table {{ width: 100%; border-collapse: collapse; background: white; border: 1px solid #d8dee4; border-radius: 8px; overflow: hidden; }}
    th, td {{ padding: 8px 10px; border-bottom: 1px solid #eaeef2; text-align: right; }}
    th:first-child, td:first-child {{ text-align: left; }}
    th {{ background: #eef2f5; color: #172026; }}
    code {{ background: #eef2f5; padding: 2px 5px; border-radius: 5px; }}
  </style>
</head>
<body>
  <header>
    <h1>{html.escape(title)}</h1>
    <div>Input CT, FE mesh, segmentation, and mapped element BV/TV outputs.</div>
  </header>
  <main>
    <h2>Inputs And Outputs</h2>
    <section class="meta">
      <div class="box"><strong>CT volume</strong>{html.escape(str(ct_path))}</div>
      <div class="box"><strong>FE mesh</strong>{html.escape(str(mesh_path))}</div>
      <div class="box"><strong>Mapping CSV</strong>{html.escape(str(results_path))}</div>
      <div class="box"><strong>Volume shape</strong>{volume_shape}</div>
      <div class="box"><strong>Voxel spacing mm</strong>{spacing[0]:.4g}, {spacing[1]:.4g}, {spacing[2]:.4g}</div>
      <div class="box"><strong>Threshold</strong>{threshold:.6g}</div>
    </section>

    <h2>Summary</h2>
    <section class="stats">
      <div class="box"><strong>Elements</strong>{len(rows)}</div>
      <div class="box"><strong>Nonzero BV/TV</strong>{nonzero}</div>
      <div class="box"><strong>BV/TV min mean max</strong>{np.nanmin(bvtv):.4f} / {np.nanmean(bvtv):.4f} / {np.nanmax(bvtv):.4f}</div>
      <div class="box"><strong>Full BV/TV elements</strong>{full}</div>
      <div class="box"><strong>Young modulus min max MPa</strong>{np.nanmin(finite_modulus):.2f} / {np.nanmax(finite_modulus):.2f}</div>
      <div class="box"><strong>Sample voxels min max</strong>{int(sample_voxels.min())} / {int(sample_voxels.max())}</div>
    </section>

    <h2>Input CT Slices</h2>
    <section class="grid">
      {img_tag("ct_axial", "Input CT axial slice")}
      {img_tag("ct_coronal", "Input CT coronal slice")}
      {img_tag("ct_sagittal", "Input CT sagittal slice")}
    </section>

    <h2>Segmentation Used For BV/TV</h2>
    <section class="grid">
      {img_tag("mask_axial", "Segmented axial slice")}
      {img_tag("mask_coronal", "Segmented coronal slice")}
      {img_tag("mask_sagittal", "Segmented sagittal slice")}
    </section>

    <h2>Mapped Output On CT</h2>
    <p>Grey markers indicate nearby sampled FE element centroids with BV/TV = 0. Colored markers indicate positive BV/TV values.</p>
    <section class="grid">
      {img_tag("overlay_axial", "Element centroids near axial slice colored by BV/TV")}
      {img_tag("overlay_coronal", "Element centroids near coronal slice colored by BV/TV")}
      {img_tag("overlay_sagittal", "Element centroids near sagittal slice colored by BV/TV")}
    </section>

    <h2>Output Distributions</h2>
    <section class="grid">
      {img_tag("scatter3d", "3D FE element centroid cloud colored by BV/TV")}
      {img_tag("hist_bvtv", "BV/TV histogram")}
      {img_tag("hist_modulus", "Young modulus histogram")}
    </section>

    <h2>Highest BV/TV Elements</h2>
    <table>
      <thead><tr><th>Element</th><th>BV/TV</th><th>E MPa</th><th>Sample voxels</th></tr></thead>
      <tbody>{table_rows}</tbody>
    </table>
  </main>
</body>
</html>
"""
    path.write_text(html_text, encoding="utf-8")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Create an HTML report for CT-to-FE BV/TV results.")
    parser.add_argument("--ct", required=True, help="Input CT volume, .nii or .nii.gz")
    parser.add_argument("--mesh", required=True, help="Input Abaqus .inp mesh")
    parser.add_argument("--results", required=True, help="CSV output from map_bvtv_to_fe.py")
    parser.add_argument("--out-dir", required=True, help="Output directory for report and PNG assets")
    parser.add_argument("--title", default="CT-to-FE BV/TV Visual Report", help="Report title")
    parser.add_argument("--slice-tolerance", type=float, default=3.0, help="Overlay tolerance in voxel units")
    parser.add_argument("--window-min", type=float, default=None, help="CT display window minimum")
    parser.add_argument("--window-max", type=float, default=None, help="CT display window maximum")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    out_dir = Path(args.out_dir)
    assets_dir = out_dir / "assets"
    assets_dir.mkdir(parents=True, exist_ok=True)

    volume = read_nifti(args.ct)
    parse_abaqus_inp(args.mesh)
    rows = read_results(args.results)
    threshold = rows[0].threshold
    mask, _ = segment_volume(volume.data, threshold=threshold, bone_is_high=True)

    inverse_affine = np.linalg.inv(volume.affine)
    world_centroids = np.vstack([row.centroid for row in rows])
    points_h = np.vstack([world_centroids.T, np.ones(world_centroids.shape[0])])
    points_ijk = (inverse_affine @ points_h)[:3].T
    bvtv = np.array([row.bvtv for row in rows], dtype=float)
    modulus = np.array([row.youngs_modulus for row in rows], dtype=float)

    finite = volume.data[np.isfinite(volume.data)]
    vmin = float(args.window_min) if args.window_min is not None else float(np.percentile(finite, 1.0))
    vmax = float(args.window_max) if args.window_max is not None else float(np.percentile(finite, 99.5))

    slice_indices = {
        0: representative_slice(points_ijk, 0, volume.data.shape, bvtv=bvtv),
        1: representative_slice(points_ijk, 1, volume.data.shape, bvtv=bvtv),
        2: representative_slice(points_ijk, 2, volume.data.shape, bvtv=bvtv),
    }

    images: dict[str, str] = {}
    for axis, name in [(2, "axial"), (1, "coronal"), (0, "sagittal")]:
        index = slice_indices[axis]
        ct_name = f"ct_{name}.png"
        mask_name = f"mask_{name}.png"
        overlay_name = f"overlay_{name}.png"
        save_slice(
            assets_dir / ct_name,
            volume.data,
            axis,
            index,
            f"Input CT {plane_label(axis)} slice {index}",
            vmin,
            vmax,
        )
        save_mask_slice(
            assets_dir / mask_name,
            mask,
            axis,
            index,
            f"Segmented {plane_label(axis)} slice {index}",
        )
        save_overlay(
            assets_dir / overlay_name,
            volume.data,
            points_ijk,
            bvtv,
            axis,
            index,
            f"BV/TV overlay on {plane_label(axis)} slice {index}",
            vmin,
            vmax,
            slice_tolerance=args.slice_tolerance,
        )
        images[f"ct_{name}"] = f"assets/{ct_name}"
        images[f"mask_{name}"] = f"assets/{mask_name}"
        images[f"overlay_{name}"] = f"assets/{overlay_name}"

    save_3d_scatter(assets_dir / "bvtv_3d_scatter.png", points_ijk, bvtv)
    save_histogram(assets_dir / "bvtv_histogram.png", bvtv, "BV/TV distribution", "BV/TV")
    save_histogram(
        assets_dir / "youngs_modulus_histogram.png",
        modulus[np.isfinite(modulus)],
        "Young modulus distribution",
        "E (MPa)",
    )
    images["scatter3d"] = "assets/bvtv_3d_scatter.png"
    images["hist_bvtv"] = "assets/bvtv_histogram.png"
    images["hist_modulus"] = "assets/youngs_modulus_histogram.png"

    report_path = out_dir / "report.html"
    write_html_report(
        report_path,
        args.title,
        Path(args.ct),
        Path(args.mesh),
        Path(args.results),
        images,
        rows,
        volume.data.shape,
        volume.spacing,
        threshold,
    )

    print(f"Wrote report: {report_path}")
    print(f"Wrote assets: {assets_dir}")
    print(f"Slice indices: sagittal={slice_indices[0]}, coronal={slice_indices[1]}, axial={slice_indices[2]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
