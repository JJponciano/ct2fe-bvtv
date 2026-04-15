# CT-to-FE BV/TV Mapping

**Implementation author:** Dr. Jean-Jacques Ponciano

This repository provides a non-official Python implementation of a computed-tomography-to-finite-element preprocessing workflow for assigning local bone volume fraction (`BV/TV`) values to finite element (FE) meshes.

The implementation is inspired by the material mapping procedure described in the referenced paper by Vautrin et al. The objective is to provide a transparent, reproducible, and inspectable research codebase for:

- segmenting volumetric CT data;
- registering an FE mesh to a CT coordinate system;
- computing element-wise local `BV/TV` from segmented CT voxels;
- converting `BV/TV` values into elastic material properties;
- exporting tabular and Abaqus-compatible material assignment files;
- visualizing input CT data, segmentation, FE element locations, and mapped output values;
- computing an article-style peri-implant `BV/TV` value in a hollow cylindrical region.

## Disclaimer

This repository is a **non-official implementation**. It is not authored, reviewed, endorsed, or validated by the authors of the original paper. The code is intended for research, reproducibility, and methodological exploration. It does not reproduce the complete nonlinear finite element model of the paper and should not be used as clinical or regulatory software.

## Reference Paper

The implementation is based on the method described in:

Vautrin, A., Thierrin, R., Wili, P., Voumard, B., Klingler, S., Chappuis, V., Varga, P., & Zysset, P. (2024). *Homogenized finite element simulations can predict the primary stability of dental implants in human jawbone*. Journal of the Mechanical Behavior of Biomedical Materials, 158, 106688. <https://doi.org/10.1016/j.jmbbm.2024.106688>

When using this repository, cite both this software repository and the original article.

## Scope Of The Implementation

The repository implements the preprocessing part of the workflow:

```text
CT volume -> segmentation -> CT/FE registration -> element-wise BV/TV -> FE material export
```

The default material mapping follows the relation reported in the article:

```text
E = 8534.64 * (BV/TV)^1.63 MPa
nu = 0.246
```

The default local mapping sphere diameter is `1.25 mm`, matching the article. The demonstration uses a larger sphere because the public demonstration CT has much coarser voxel spacing than the micro-CT data used in the article.

## Repository Structure

```text
.
├── README.md
├── CITATION.cff
├── LICENSE
├── requirements.txt
├── environment.yml
├── pyproject.toml
├── data/
│   ├── README.md
│   ├── CT_Abdo.nii.gz
│   ├── demo_abdomen_mesh.inp
│   └── demo_registration_landmarks.json
├── docs/
│   ├── DATA_REQUIREMENTS.md
│   ├── IMPLEMENTATION_NOTES.md
│   └── REPRODUCIBILITY.md
├── examples/
│   ├── run_demo.sh
│   ├── landmarks_template.json
│   └── transform_template.json
├── outputs/
│   ├── README.md
│   ├── demo_abdomen_bvtv.csv
│   ├── demo_abdomen_bvtv.materials.inp
│   ├── demo_abdomen_bvtv_otsu.csv
│   ├── demo_abdomen_bvtv_otsu.materials.inp
│   ├── demo_abdomen_bvtv_registered.csv
│   ├── demo_abdomen_bvtv_registered.materials.inp
│   ├── demo_peri_implant_bvtv.csv
│   ├── verification_report.md
│   ├── visualization_300/
│   └── visualization_otsu/
├── scripts/
│   ├── download_demo_data.py
│   ├── make_demo_inp.py
│   ├── map_bvtv_to_fe.py
│   ├── verify_pipeline.py
│   └── visualize_bvtv.py
└── tests/
    └── test_mapper_synthetic.py
```

## Component Roles

- `scripts/map_bvtv_to_fe.py`: main pipeline. Reads CT data and an Abaqus mesh, computes element-wise `BV/TV`, writes CSV and Abaqus material include files, and optionally computes peri-implant `BV/TV`.
- `scripts/make_demo_inp.py`: creates a regular Abaqus C3D8 demonstration mesh in the CT coordinate system.
- `scripts/visualize_bvtv.py`: creates HTML/PNG reports showing CT slices, segmentation, mapped FE element centroids, histograms, and 3D centroid plots.
- `scripts/verify_pipeline.py`: performs deterministic verification checks on synthetic cases and demonstration outputs.
- `scripts/download_demo_data.py`: downloads the public demonstration CT volume.
- `tests/test_mapper_synthetic.py`: minimal regression test for expected full-bone and empty-bone `BV/TV` values.
- `examples/run_demo.sh`: runs the complete demonstration and verification workflow.

## Installation

### Hardware

The demonstration data can be processed on a standard laptop or workstation. No GPU is required. For typical demonstration use, 2-4 GB of available RAM is sufficient. Larger micro-CT volumes may require substantially more memory depending on voxel count and mesh size.

### Python Environment

Python `3.10` or newer is recommended. The repository was tested with Python `3.11`.

Using `venv`:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Using Conda:

```bash
conda env create -f environment.yml
conda activate ct2fe-bvtv
```

The required Python packages are:

- `numpy`
- `scikit-image`
- `matplotlib`

No stochastic algorithm is used in the current implementation. No random seed is required to reproduce the demonstration outputs.

## Quick Reproduction

Run the full demonstration:

```bash
bash examples/run_demo.sh
```

This command downloads the demonstration CT if needed, generates the demonstration mesh, runs fixed-threshold and Otsu-threshold mappings, computes a peri-implant demonstration `BV/TV`, creates visual reports, and runs verification checks.

Expected successful terminal output includes:

```text
Synthetic BV/TV mapping test passed
PASS: Synthetic full/empty bone test passed
PASS: Anisotropic affine sphere test passed
PASS: Rigid landmark registration test passed
PASS: Peri-implant hollow-cylinder BV/TV test passed
PASS: threshold-300: BV/TV range, sampled voxel counts, nu, E formula, and Abaqus modulus all pass
PASS: otsu: BV/TV range, sampled voxel counts, nu, E formula, and Abaqus modulus all pass
```

The verification report is written to:

```text
outputs/verification_report.md
```

## Step-By-Step Usage

### 1. Download The Demonstration CT

```bash
python3 scripts/download_demo_data.py
```

The demonstration CT is a public low-resolution abdominal CT image from the NiiVue sample image collection. It is used only to exercise the software. It is not dental micro-CT data and is not expected to reproduce the numerical values of the referenced paper.

### 2. Generate A Demonstration FE Mesh

```bash
python3 scripts/make_demo_inp.py \
  --ct data/CT_Abdo.nii.gz \
  --out data/demo_abdomen_mesh.inp \
  --elements 8,8,6 \
  --anchor-threshold 300 \
  --padding-voxels 4
```

Expected output:

```text
data/demo_abdomen_mesh.inp
```

This file is a regular Abaqus C3D8 mesh generated in the CT coordinate system.

### 3. Map BV/TV With A Fixed Threshold

```bash
python3 scripts/map_bvtv_to_fe.py \
  --ct data/CT_Abdo.nii.gz \
  --mesh data/demo_abdomen_mesh.inp \
  --out-prefix outputs/demo_abdomen_bvtv \
  --threshold 300 \
  --sphere-diameter-mm 8
```

Expected outputs:

```text
outputs/demo_abdomen_bvtv.csv
outputs/demo_abdomen_bvtv.materials.inp
```

### 4. Map BV/TV With Otsu Thresholding

```bash
python3 scripts/map_bvtv_to_fe.py \
  --ct data/CT_Abdo.nii.gz \
  --mesh data/demo_abdomen_mesh.inp \
  --out-prefix outputs/demo_abdomen_bvtv_otsu \
  --threshold otsu \
  --sphere-diameter-mm 8
```

Expected outputs:

```text
outputs/demo_abdomen_bvtv_otsu.csv
outputs/demo_abdomen_bvtv_otsu.materials.inp
```

### 5. Use Registration Landmarks

If the mesh is not already in the CT/world coordinate system, provide corresponding landmarks:

```bash
python3 scripts/map_bvtv_to_fe.py \
  --ct data/CT_Abdo.nii.gz \
  --mesh data/demo_abdomen_mesh.inp \
  --out-prefix outputs/demo_abdomen_bvtv_registered \
  --threshold 300 \
  --sphere-diameter-mm 8 \
  --landmarks-json data/demo_registration_landmarks.json
```

The landmark JSON format is:

```json
{
  "mesh": [[x1, y1, z1], [x2, y2, z2], [x3, y3, z3]],
  "ct": [[x1, y1, z1], [x2, y2, z2], [x3, y3, z3]]
}
```

The estimated rigid transform maps FE mesh coordinates into CT/world coordinates.

Alternatively, provide a 4x4 transform:

```bash
--transform-json examples/transform_template.json
```

### 6. Compute Peri-Implant BV/TV

The article also reports `BV/TV` in a hollow cylindrical peri-implant region. This can be computed by providing the implant or pilot-hole centerline:

```bash
python3 scripts/map_bvtv_to_fe.py \
  --ct data/CT_Abdo.nii.gz \
  --mesh data/demo_abdomen_mesh.inp \
  --out-prefix outputs/demo_abdomen_bvtv_registered \
  --threshold 300 \
  --sphere-diameter-mm 8 \
  --peri-implant-out outputs/demo_peri_implant_bvtv.csv \
  --implant-centerline-start=-12.7902,-71.5878,-31.8359 \
  --implant-centerline-end=-12.7902,-71.5878,-13.9062 \
  --implant-diameter-mm 4.0 \
  --peri-implant-shell-thickness-mm 1.0
```

Expected output:

```text
outputs/demo_peri_implant_bvtv.csv
```

### 7. Generate Visual Reports

```bash
python3 scripts/visualize_bvtv.py \
  --ct data/CT_Abdo.nii.gz \
  --mesh data/demo_abdomen_mesh.inp \
  --results outputs/demo_abdomen_bvtv.csv \
  --out-dir outputs/visualization_300 \
  --title "CT-to-FE BV/TV Report - threshold 300 HU"
```

Open:

```text
outputs/visualization_300/report.html
```

The report contains CT slices, segmentation slices, `BV/TV` overlays, histograms, and a 3D centroid plot.

### 8. Verify The Pipeline

```bash
python3 tests/test_mapper_synthetic.py
python3 scripts/verify_pipeline.py
```

The verification script checks:

- expected synthetic full-bone and empty-bone cases;
- anisotropic affine voxel handling;
- rigid landmark registration;
- hollow-cylinder peri-implant `BV/TV`;
- CSV row count versus FE element count;
- `BV/TV` range `[0, 1]`;
- material relation `E = 8534.64 * BV/TV^1.63`;
- Abaqus material export consistency;
- selected element recomputation from CT and mesh.

## Expected Data Flow

```text
Input CT volume
    |
    v
CT segmentation
    |
    v
Input FE mesh ---- optional landmarks/transform ----> mesh in CT/world coordinates
    |
    v
Element centroids
    |
    v
Spherical CT sampling around each centroid
    |
    v
Element-wise BV/TV
    |
    v
Material relation
    |
    v
CSV output + Abaqus material include + visualization report
```

## Expected Input Formats

### CT Volume

Currently supported:

- NIfTI-1 `.nii`
- NIfTI-1 `.nii.gz`

The NIfTI affine is used to map voxel indices to physical CT/world coordinates. If your data are DICOM, TIFF stack, MHD/RAW, or Scanco export, convert them to NIfTI before using this repository or extend the reader.

### FE Mesh

Currently supported:

- Abaqus `.inp` files containing `*Node` and `*Element` blocks.

The parser expects node coordinates in physical units, preferably millimetres. The mesh must represent the bone domain to which `BV/TV`-dependent material properties should be assigned.

### Registration

At least one of the following must be true:

- CT and FE mesh are already in the same coordinate system;
- a 4x4 mesh-to-CT transform is available;
- corresponding landmarks are available in mesh and CT/world coordinates.

## Adapting The Pipeline To Custom Data

For a new dataset, modify:

- `--ct`: path to the CT volume;
- `--mesh`: path to the Abaqus mesh;
- `--out-prefix`: output file prefix;
- `--threshold`: numeric threshold or `otsu`;
- `--sphere-diameter-mm`: local mapping sphere diameter;
- `--landmarks-json` or `--transform-json`: registration information when needed;
- `--implant-centerline-start`, `--implant-centerline-end`, `--implant-diameter-mm`: peri-implant `BV/TV` parameters when needed.

For high-resolution micro-CT data comparable to the paper, `--sphere-diameter-mm 1.25` is appropriate. For coarse clinical CT data, a larger sphere may be needed to avoid undersampling.

## Reproducibility Requirements

To reproduce the demonstration exactly:

1. Use Python `3.10` or newer.
2. Install dependencies from `requirements.txt` or `environment.yml`.
3. Use the public `CT_Abdo.nii.gz` demonstration volume downloaded by `scripts/download_demo_data.py`.
4. Run `bash examples/run_demo.sh`.
5. Confirm that `scripts/verify_pipeline.py` reports all checks as `PASS`.

No random seed is required because the pipeline is deterministic.

## Limitations And Assumptions

### Technical Limitations

- The CT reader currently supports NIfTI-1 files only.
- The mesh reader currently supports basic Abaqus `.inp` `*Node` and `*Element` blocks.
- The Abaqus export creates per-element elastic material definitions. It does not create a complete simulation model.
- Large micro-CT volumes may require memory optimization beyond the current implementation.
- Registration is rigid by default. Optional uniform scaling can be enabled, but deformable registration is not implemented.

### Methodological Limitations

- The repository does not include the nonlinear Abaqus UMAT used in the paper.
- Plasticity, damage evolution, contact, loading, and boundary conditions from the full paper are not reproduced.
- The demonstration data are clinical abdominal CT data, not jawbone micro-CT data.
- The demonstration mesh is a generated regular grid, not a dental implant/jawbone model.
- Quantitative output from the demonstration should not be compared to the paper's reported experimental or FE results.

### Differences From The Original Paper

- The original work used 24.6 micrometre isotropic micro-CT data; the demonstration uses approximately 1.49 mm clinical CT voxels.
- The original work used a sample-specific dental implant/jawbone FE setup; the demonstration uses a generated grid.
- The original work used a nonlinear constitutive model with plasticity and damage; this repository maps elastic properties only.
- The original work estimated sample-specific registration from anatomical/experimental landmarks; this repository provides a general landmark-based registration interface.

### Potential Sources Of Error

- Incorrect CT/FE registration.
- Wrong voxel spacing or affine metadata.
- Inappropriate segmentation threshold.
- Mesh units inconsistent with CT physical units.
- Sampling sphere too small relative to image resolution.
- Inclusion of non-bone elements in the mapping.
- Partial volume effects in coarse CT data.

## Citation

If using this repository, cite:

```text
Dr. Jean-Jacques Ponciano. CT-to-FE BV/TV Mapping: A Non-Official Implementation. 2026.
```

Also cite the original paper:

```text
Vautrin A, Thierrin R, Wili P, Voumard B, Klingler S, Chappuis V, Varga P, Zysset P.
Homogenized finite element simulations can predict the primary stability of dental implants in human jawbone.
Journal of the Mechanical Behavior of Biomedical Materials. 2024;158:106688.
doi:10.1016/j.jmbbm.2024.106688
```

## License

This repository is distributed under the MIT License. See `LICENSE`.
