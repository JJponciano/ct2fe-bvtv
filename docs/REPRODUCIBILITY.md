# Reproducibility

**Implementation author:** Dr. Jean-Jacques Ponciano

This document describes how to reproduce the demonstration outputs.

The workflow is a non-official implementation of the preprocessing concepts described in the referenced paper.

## Environment

Create a Python environment using either:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

or:

```bash
conda env create -f environment.yml
conda activate ct2fe-bvtv
```

## Full Demonstration

Run:

```bash
bash examples/run_demo.sh
```

The workflow will:

1. download the demonstration CT if missing;
2. create the demonstration FE mesh;
3. run fixed-threshold mapping;
4. run Otsu-threshold mapping;
5. run mapping through the registration interface;
6. compute a peri-implant demonstration `BV/TV`;
7. create HTML reports;
8. run synthetic and output consistency checks.

## High-Quality Offline Demonstration

For a longer, higher-resolution demonstration, run:

```bash
bash examples/run_high_quality_offline_demo.sh
```

This offline script expects `data/CT_Abdo.nii.gz` to already exist locally. It uses a finer default mesh (`32,24,32`, or `24,576` C3D8 elements), high-DPI visualization assets, and full-element final FE mesh rendering. Outputs are written under `outputs/high_quality/`.

## Determinism

The current implementation is deterministic. No random sampling or stochastic optimization is used. Fixed seeds are therefore not required.

## Expected Files

After a successful run, the following files should exist:

```text
outputs/demo_abdomen_bvtv.csv
outputs/demo_abdomen_bvtv.materials.inp
outputs/demo_abdomen_bvtv.mapped.inp
outputs/demo_abdomen_bvtv_otsu.csv
outputs/demo_abdomen_bvtv_otsu.materials.inp
outputs/demo_abdomen_bvtv_otsu.mapped.inp
outputs/demo_abdomen_bvtv_registered.csv
outputs/demo_abdomen_bvtv_registered.materials.inp
outputs/demo_abdomen_bvtv_registered.mapped.inp
outputs/demo_peri_implant_bvtv.csv
outputs/visualization_300/report.html
outputs/visualization_300/assets/mapped_fe_bvtv_3d.png
outputs/visualization_300/assets/mapped_fe_positive_bvtv_3d.png
outputs/visualization_300/assets/mapped_fe_modulus_3d.png
outputs/visualization_otsu/report.html
outputs/visualization_otsu/assets/mapped_fe_bvtv_3d.png
outputs/visualization_otsu/assets/mapped_fe_positive_bvtv_3d.png
outputs/visualization_otsu/assets/mapped_fe_modulus_3d.png
outputs/verification_report.md
```

## Verification

Run:

```bash
python3 scripts/verify_pipeline.py
```

Expected output includes:

```text
PASS: Synthetic full/empty bone test passed
PASS: Anisotropic affine sphere test passed
PASS: Rigid landmark registration test passed
PASS: Peri-implant hollow-cylinder BV/TV test passed
```

The verification report is written to:

```text
outputs/verification_report.md
```
