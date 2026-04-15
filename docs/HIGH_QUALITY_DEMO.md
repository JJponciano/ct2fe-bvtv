# High-Quality Offline Demonstration

**Implementation author:** Dr. Jean-Jacques Ponciano

This document describes the offline high-quality demonstration workflow. It is intended for visual inspection of the final mapped finite element input with a substantially finer mesh than the quick demonstration.

The workflow is offline: it does not download data. It requires the demonstration CT volume to already exist at `data/CT_Abdo.nii.gz`.

## Run

```bash
bash examples/run_high_quality_offline_demo.sh
```

The script writes all high-quality outputs under:

```text
outputs/high_quality/
```

This directory is ignored by Git because the generated Abaqus inputs, CSV files, and high-resolution PNG assets can be large.

## Default Quality Settings

The default high-quality settings are:

```text
ELEMENTS=32,24,32
SPHERE_DIAMETER_MM=8
IMAGE_DPI=240
MAX_ELEMENTS_RENDERED=999999
```

This generates `24,576` C3D8 finite elements and renders all elements in the final mesh visualization. The quick demonstration uses `16,12,16`, or `3,072` elements.

## Heavier Settings

For a longer run with smaller elements, use:

```bash
ELEMENTS=48,36,48 IMAGE_DPI=300 bash examples/run_high_quality_offline_demo.sh
```

This generates `82,944` C3D8 elements. It can take substantially longer and produces larger output files. The visual and numerical quality is still limited by the resolution and anatomy of the public demonstration CT; this dataset is not the dental micro-CT dataset used in the referenced article.

## Main Outputs

After a successful run, inspect:

```text
outputs/high_quality/visualization_threshold_300/report.html
outputs/high_quality/visualization_otsu/report.html
outputs/high_quality/demo_abdomen_bvtv_hq.mapped.inp
outputs/high_quality/demo_abdomen_bvtv_hq_otsu.mapped.inp
outputs/high_quality/demo_abdomen_bvtv_hq_registered.mapped.inp
outputs/high_quality/verification_report_hq.md
```

The HTML reports include:

- CT slices;
- segmentation slices;
- mapped centroids on CT slices;
- element-wise `BV/TV` and Young's modulus histograms;
- the final mapped FE mesh colored by `BV/TV`;
- the final mapped FE mesh showing only elements with `BV/TV > 0`;
- the final mapped FE mesh colored by Young's modulus.

## Custom Data

To use another local CT volume and output directory:

```bash
CT_PATH=/path/to/ct.nii.gz \
OUT_ROOT=outputs/my_high_quality_case \
ELEMENTS=48,36,48 \
bash examples/run_high_quality_offline_demo.sh
```

For a real study, replace the demonstration mesh generation step with a study-specific Abaqus mesh and provide study-specific registration landmarks or a 4x4 transform. The included high-quality script is a complete demonstration pipeline, not a substitute for a validated anatomical FE model.
