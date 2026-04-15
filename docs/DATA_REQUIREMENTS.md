# Data Requirements

**Implementation author:** Dr. Jean-Jacques Ponciano

This document describes the data required to apply the CT-to-FE `BV/TV` mapping workflow to a custom study.

The repository is a non-official implementation of the referenced paper and should be treated as independent research software.

## Required Inputs

### 1. CT Volume

The current implementation accepts NIfTI-1 files:

- `.nii`
- `.nii.gz`

The volume must contain valid voxel spacing and affine metadata. The affine is used to convert voxel indices to physical CT/world coordinates.

For other source formats, convert to NIfTI before running the pipeline:

- DICOM stack
- TIFF stack
- MHD/RAW
- Scanco export

The following metadata should be known:

- voxel size;
- physical units;
- intensity scale;
- whether intensity values are calibrated;
- segmentation threshold, if predetermined.

## 2. FE Mesh

The current implementation accepts Abaqus `.inp` files containing:

- `*Node` block;
- `*Element` block;
- node coordinates in physical units, preferably millimetres.

The mesh should represent the bone domain to which local `BV/TV` values will be assigned. If the input model contains implant, embedding, or loading components, those regions should be separated by element sets or provided as a bone-only mesh for mapping.

## 3. Registration Information

The CT volume and FE mesh must be expressed in the same physical coordinate frame before local `BV/TV` values can be assigned.

The workflow supports three cases:

1. CT and FE mesh are already aligned.
2. A 4x4 mesh-to-CT transform is provided.
3. Corresponding landmarks are provided in mesh and CT/world coordinates.

Landmarks should be non-collinear and should capture the relevant anatomical and experimental orientation. Examples include:

- pilot-hole or implant centerline;
- bone level;
- oral marker;
- mesial marker;
- loading direction marker.

## 4. Implant Geometry For Peri-Implant BV/TV

To compute peri-implant `BV/TV` in a hollow cylindrical region, provide:

- implant or pilot-hole centerline start point;
- implant or pilot-hole centerline end point;
- implant diameter;
- shell thickness around the implant diameter.

The referenced paper used a 1 mm thick hollow cylindrical shell.

## Output Data

The mapper produces:

- element-wise CSV with centroid coordinates, sampled voxel count, `BV/TV`, Young's modulus, Poisson's ratio, and threshold;
- Abaqus material include file;
- optional peri-implant `BV/TV` CSV;
- optional HTML/PNG visualization report.

## Data Quality Considerations

The most important sources of error are:

- CT/FE misregistration;
- incorrect voxel spacing;
- inconsistent mesh and CT units;
- inappropriate segmentation threshold;
- undersampling due to a sphere diameter smaller than the voxel spacing;
- mapping non-bone elements;
- partial-volume effects in low-resolution images.
