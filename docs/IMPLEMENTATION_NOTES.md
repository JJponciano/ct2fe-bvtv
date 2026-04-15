# Implementation Notes

**Implementation author:** Dr. Jean-Jacques Ponciano

This document summarizes the implemented method and its relationship to the referenced paper.

The repository is a non-official implementation and does not constitute an endorsed reproduction by the authors of the referenced paper.

## Implemented Workflow

The central operation is local assignment of CT-derived `BV/TV` values to FE elements:

1. Read the CT volume.
2. Segment the CT image by a numeric threshold or Otsu thresholding.
3. Read the FE mesh.
4. Optionally register mesh coordinates into CT/world coordinates.
5. Compute the centroid of each FE element.
6. Sample segmented voxels within a sphere centred at the element centroid.
7. Compute local `BV/TV` as the mean binary voxel value in that sphere.
8. Convert `BV/TV` to elastic modulus:

```text
E = 8534.64 * (BV/TV)^1.63 MPa
nu = 0.246
```

9. Export CSV and Abaqus material definitions.

## Registration

Rigid registration is estimated from corresponding landmarks using the Kabsch algorithm. The estimated transform maps FE mesh coordinates into CT/world coordinates.

Optional uniform scaling can be enabled, but it is disabled by default because FE and CT data should usually share physical units.

## Peri-Implant BV/TV

The implementation can compute `BV/TV` in a hollow cylindrical shell around an implant or pilot-hole centerline. This follows the article's description of a peri-implant hollow cylindrical region. The default shell thickness used in the examples is 1 mm.

## Abaqus Export

The generated `.materials.inp` file contains one element set and one elastic material per element. This is simple and explicit, but it is not the most compact representation for very large meshes.

For large production models, users may prefer to bin `BV/TV` values into material classes or use field variables in Abaqus.

## Differences From The Original Study

This repository implements preprocessing and material mapping, not the complete study. In particular, it does not include:

- the original jawbone micro-CT data;
- the original implant/jawbone FE meshes;
- the Abaqus UMAT;
- nonlinear plasticity and damage simulation;
- insertion damage modelling;
- experimental loading protocol;
- experimental validation against ultimate load.

The demonstration data are intentionally used only to test the computational pipeline.
