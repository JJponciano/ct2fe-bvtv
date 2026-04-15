# Verification Report

## Automated Checks

- PASS: Synthetic full/empty bone test passed
- PASS: Anisotropic affine sphere test passed
- PASS: Rigid landmark registration test passed
- PASS: Peri-implant hollow-cylinder BV/TV test passed
- PASS: threshold-300: row count matches mesh elements (3072)
- PASS: threshold-300: Abaqus material count matches CSV rows (3072)
- PASS: threshold-300: self-contained mapped Abaqus input contains the source FE mesh (3072 elements)
- PASS: threshold-300: self-contained mapped Abaqus input material sections match the material include (3072 sections)
- PASS: threshold-300: threshold is consistent (300)
- PASS: threshold-300: BV/TV range, sampled voxel counts, nu, E formula, and Abaqus modulus all pass
- PASS: threshold-300: selected elements recompute exactly from CT and mesh
- PASS: otsu: row count matches mesh elements (3072)
- PASS: otsu: Abaqus material count matches CSV rows (3072)
- PASS: otsu: self-contained mapped Abaqus input contains the source FE mesh (3072 elements)
- PASS: otsu: self-contained mapped Abaqus input material sections match the material include (3072 sections)
- PASS: otsu: threshold is consistent (-441.004)
- PASS: otsu: BV/TV range, sampled voxel counts, nu, E formula, and Abaqus modulus all pass
- PASS: otsu: selected elements recompute exactly from CT and mesh
- PASS: registered-threshold-300: row count matches mesh elements (3072)
- PASS: registered-threshold-300: Abaqus material count matches CSV rows (3072)
- PASS: registered-threshold-300: self-contained mapped Abaqus input contains the source FE mesh (3072 elements)
- PASS: registered-threshold-300: self-contained mapped Abaqus input material sections match the material include (3072 sections)
- PASS: registered-threshold-300: threshold is consistent (300)
- PASS: registered-threshold-300: BV/TV range, sampled voxel counts, nu, E formula, and Abaqus modulus all pass
- PASS: registered-threshold-300: selected elements recompute exactly from CT and mesh
- PASS: peri-implant output is valid (BV/TV=0.541667, sample_voxels=96)
- PASS: visualization-threshold-300: visual report contains the final mapped FE mesh section
- PASS: visualization-threshold-300: final FE mesh BV/TV and modulus visualization assets exist
- PASS: visualization-otsu: visual report contains the final mapped FE mesh section
- PASS: visualization-otsu: final FE mesh BV/TV and modulus visualization assets exist

## Agreement With The Article

- PASS: The core mapping operation follows the article: segment CT, compute element centroid, sample a spherical region, average binary voxels to get local BV/TV.
- PASS: Otsu thresholding is implemented, matching the article's global Otsu segmentation step.
- PASS: The default sampling sphere diameter is 1.25 mm, matching the article.
- PASS: The elastic relation uses the article constants E0 = 8534.64 MPa, k = 1.63, and nu = 0.246.
- PASS: Optional rigid landmark registration is implemented for mapping mesh coordinates into CT/world coordinates.
- PASS: Self-contained mapped Abaqus input files are generated and verified to contain the source FE mesh plus element-wise material sections.
- PASS: Article-style peri-implant BV/TV in a 1 mm hollow cylindrical shell can be computed from an implant/pilot-hole centerline.
- LIMITATION: The demo data are public clinical abdomen CT data, not the article's 24.6 micrometer jawbone micro-CT scans.
- LIMITATION: The demo mesh is a generated regular Abaqus grid, not the article's dental implant/jawbone geometry.
- LIMITATION: The demo mapping does not use real experimental landmarks; when no landmarks or transform are provided, the pipeline assumes CT and FE mesh coordinates are already aligned.
- LIMITATION: The nonlinear Abaqus UMAT plasticity and damage model from the article is not reproduced; the export is elastic material assignment by element.
- LIMITATION: The demo uses an 8 mm sphere because the public CT voxel spacing is about 1.49 mm. The article's 1.25 mm sphere is appropriate for their high-resolution micro-CT, not for this low-resolution demo CT.

Conclusion: the implementation is algorithmically consistent with the article's BV/TV-to-element material mapping step, but the demo results are not quantitatively comparable to the paper's published jawbone/implant results.
