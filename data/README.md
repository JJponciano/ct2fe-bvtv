# Data Directory

**Implementation author:** Dr. Jean-Jacques Ponciano

This directory contains input data used by the demonstration workflow.

## Demonstration CT

`CT_Abdo.nii.gz` is a public low-resolution CT demonstration volume downloaded from the NiiVue sample image collection:

<https://neurolabusc.github.io/niivue-images/>

The file is used only to exercise the software pipeline. It is not dental micro-CT data and is not intended to reproduce the numerical results of the referenced dental implant paper.

If the file is not present, run:

```bash
python3 scripts/download_demo_data.py
```

## Generated Demonstration Mesh

`demo_abdomen_mesh.inp` is generated from the CT bounds by:

```bash
python3 scripts/make_demo_inp.py \
  --ct data/CT_Abdo.nii.gz \
  --out data/demo_abdomen_mesh.inp
```

It is a regular Abaqus C3D8 grid used for software validation and visualization.

## Registration Example

`demo_registration_landmarks.json` is an identity registration example because the demonstration mesh is generated directly in the CT/world coordinate frame.
