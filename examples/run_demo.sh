#!/usr/bin/env bash
set -euo pipefail

python3 scripts/download_demo_data.py

python3 scripts/make_demo_inp.py \
  --ct data/CT_Abdo.nii.gz \
  --out data/demo_abdomen_mesh.inp \
  --elements 16,12,16 \
  --anchor-threshold 300 \
  --padding-voxels 4

python3 scripts/map_bvtv_to_fe.py \
  --ct data/CT_Abdo.nii.gz \
  --mesh data/demo_abdomen_mesh.inp \
  --out-prefix outputs/demo_abdomen_bvtv \
  --threshold 300 \
  --sphere-diameter-mm 8

python3 scripts/map_bvtv_to_fe.py \
  --ct data/CT_Abdo.nii.gz \
  --mesh data/demo_abdomen_mesh.inp \
  --out-prefix outputs/demo_abdomen_bvtv_otsu \
  --threshold otsu \
  --sphere-diameter-mm 8

python3 scripts/map_bvtv_to_fe.py \
  --ct data/CT_Abdo.nii.gz \
  --mesh data/demo_abdomen_mesh.inp \
  --out-prefix outputs/demo_abdomen_bvtv_registered \
  --threshold 300 \
  --sphere-diameter-mm 8 \
  --landmarks-json data/demo_registration_landmarks.json \
  --peri-implant-out outputs/demo_peri_implant_bvtv.csv \
  --implant-centerline-start=-12.7902,-71.5878,-31.8359 \
  --implant-centerline-end=-12.7902,-71.5878,-13.9062 \
  --implant-diameter-mm 4.0 \
  --peri-implant-shell-thickness-mm 1.0

python3 scripts/visualize_bvtv.py \
  --ct data/CT_Abdo.nii.gz \
  --mesh data/demo_abdomen_mesh.inp \
  --mapped-inp outputs/demo_abdomen_bvtv.mapped.inp \
  --results outputs/demo_abdomen_bvtv.csv \
  --out-dir outputs/visualization_300 \
  --title "CT-to-FE BV/TV Report - threshold 300 HU"

python3 scripts/visualize_bvtv.py \
  --ct data/CT_Abdo.nii.gz \
  --mesh data/demo_abdomen_mesh.inp \
  --mapped-inp outputs/demo_abdomen_bvtv_otsu.mapped.inp \
  --results outputs/demo_abdomen_bvtv_otsu.csv \
  --out-dir outputs/visualization_otsu \
  --title "CT-to-FE BV/TV Report - Otsu threshold"

python3 tests/test_mapper_synthetic.py
python3 scripts/verify_pipeline.py
