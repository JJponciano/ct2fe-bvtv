#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

# Offline high-quality demonstration.
# It does not download data. It expects the CT volume to already exist locally.
#
# Typical use:
#   bash examples/run_high_quality_offline_demo.sh
#
# Optional overrides:
#   ELEMENTS=48,36,48 IMAGE_DPI=300 bash examples/run_high_quality_offline_demo.sh

CT_PATH="${CT_PATH:-data/CT_Abdo.nii.gz}"
LANDMARKS_PATH="${LANDMARKS_PATH:-data/demo_registration_landmarks.json}"
OUT_ROOT="${OUT_ROOT:-outputs/high_quality}"
MESH_PATH="${MESH_PATH:-${OUT_ROOT}/demo_abdomen_mesh_hq.inp}"

ELEMENTS="${ELEMENTS:-32,24,32}"
FIXED_THRESHOLD="${FIXED_THRESHOLD:-300}"
SPHERE_DIAMETER_MM="${SPHERE_DIAMETER_MM:-8}"
ANCHOR_THRESHOLD="${ANCHOR_THRESHOLD:-300}"
PADDING_VOXELS="${PADDING_VOXELS:-4}"

IMAGE_DPI="${IMAGE_DPI:-240}"
MAX_ELEMENTS_RENDERED="${MAX_ELEMENTS_RENDERED:-999999}"

FIXED_PREFIX="${OUT_ROOT}/demo_abdomen_bvtv_hq"
OTSU_PREFIX="${OUT_ROOT}/demo_abdomen_bvtv_hq_otsu"
REGISTERED_PREFIX="${OUT_ROOT}/demo_abdomen_bvtv_hq_registered"
PERI_IMPLANT_OUT="${OUT_ROOT}/demo_peri_implant_bvtv_hq.csv"
VERIFY_REPORT="${OUT_ROOT}/verification_report_hq.md"

if [[ ! -f "$CT_PATH" ]]; then
  echo "Missing CT volume: ${CT_PATH}" >&2
  echo "This script is offline and will not download data." >&2
  echo "Run 'python3 scripts/download_demo_data.py' once, then rerun this script." >&2
  exit 1
fi

if [[ ! -f "$LANDMARKS_PATH" ]]; then
  echo "Missing landmarks file: ${LANDMARKS_PATH}" >&2
  exit 1
fi

mkdir -p "$OUT_ROOT"

run_step() {
  echo
  echo "== $*"
  "$@"
}

echo "High-quality offline CT-to-FE BV/TV demonstration"
echo "CT_PATH=${CT_PATH}"
echo "OUT_ROOT=${OUT_ROOT}"
echo "ELEMENTS=${ELEMENTS}"
echo "SPHERE_DIAMETER_MM=${SPHERE_DIAMETER_MM}"
echo "IMAGE_DPI=${IMAGE_DPI}"
echo "MAX_ELEMENTS_RENDERED=${MAX_ELEMENTS_RENDERED}"

run_step python3 scripts/make_demo_inp.py \
  --ct "$CT_PATH" \
  --out "$MESH_PATH" \
  --elements "$ELEMENTS" \
  --anchor-threshold "$ANCHOR_THRESHOLD" \
  --padding-voxels "$PADDING_VOXELS"

run_step python3 scripts/map_bvtv_to_fe.py \
  --ct "$CT_PATH" \
  --mesh "$MESH_PATH" \
  --out-prefix "$FIXED_PREFIX" \
  --threshold "$FIXED_THRESHOLD" \
  --sphere-diameter-mm "$SPHERE_DIAMETER_MM"

run_step python3 scripts/map_bvtv_to_fe.py \
  --ct "$CT_PATH" \
  --mesh "$MESH_PATH" \
  --out-prefix "$OTSU_PREFIX" \
  --threshold otsu \
  --sphere-diameter-mm "$SPHERE_DIAMETER_MM"

run_step python3 scripts/map_bvtv_to_fe.py \
  --ct "$CT_PATH" \
  --mesh "$MESH_PATH" \
  --out-prefix "$REGISTERED_PREFIX" \
  --threshold "$FIXED_THRESHOLD" \
  --sphere-diameter-mm "$SPHERE_DIAMETER_MM" \
  --landmarks-json "$LANDMARKS_PATH" \
  --peri-implant-out "$PERI_IMPLANT_OUT" \
  --implant-centerline-start=-12.7902,-71.5878,-31.8359 \
  --implant-centerline-end=-12.7902,-71.5878,-13.9062 \
  --implant-diameter-mm 4.0 \
  --peri-implant-shell-thickness-mm 1.0

run_step python3 scripts/visualize_bvtv.py \
  --ct "$CT_PATH" \
  --mesh "$MESH_PATH" \
  --mapped-inp "${FIXED_PREFIX}.mapped.inp" \
  --results "${FIXED_PREFIX}.csv" \
  --out-dir "${OUT_ROOT}/visualization_threshold_${FIXED_THRESHOLD}" \
  --title "High-quality CT-to-FE BV/TV Report - threshold ${FIXED_THRESHOLD} HU" \
  --image-dpi "$IMAGE_DPI" \
  --max-elements-rendered "$MAX_ELEMENTS_RENDERED"

run_step python3 scripts/visualize_bvtv.py \
  --ct "$CT_PATH" \
  --mesh "$MESH_PATH" \
  --mapped-inp "${OTSU_PREFIX}.mapped.inp" \
  --results "${OTSU_PREFIX}.csv" \
  --out-dir "${OUT_ROOT}/visualization_otsu" \
  --title "High-quality CT-to-FE BV/TV Report - Otsu threshold" \
  --image-dpi "$IMAGE_DPI" \
  --max-elements-rendered "$MAX_ELEMENTS_RENDERED"

run_step python3 tests/test_mapper_synthetic.py

run_step python3 scripts/verify_pipeline.py \
  --ct "$CT_PATH" \
  --mesh "$MESH_PATH" \
  --csv-300 "${FIXED_PREFIX}.csv" \
  --materials-300 "${FIXED_PREFIX}.materials.inp" \
  --mapped-300 "${FIXED_PREFIX}.mapped.inp" \
  --csv-otsu "${OTSU_PREFIX}.csv" \
  --materials-otsu "${OTSU_PREFIX}.materials.inp" \
  --mapped-otsu "${OTSU_PREFIX}.mapped.inp" \
  --csv-registered "${REGISTERED_PREFIX}.csv" \
  --materials-registered "${REGISTERED_PREFIX}.materials.inp" \
  --mapped-registered "${REGISTERED_PREFIX}.mapped.inp" \
  --peri-implant "$PERI_IMPLANT_OUT" \
  --visualization-300 "${OUT_ROOT}/visualization_threshold_${FIXED_THRESHOLD}" \
  --visualization-otsu "${OUT_ROOT}/visualization_otsu" \
  --report "$VERIFY_REPORT"

echo
echo "High-quality offline demonstration completed."
echo "Main fixed-threshold report: ${OUT_ROOT}/visualization_threshold_${FIXED_THRESHOLD}/report.html"
echo "Main Otsu report: ${OUT_ROOT}/visualization_otsu/report.html"
echo "Final mapped Abaqus input: ${FIXED_PREFIX}.mapped.inp"
echo "Verification report: ${VERIFY_REPORT}"
