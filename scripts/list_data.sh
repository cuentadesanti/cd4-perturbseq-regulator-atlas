#!/usr/bin/env bash
# Lista el contenido del bucket público del dataset (sin credenciales).
set -euo pipefail

BUCKET="s3://genome-scale-tcell-perturb-seq/marson2025_data/"

echo "== Contenido de ${BUCKET} =="
aws s3 ls --no-sign-request --recursive --human-readable "${BUCKET}"
