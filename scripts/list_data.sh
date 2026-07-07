#!/usr/bin/env bash
# Lists the contents of the dataset's public bucket (no credentials).
set -euo pipefail

BUCKET="s3://genome-scale-tcell-perturb-seq/marson2025_data/"

echo "== Contents of ${BUCKET} =="
aws s3 ls --no-sign-request --recursive --human-readable "${BUCKET}"
