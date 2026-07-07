#!/usr/bin/env bash
#
# Selective download of the "Genome-scale T cell Perturb-seq" dataset (Marson 2025)
# Source: public CZI Virtual Cells Platform bucket (no credentials required).
#
# SIZE WARNING:
#   - The 12 cell-level files (D*_*.assigned_guide.h5ad) are ~120-170 GB EACH
#     (~1.7 TB total). They do NOT fit on most laptops.
#   - To get started, use the "light" group (CSV tables + DE stats + pseudobulk),
#     which is a few GB.
#
# Usage:
#   scripts/download.sh light                 # tables + DE stats + pseudobulk (~90 GB)
#   scripts/download.sh tables                # supplementary CSV tables only (~15 MB)
#   scripts/download.sh de                    # DE results (h5ad/h5mu) (~63 GB)
#   scripts/download.sh pseudobulk            # pseudobulk_merged.h5ad (~45 GB)
#   scripts/download.sh cell D1_Rest          # one specific cell-level file (~140 GB)
#   scripts/download.sh all                   # EVERYTHING (~1.8 TB) -- asks for confirmation
#   scripts/download.sh ls                    # list only, do not download
#
# Destination: ./data/  (configurable with DEST=...)
set -euo pipefail

BUCKET="s3://genome-scale-tcell-perturb-seq/marson2025_data"
DEST="${DEST:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/data}"
AWS=(aws s3 cp --no-sign-request)
AWS_SYNC=(aws s3 sync --no-sign-request)

if ! command -v aws >/dev/null 2>&1; then
  echo "ERROR: 'aws' is not installed. Install it with:  brew install awscli" >&2
  exit 1
fi

mkdir -p "${DEST}"
CMD="${1:-help}"

dl_tables() {
  echo ">> Downloading supplementary tables (CSV) to ${DEST}/suppl_tables/"
  "${AWS_SYNC[@]}" "${BUCKET}/suppl_tables/" "${DEST}/suppl_tables/"
}

dl_de() {
  echo ">> Downloading DE results to ${DEST}/"
  "${AWS[@]}" "${BUCKET}/GWCD4i.DE_stats.h5ad"           "${DEST}/"
  "${AWS[@]}" "${BUCKET}/GWCD4i.DE_stats.by_guide.h5mu"  "${DEST}/"
  "${AWS[@]}" "${BUCKET}/GWCD4i.DE_stats.by_donors.h5mu" "${DEST}/"
}

dl_pseudobulk() {
  echo ">> Downloading pseudobulk to ${DEST}/"
  "${AWS[@]}" "${BUCKET}/GWCD4i.pseudobulk_merged.h5ad" "${DEST}/"
}

case "${CMD}" in
  ls)
    aws s3 ls --no-sign-request --recursive --human-readable "${BUCKET}/"
    ;;
  tables)
    dl_tables
    ;;
  de)
    dl_de
    ;;
  pseudobulk)
    dl_pseudobulk
    ;;
  light)
    dl_tables; dl_pseudobulk; dl_de
    ;;
  cell)
    NAME="${2:-}"
    if [[ -z "${NAME}" ]]; then
      echo "Usage: $0 cell <D1_Rest|D1_Stim8hr|D1_Stim48hr|D2_...|D4_Stim48hr>" >&2
      exit 1
    fi
    KEY="${BUCKET}/${NAME}.assigned_guide.h5ad"
    echo ">> Downloading cell-level file (~120-170 GB): ${KEY}"
    "${AWS[@]}" "${KEY}" "${DEST}/"
    ;;
  all)
    read -r -p "You are about to download ~1.8 TB. Are you sure? (yes/no) " ans
    [[ "${ans}" == "yes" ]] || { echo "Cancelled."; exit 0; }
    "${AWS_SYNC[@]}" "${BUCKET}/" "${DEST}/"
    ;;
  help|*)
    sed -n '2,26p' "${BASH_SOURCE[0]}"
    ;;
esac

echo "Done."
