#!/usr/bin/env bash
#
# Descarga selectiva del dataset "Genome-scale T cell Perturb-seq" (Marson 2025)
# Fuente: bucket público de la CZI Virtual Cells Platform (no requiere credenciales).
#
# ADVERTENCIA DE TAMAÑO:
#   - Los 12 archivos cell-level (D*_*.assigned_guide.h5ad) pesan ~120-170 GB CADA UNO
#     (~1.7 TB en total). NO caben en la mayoría de laptops.
#   - Para empezar/hackear usa el grupo "light" (tablas CSV + DE stats + pseudobulk),
#     que es de unos pocos GB.
#
# Uso:
#   scripts/download.sh light                 # tablas + DE stats + pseudobulk (~90 GB)
#   scripts/download.sh tables                # solo tablas CSV suplementarias (~15 MB)
#   scripts/download.sh de                    # resultados de DE (h5ad/h5mu) (~63 GB)
#   scripts/download.sh pseudobulk            # pseudobulk_merged.h5ad (~45 GB)
#   scripts/download.sh cell D1_Rest          # un archivo cell-level concreto (~140 GB)
#   scripts/download.sh all                   # TODO (~1.8 TB) -- pide confirmación
#   scripts/download.sh ls                    # solo listar, no descargar
#
# Destino: ./data/  (configurable con DEST=...)
set -euo pipefail

BUCKET="s3://genome-scale-tcell-perturb-seq/marson2025_data"
DEST="${DEST:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/data}"
AWS=(aws s3 cp --no-sign-request)
AWS_SYNC=(aws s3 sync --no-sign-request)

if ! command -v aws >/dev/null 2>&1; then
  echo "ERROR: 'aws' no está instalado. Instálalo con:  brew install awscli" >&2
  exit 1
fi

mkdir -p "${DEST}"
CMD="${1:-help}"

dl_tables() {
  echo ">> Descargando tablas suplementarias (CSV) a ${DEST}/suppl_tables/"
  "${AWS_SYNC[@]}" "${BUCKET}/suppl_tables/" "${DEST}/suppl_tables/"
}

dl_de() {
  echo ">> Descargando resultados de DE a ${DEST}/"
  "${AWS[@]}" "${BUCKET}/GWCD4i.DE_stats.h5ad"           "${DEST}/"
  "${AWS[@]}" "${BUCKET}/GWCD4i.DE_stats.by_guide.h5mu"  "${DEST}/"
  "${AWS[@]}" "${BUCKET}/GWCD4i.DE_stats.by_donors.h5mu" "${DEST}/"
}

dl_pseudobulk() {
  echo ">> Descargando pseudobulk a ${DEST}/"
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
      echo "Uso: $0 cell <D1_Rest|D1_Stim8hr|D1_Stim48hr|D2_...|D4_Stim48hr>" >&2
      exit 1
    fi
    KEY="${BUCKET}/${NAME}.assigned_guide.h5ad"
    echo ">> Descargando archivo cell-level (~120-170 GB): ${KEY}"
    "${AWS[@]}" "${KEY}" "${DEST}/"
    ;;
  all)
    read -r -p "Vas a descargar ~1.8 TB. ¿Seguro? (yes/no) " ans
    [[ "${ans}" == "yes" ]] || { echo "Cancelado."; exit 0; }
    "${AWS_SYNC[@]}" "${BUCKET}/" "${DEST}/"
    ;;
  help|*)
    sed -n '2,26p' "${BASH_SOURCE[0]}"
    ;;
esac

echo "Listo."
