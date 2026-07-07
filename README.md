# Genome-scale T cell Perturb-seq — hackathon

Pipeline para trabajar con el dataset **Primary Human CD4+ T Cell Perturb-seq (Marson Lab, 2025)**
publicado en la [CZI Virtual Cells Platform](https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq).

Perturb-seq de escala genómica en células T CD4+ primarias humanas: 4 donantes × 3 condiciones
(Rest, Stim8hr, Stim48hr).

## Fuente de datos

Bucket S3 **público** (no requiere credenciales, se accede con `--no-sign-request`):

```
s3://genome-scale-tcell-perturb-seq/marson2025_data/
```

- Código de análisis original: https://github.com/emdann/GWT_perturbseq_analysis_2025
- Preprint: https://www.biorxiv.org/content/10.64898/2025.12.23.696273v1
- Raw / cellranger: SRA `SRP643211` / GEO `GSE314342`

## ⚠️ Tamaño

| Grupo | Contenido | Tamaño aprox. |
|-------|-----------|---------------|
| `tables` | CSV suplementarias | ~15 MB |
| `de` | `GWCD4i.DE_stats.*` (h5ad/h5mu) | ~63 GB |
| `pseudobulk` | `GWCD4i.pseudobulk_merged.h5ad` | ~45 GB |
| `cell` (c/u) | `D*_*.assigned_guide.h5ad` | ~120–170 GB por archivo |
| **all** | **todo** | **~1.8 TB** |

Los 12 archivos cell-level suman ~1.7 TB. **No caben en una laptop.** Para el hackathon
empieza por `tables` y, si necesitas expresión, `pseudobulk`.

## Setup

```bash
# 1. AWS CLI (ya instalado en esta máquina vía Homebrew)
aws --version

# 2. Entorno Python
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

## Descargar datos

```bash
scripts/download.sh ls                 # listar todo el bucket (sin descargar)
scripts/download.sh tables             # solo tablas CSV (~15 MB) -- recomendado para arrancar
scripts/download.sh pseudobulk         # matriz pseudobulk (~45 GB)
scripts/download.sh de                 # resultados de DE (~63 GB)
scripts/download.sh cell D1_Rest       # un archivo cell-level concreto (~140 GB)
scripts/download.sh all                # TODO (~1.8 TB), pide confirmación
```

Los datos se guardan en `./data/` (ignorado por git). Cambia el destino con `DEST=/ruta`.

## Inspeccionar

```bash
python scripts/inspect.py data/GWCD4i.DE_stats.h5ad
```

## Estructura del repo

```
.
├── data/            # datos descargados (git-ignored)
├── metadata/        # notas / metadatos versionables
├── notebooks/       # análisis exploratorio
├── scripts/
│   ├── download.sh  # descarga selectiva desde S3
│   ├── list_data.sh # listar el bucket
│   └── inspect.py   # inspección rápida de h5ad/h5mu
├── requirements.txt
└── README.md
```

## Esquema de los datos

Ver descripción completa de columnas (`.obs`, `.var`, `.layers`) en el readme oficial del bucket:

```bash
aws s3 cp --no-sign-request \
  s3://genome-scale-tcell-perturb-seq/marson2025_data/data_sharing_readme.md -
```
