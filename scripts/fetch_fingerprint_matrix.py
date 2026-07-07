#!/usr/bin/env python3
"""Descarga y cachea la matriz de huellas transcriptómicas (layers/log_fc).

El h5ad remoto (17 GB) guarda `layers/log_fc` como matriz DENSA, CONTIGUA y SIN
comprimir de shape (33983 perturbación×condición) x (10282 genes medidos). Cada fila
es la *huella* de una perturbación: el log fold-change de cada gen medido downstream.

Como es contigua, una lectura secuencial rinde ~10 MB/s -> la matriz completa
(~2.8 GB en float64) baja en ~5 min. La cacheamos en float32 (1.4 GB) para que el
análisis downstream (similitud / PCA / clustering) corra instantáneo y offline.

    pip install h5py s3fs fsspec
    python scripts/fetch_fingerprint_matrix.py

Salidas (en la copia principal, data/ está gitignored):
    data/cache/log_fc.f32.npy        (33983, 10282) float32
    data/cache/fingerprint_obs.csv   metadata por fila (perturbación, condición, QC)
    data/cache/fingerprint_var.csv   nombres de genes medidos (columnas)
"""
import time
from pathlib import Path
import numpy as np
import pandas as pd

S3_URL = "s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.DE_stats.h5ad"
# data/ está gitignored y vive en la copia principal, no en el worktree
CACHE = Path("/Users/cuentadesanti/code/hackaton/data/cache")
BLOCK = 2000  # filas por lectura secuencial

# columnas obs útiles para el análisis (QC, reproducibilidad, conteos)
OBS_COLS = [
    "index", "target_contrast_gene_name", "culture_condition", "target_contrast",
    "n_cells_target", "n_up_genes", "n_down_genes", "n_total_de_genes", "n_downstream",
    "ontarget_effect_size", "ontarget_significant", "offtarget_flag", "distal_offtarget_flag",
    "low_target_gex", "neighboring_gene_KD", "single_guide_estimate", "n_guides",
    "guide_correlation_all", "guide_correlation_signif",
    "donor_correlation_all_mean", "donor_correlation_hits_mean",
]


def read_col(grp, name, h5py):
    node = grp[name]
    if isinstance(node, h5py.Group):  # categórica AnnData
        cats = [c.decode() if isinstance(c, bytes) else c for c in node["categories"][:]]
        codes = node["codes"][:]
        return np.array([cats[i] if i >= 0 else None for i in codes], dtype=object)
    arr = node[:]
    if arr.dtype.kind in ("S", "O"):
        return np.array([x.decode() if isinstance(x, bytes) else x for x in arr], dtype=object)
    return arr


def main():
    import h5py, fsspec
    CACHE.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    f = fsspec.open(S3_URL, anon=True, default_cache_type="readahead").open()
    h5 = h5py.File(f, "r")
    print(f"[{time.time()-t0:5.1f}s] h5ad abierto (streaming)")

    # --- obs (metadata por fila) ---
    obs = h5["obs"]
    avail = [c for c in OBS_COLS if c in obs]
    obs_df = pd.DataFrame({c: read_col(obs, c, h5py) for c in avail})
    obs_df.to_csv(CACHE / "fingerprint_obs.csv", index=False)
    print(f"[{time.time()-t0:5.1f}s] obs -> {obs_df.shape} guardado")

    # --- var (genes medidos = columnas) ---
    var = h5["var"]
    vkey = "gene_name" if "gene_name" in var else "_index"
    gene_names = read_col(var, vkey, h5py)
    gene_ids = read_col(var, "gene_ids", h5py) if "gene_ids" in var else gene_names
    pd.DataFrame({"gene_name": gene_names, "gene_id": gene_ids}).to_csv(
        CACHE / "fingerprint_var.csv", index=False)
    print(f"[{time.time()-t0:5.1f}s] var -> {len(gene_names)} genes guardado")

    # --- log_fc (matriz densa, lectura secuencial por bloques) ---
    lf = h5["layers"]["log_fc"]
    n, m = lf.shape
    out = np.empty((n, m), dtype=np.float32)
    for i in range(0, n, BLOCK):
        j = min(i + BLOCK, n)
        out[i:j] = lf[i:j, :].astype(np.float32)
        mb = out[:j].nbytes / 1e6
        print(f"[{time.time()-t0:5.1f}s]   filas {j:5d}/{n}  ({mb:6.0f} MB float32)")
    np.save(CACHE / "log_fc.f32.npy", out)
    print(f"[{time.time()-t0:5.1f}s] log_fc -> {out.shape} float32 guardado "
          f"({out.nbytes/1e9:.2f} GB)")
    print(f"OK · total {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
