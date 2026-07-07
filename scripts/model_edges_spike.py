#!/usr/bin/env python3
"""Model 1 — SPIKE (controlled experiment, STRICTLY OPTIONAL).

Goal: check whether we can read specific rows of log_fc/lfcSE from the 17 GB h5ad
by *slice* directly from S3, WITHOUT downloading the file (only ~9.8 GB of disk).

It assumes nothing about the layout: it inspects chunks and MEASURES bytes and time
per row. If the libraries are missing, or the read is slow/fragile, it reports
'NOT VIABLE' and exits successfully (exit 0) — the official deliverable does NOT depend on this.

    pip install h5py s3fs fsspec
    python scripts/model_edges_spike.py
"""
import sys
import time

S3_URL = "s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.DE_stats.h5ad"
SLOW_SECONDS_PER_ROW = 20.0     # threshold: if a row takes longer, it is declared slow
N_ROWS_TEST = 3


def verdict(msg, viable):
    print("\n" + "=" * 64)
    print(f"  VERDICT: {'VIABLE' if viable else 'NOT VIABLE'} — {msg}")
    print("=" * 64)
    if not viable:
        print("  → The official deliverable is Model 2 + documentation (unaffected).")
    sys.exit(0)   # never fatal: it is optional


def main():
    try:
        import h5py, fsspec, s3fs  # noqa: F401
    except Exception as e:
        verdict(f"missing libraries ({type(e).__name__}: {e}). Install: pip install h5py s3fs fsspec", False)

    print(f"== Model 1 · remote-access spike ==\n  source: {S3_URL}")
    import h5py, fsspec
    try:
        t0 = time.time()
        f = fsspec.open(S3_URL, anon=True).open()
        h5 = h5py.File(f, "r")
        print(f"  h5ad opened in {time.time()-t0:.1f}s (streaming, no download)")
    except Exception as e:
        verdict(f"could not open the remote h5ad ({type(e).__name__}: {e})", False)

    # inspect the layers layout
    try:
        layers = h5["layers"]
        names = list(layers.keys())
        print(f"  available layers: {names}")
        need = [n for n in ("log_fc", "lfcSE") if n in names]
        if len(need) < 2:
            verdict(f"missing log_fc/lfcSE layers (have {names})", False)
        for n in need:
            d = layers[n]
            print(f"    {n:8s} shape={d.shape} dtype={d.dtype} chunks={d.chunks} "
                  f"compression={d.compression}")
        d0 = layers[need[0]]
        n_obs, n_vars = d0.shape
        row_bytes = n_vars * d0.dtype.itemsize
        print(f"  matrix: {n_obs:,} obs × {n_vars:,} genes · "
              f"~{row_bytes/1e6:.2f} MB per row per layer")
        if d0.chunks and d0.chunks[0] > 64:
            print(f"  [warning] chunk of {d0.chunks[0]} rows: reading 1 row may drag in the whole block")
    except Exception as e:
        verdict(f"could not inspect the layers ({type(e).__name__}: {e})", False)

    # measure per-slice read of N rows
    try:
        idx = sorted(set(int(i) for i in
                     [0, n_obs // 2, n_obs - 1][:N_ROWS_TEST]))
        per_row = []
        for i in idx:
            t = time.time()
            a = layers[need[0]][i, :]
            b = layers[need[1]][i, :]
            dt = time.time() - t
            per_row.append(dt)
            print(f"  row {i:>6}: read in {dt:5.1f}s · "
                  f"log_fc[:3]={a[:3].round(3)} lfcSE[:3]={b[:3].round(3)}")
        avg = sum(per_row) / len(per_row)
        total_mb = len(idx) * 2 * row_bytes / 1e6
        print(f"  mean {avg:.1f}s/row · ~{total_mb:.1f} MB transferred total")
        h5.close()
    except Exception as e:
        verdict(f"the per-slice read failed ({type(e).__name__}: {e})", False)

    if avg > SLOW_SECONDS_PER_ROW:
        verdict(f"slow read ({avg:.1f}s/row > {SLOW_SECONDS_PER_ROW}s). "
                f"Scaling to ~150 rows would be impractical.", False)
    else:
        est = avg * 150
        verdict(f"per-slice read OK ({avg:.1f}s/row). "
                f"Scaling to 150 regulators ≈ {est:.0f}s. model_edges.py can be run.", True)


if __name__ == "__main__":
    main()
