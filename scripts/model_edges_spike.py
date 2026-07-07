#!/usr/bin/env python3
"""Modelo 1 — SPIKE (experimento controlado, ESTRICTAMENTE OPCIONAL).

Objetivo: comprobar si podemos leer por *slice* filas concretas de log_fc/lfcSE del
h5ad de 17 GB directamente desde S3, SIN descargar el archivo (solo hay ~9.8 GB de disco).

NO asume nada sobre el layout: inspecciona chunks y MIDE bytes y tiempo por fila.
Si las librerías faltan, o la lectura es lenta/frágil, reporta 'INVIABLE' y termina
con éxito (exit 0) — el entregable oficial NO depende de esto.

    pip install h5py s3fs fsspec
    python scripts/model_edges_spike.py
"""
import sys
import time

S3_URL = "s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.DE_stats.h5ad"
SLOW_SECONDS_PER_ROW = 20.0     # umbral: si una fila tarda más, se declara lento
N_ROWS_TEST = 3


def verdict(msg, viable):
    print("\n" + "=" * 64)
    print(f"  VEREDICTO: {'VIABLE' if viable else 'INVIABLE'} — {msg}")
    print("=" * 64)
    if not viable:
        print("  → El entregable oficial es Modelo 2 + documentación (no se ve afectado).")
    sys.exit(0)   # nunca fatal: es opcional


def main():
    try:
        import h5py, fsspec, s3fs  # noqa: F401
    except Exception as e:
        verdict(f"faltan librerías ({type(e).__name__}: {e}). Instala: pip install h5py s3fs fsspec", False)

    print(f"== Modelo 1 · spike de acceso remoto ==\n  fuente: {S3_URL}")
    import h5py, fsspec
    try:
        t0 = time.time()
        f = fsspec.open(S3_URL, anon=True).open()
        h5 = h5py.File(f, "r")
        print(f"  h5ad abierto en {time.time()-t0:.1f}s (streaming, sin descargar)")
    except Exception as e:
        verdict(f"no se pudo abrir el h5ad remoto ({type(e).__name__}: {e})", False)

    # inspección del layout de los layers
    try:
        layers = h5["layers"]
        names = list(layers.keys())
        print(f"  layers disponibles: {names}")
        need = [n for n in ("log_fc", "lfcSE") if n in names]
        if len(need) < 2:
            verdict(f"faltan layers log_fc/lfcSE (hay {names})", False)
        for n in need:
            d = layers[n]
            print(f"    {n:8s} shape={d.shape} dtype={d.dtype} chunks={d.chunks} "
                  f"compression={d.compression}")
        d0 = layers[need[0]]
        n_obs, n_vars = d0.shape
        row_bytes = n_vars * d0.dtype.itemsize
        print(f"  matriz: {n_obs:,} obs × {n_vars:,} genes · "
              f"~{row_bytes/1e6:.2f} MB por fila por layer")
        if d0.chunks and d0.chunks[0] > 64:
            print(f"  [aviso] chunk de {d0.chunks[0]} filas: leer 1 fila puede arrastrar el bloque entero")
    except Exception as e:
        verdict(f"no se pudo inspeccionar los layers ({type(e).__name__}: {e})", False)

    # medir lectura por slice de N filas
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
            print(f"  fila {i:>6}: leída en {dt:5.1f}s · "
                  f"log_fc[:3]={a[:3].round(3)} lfcSE[:3]={b[:3].round(3)}")
        avg = sum(per_row) / len(per_row)
        total_mb = len(idx) * 2 * row_bytes / 1e6
        print(f"  media {avg:.1f}s/fila · ~{total_mb:.1f} MB transferidos en total")
        h5.close()
    except Exception as e:
        verdict(f"la lectura por slice falló ({type(e).__name__}: {e})", False)

    if avg > SLOW_SECONDS_PER_ROW:
        verdict(f"lectura lenta ({avg:.1f}s/fila > {SLOW_SECONDS_PER_ROW}s). "
                f"Escalar a ~150 filas sería impráctico.", False)
    else:
        est = avg * 150
        verdict(f"lectura por slice OK ({avg:.1f}s/fila). "
                f"Escalar a 150 reguladores ≈ {est:.0f}s. Se puede correr model_edges.py.", True)


if __name__ == "__main__":
    main()
