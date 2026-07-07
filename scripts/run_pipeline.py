#!/usr/bin/env python3
"""Pipeline reproducible del core (solo CSV locales, sin descargas ni h5ad).

    python scripts/run_pipeline.py

1. verifica que existan los CSV de entrada
2. corre EDA → Modelo 2 (EB) → reporte
3. verifica que existan los outputs esperados
4. imprime un resumen conciso de éxito
"""
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PY = sys.executable

INPUTS = [
    "data/suppl_tables/DE_stats.suppl_table.csv",
    "data/suppl_tables/sgrna_library_metadata.suppl_table.csv",
    "data/suppl_tables/sample_metadata.suppl_table.csv",
]
STEPS = [
    ("EDA 80/20", "scripts/eda.py"),
    ("Modelo 2 · EB", "scripts/model_hubs.py"),
    ("Ranking audit", "scripts/audit_ranking.py"),
    ("Reporte", "scripts/build_report.py"),
]
OUTPUTS = [
    "docs/tables/hub_ranking_bayes.csv",
    "docs/tables/top_regulators_for_review.csv",
    "docs/tables/top_robust_regulators.csv",
    "docs/tables/ranking_baseline_comparison.csv",
    "docs/tables/hub_ranking_stability.csv",
    "docs/tables/top_global_regulators.csv",
    "docs/tables/top_condition_specific_regulators.csv",
    "docs/tables/regulator_classes.csv",
    "docs/figures/07_hub_posterior_ranking.png",
    "docs/figures/08_kd_gate_changes_ranking.png",
    "docs/figures/10_global_vs_context_specific.png",
    "docs/figures/11_ranking_stability.png",
    "docs/figures/00_pipeline_overview.png",
    "docs/report.md",
]


def fail(msg):
    print(f"\n✗ FALLO: {msg}")
    sys.exit(1)


def main():
    print("=" * 60 + "\n  PIPELINE — core reproducible (solo CSV local)\n" + "=" * 60)

    print("\n[1/4] Verificando inputs…")
    missing = [p for p in INPUTS if not (ROOT / p).exists()]
    if missing:
        fail("faltan CSV de entrada:\n   " + "\n   ".join(missing) +
             "\n   → descárgalos con: scripts/download.sh tables")
    for p in INPUTS:
        print(f"   ✓ {p}")

    print("\n[2/4] Ejecutando pasos…")
    t0 = time.time()
    for name, script in STEPS:
        print(f"   → {name} ({script})")
        r = subprocess.run([PY, str(ROOT / script)], cwd=ROOT,
                           capture_output=True, text=True)
        if r.returncode != 0:
            print(r.stdout[-1500:]); print(r.stderr[-1500:])
            fail(f"{script} devolvió código {r.returncode}")
    elapsed = time.time() - t0

    print("\n[3/4] Verificando outputs…")
    missing = [p for p in OUTPUTS if not (ROOT / p).exists()]
    if missing:
        fail("faltan outputs esperados:\n   " + "\n   ".join(missing))
    for p in OUTPUTS:
        kb = (ROOT / p).stat().st_size / 1024
        print(f"   ✓ {p:<48} {kb:7.1f} KB")

    print("\n[4/4] Resumen")
    import pandas as pd
    rank = pd.read_csv(ROOT / "docs/tables/hub_ranking_bayes.csv")
    top3 = ", ".join(rank.head(3)["target_contrast_gene_name"])
    print(f"   reguladores rankeados : {len(rank):,}")
    print(f"   top 3 (EB)            : {top3}")
    edges = ROOT / "docs/tables/robust_edges.csv"
    if edges.exists():
        import pandas as pd
        print(f"   edges bonus (Modelo 1): {len(pd.read_csv(edges)):,}  (docs/tables/robust_edges.csv)")
    else:
        print("   edges bonus (Modelo 1): no generados (opcional · make edges)")
    print(f"\n✓ PIPELINE OK en {elapsed:.0f}s — reporte en docs/report.md")


if __name__ == "__main__":
    main()
