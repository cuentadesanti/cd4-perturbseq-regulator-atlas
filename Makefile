PY := python3

.PHONY: all eda model audit report spike edges eda-edges repro-meta fingerprints spectral class-programs specificity-control disease-overlap module-gwas convergence-extras convergence-figures api clean pipeline help

help:
	@echo "Targets:"
	@echo "  make eda      - 80/20 EDA (figures + tables)            [local CSV only]"
	@echo "  make model    - Model 2 EB: regulator ranking           [local CSV only]"
	@echo "  make audit    - regulator ranking audit (baselines, stability, classes)"
	@echo "  make report   - consolidated report + overview figure"
	@echo "  make all      - full verified pipeline (run_pipeline.py)"
	@echo "  make spike    - Model 1: remote-access spike (OPTIONAL)"
	@echo "  make edges    - Model 1: effect network, bonus (OPTIONAL, remote)"
	@echo "  make eda-edges- EDA of the effect network (uses robust_edges.csv if present)"
	@echo "  make repro-meta - extract reproducibility .obs from the h5ad (OPTIONAL, remote)"
	@echo "  make fingerprints - transcriptional programs (PCA/similarity/clusters) [remote zscore or log_fc cache]"
	@echo "  make spectral - spectral sanity check on the program assignments (after fingerprints)"
	@echo "  make class-programs - balanced 30-regulator panel: distinct classes → distinct programs? (after fingerprints)"
	@echo "  make convergence-extras - specificity control + disease overlap (fully offline)"
	@echo "  make module-gwas - autoimmune GWAS overlap via Open Targets (NEEDS NETWORK; committed table is the record)"
	@echo "  make convergence-figures - re-render figs 26 + 28 from tables via the shared palette (offline)"
	@echo "  make api      - Regulator Atlas read-only API (uvicorn :8000, Swagger at /docs)"
	@echo "  make clean    - remove generated outputs (docs/figures, docs/tables, report)"

eda:
	$(PY) scripts/eda.py

model:
	$(PY) scripts/model_hubs.py

audit:
	$(PY) scripts/audit_ranking.py

report:
	$(PY) scripts/build_report.py

all: pipeline

pipeline:
	$(PY) scripts/run_pipeline.py

# --- opcional / remoto (requiere: pip install h5py s3fs fsspec) ---
spike:
	$(PY) scripts/model_edges_spike.py

edges:
	$(PY) scripts/model_edges.py --n 8

eda-edges:
	$(PY) scripts/eda_edges.py

repro-meta:
	$(PY) scripts/extract_de_obs_metadata.py

# transcriptional programs (NOT part of `make all`): fingerprint analysis.
# --matrix zscore reads a remote slice; --matrix log_fc uses the local cache if present.
fingerprints:
	$(PY) scripts/analyze_fingerprints.py --n 200 --matrix zscore --top-genes 2000

# spectral sanity check on the program assignments (needs `make fingerprints` first)
spectral:
	$(PY) scripts/spectral_sanity_check.py

# balanced 30-regulator class programs: do distinct classes hit distinct programs?
# NEEDS `make fingerprints` first (uses the cached panel + log_fc/zscore; fully offline).
class-programs:
	$(PY) scripts/analyze_class_programs.py

# convergence extras (docs/disease_and_specificity.md).
# convergence-extras is FULLY OFFLINE (specificity-control + disease-overlap).
# module-gwas is SEPARATE because it needs network (Open Targets); the committed
# module_gwas_hits.csv is the record and is NOT offline-reproducible.
specificity-control:
	$(PY) scripts/analyze_chromatin_stress_control.py
disease-overlap:
	$(PY) scripts/analyze_disease_overlap.py
convergence-extras: specificity-control disease-overlap
module-gwas:   # needs network (Open Targets) — see note above
	$(PY) scripts/analyze_module_gwas.py

# re-render the two precomputed convergence figures (26 module, 28 phase-2) from the committed
# tables + log_fc cache, threaded through the shared palette (scripts/_figstyle.py). Fully offline.
convergence-figures:
	$(PY) scripts/render_convergence_figures.py

api:
	$(PY) -m uvicorn api.main:app --reload --port 8000

clean:
	rm -f docs/figures/*.png docs/tables/*.csv docs/report.md
	@echo "generated outputs removed (input CSVs in data/ are NOT touched)"
