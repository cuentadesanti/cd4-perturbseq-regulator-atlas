PY := python3

.PHONY: all eda model audit report spike edges eda-edges repro-meta fingerprints operator operator-tensor operator-svd operator-cp operator-completion operator-donors operator-deconv spectral class-programs specificity-control disease-overlap module-gwas convergence-extras convergence-figures api clean pipeline help

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
	@echo "  make operator-completion - Step 3: out-of-panel condition extrapolation (flagship) + entry-wise sanity"
	@echo "  make operator-deconv - Step 5 STRETCH: square-block deconvolution + asymmetric subsumption (hypothesis-generating only)"
	@echo "  make operator - empirical regulatory operator (z-score): tensor + SVD + CP + completion [+donors if fetched]; see docs/OPERATOR_ANALYSIS.md"
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

# Step 0: regulator x gene x condition z-score tensor over the expanded ~800-reg panel.
# Fail-closed: refuses to cache unless the confound guard passes (|spearman(||slab||,n_cells)|<0.15).
operator-tensor:
	$(PY) scripts/build_operator_tensor.py --n-total 800 --top-genes 2000

# Step 2 CP. NB: these are the FULL defaults (sweep-to-8, 100 bootstraps) and take ~20+ min;
# the committed operator_cp_* tables were produced by a reduced local pilot
# (--max-rank 6 --stab-subsample 250 --boot-n 50). Use the full config on cluster compute.
operator-cp:
	$(PY) scripts/decompose_operator_cp.py --rank auto --scale-control rms

# Step 1: gene programs = right singular vectors of the operator matrix (z-score),
# with an ISG-anchor orientation, optional varimax rotation, offline enrichment,
# and a power gate on the left factors (needs `make operator-tensor` first).
operator-svd:
	$(PY) scripts/decompose_operator_svd.py --k 10 --tail-pct 2 --rotate

# Step 3: is the operator recoverably low-rank? 3b (flagship) holds out (regulator,
# Stim48hr) fibers preferentially for out-of-panel regulators and predicts them from
# Rest+Stim8hr via low-rank soft-impute, vs a persistence baseline (Stim48hr:=Stim8hr).
# 3a (sanity) is random entry-wise completion (elbow only; needs `make operator-tensor` first).
operator-completion:
	$(PY) scripts/operator_completion.py --max-rank 12 --holdout 0.2

# Step 4: are the gene programs donor-reproducible AS SUBSPACES?
# Principal angles between top-k gene-program subspaces of DISJOINT donor pairs only.
operator-donors:
	$(PY) scripts/operator_donor_angles.py --k 5

# Step 5 STRETCH: square-block deconvolution A = I - L^{-1} (regularized) + asymmetric
# subsumption. HYPOTHESIS-GENERATING ONLY: valid only on the square block (regulators
# that are also readout genes), linear approx to a nonlinear system, every edge carries
# the block condition number (needs `make operator-tensor` first).
operator-deconv:
	$(PY) scripts/operator_deconvolution.py --n-robust 50 --ridge 1e-2

# empirical regulatory operator umbrella: Step 0 fetch is one-time then cached; Steps 1-3 local;
# operator-donors prints NEEDS-DATA unless per-donor matrices are fetched.
operator: operator-tensor operator-svd operator-cp operator-completion operator-donors

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
