PY := python3

.PHONY: all eda model audit report spike edges eda-edges repro-meta fingerprints api clean pipeline help

help:
	@echo "Targets:"
	@echo "  make eda      - EDA 80/20 (figuras + tablas)          [solo CSV local]"
	@echo "  make model    - Modelo 2 EB: ranking de reguladores    [solo CSV local]"
	@echo "  make audit    - Regulator ranking audit (baselines, estabilidad, clases)"
	@echo "  make report   - reporte consolidado + figura overview"
	@echo "  make all      - pipeline completo verificado (run_pipeline.py)"
	@echo "  make spike    - Modelo 1: spike de acceso remoto (OPCIONAL)"
	@echo "  make edges    - Modelo 1: red de edges, bonus (OPCIONAL, remoto)"
	@echo "  make eda-edges- EDA de la red de edges (usa robust_edges.csv si existe)"
	@echo "  make repro-meta - extrae .obs de reproducibilidad del h5ad (OPCIONAL, remoto)"
	@echo "  make fingerprints - side analysis: mapa de programas (PCA/similitud/clusters) [zscore remoto o log_fc cache]"
	@echo "  make api      - Regulator Atlas API read-only (uvicorn :8000, Swagger en /docs)"
	@echo "  make clean    - borra outputs generados (docs/figures, docs/tables, report)"

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

# side analysis (NO entra en `make all`): transcriptional fingerprints.
# --matrix zscore lee slice remoto; --matrix log_fc usa el cache local si existe.
fingerprints:
	$(PY) scripts/analyze_fingerprints.py --n 200 --matrix zscore --top-genes 2000

api:
	$(PY) -m uvicorn api.main:app --reload --port 8000

clean:
	rm -f docs/figures/*.png docs/tables/*.csv docs/report.md
	@echo "outputs generados eliminados (los CSV de entrada en data/ NO se tocan)"
