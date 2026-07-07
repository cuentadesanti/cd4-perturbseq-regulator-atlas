# CD4 Perturb-seq Regulator Atlas — API

A **read-only** API over the pipeline's versioned outputs. *The pipeline produces the artifacts;
the API makes them explorable.*

> Runs **no** models · downloads **no** h5ad · touches **no** S3 · writes **nothing**.
> `make all` remains the source of truth. The API only reads `docs/tables/*.csv` at startup.

## Run

```bash
pip install -r requirements.txt          # includes fastapi + uvicorn (API section)
make all                                  # generates the CSVs the API serves
uvicorn api.main:app --reload --port 8000
```

- **UI** (Regulator Atlas): http://localhost:8000/ — the API serves the frontend on the same origin,
  so there's no file to open and no CORS to configure.
- **Swagger UI** (demo without the frontend): http://localhost:8000/docs

If port 8000 is taken, run on another (`--port 8010`) and open `http://localhost:8010/`.

If optional tables (reproducibility, edges) are missing, the endpoints still work and return
`available: false` / `null` fields.

## Endpoints

| Method | Path | What it returns |
|---|---|---|
| GET | `/health` | status + which tables loaded |
| GET | `/summary` | counts, top regulators, repro/edges availability |
| GET | `/regulators` | filterable list (`q`, `regulator_class`, `limit`, `sort_by`) |
| GET | `/regulators/{gene}` | **full profile** of a regulator |
| GET | `/regulators/{gene}/edges` | robust edges of that regulator (if any) |
| GET | `/classes/{global\|condition-specific}` | regulators by class |
| GET | `/audit/reproducibility` | audit table + coverage |
| GET | `/audit/kd-gate` | raw / KD-gated / EB comparison |
| GET | `/edges/summary` | edges per regulator |
| GET | `/edges/downstream` | most convergent downstream genes |

`sort_by` ∈ `core_rank` · `stability_frequency` · `reweighted_score`.

## Quick demo (Swagger or curl)

```bash
curl localhost:8000/summary
curl localhost:8000/regulators/SGF29          # global, survives all audits
curl localhost:8000/regulators/ZAP70          # context-specific (TCR)
curl "localhost:8000/regulators?regulator_class=global&sort_by=stability_frequency&limit=10"
curl localhost:8000/audit/reproducibility
```

## Design

- `data_loader.py` loads the CSVs into memory at startup and indexes them by gene; a gene missing from
  an optional table → `null` fields.
- `schemas.py` defines the Pydantic models (typed Swagger).
- `main.py` exposes the endpoints with CORS enabled for the local frontend.
- Per-condition profile: derived from `DE_stats.suppl_table.csv` if local (groupby only, no model).
