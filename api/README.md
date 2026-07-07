# CD4 Perturb-seq Regulator Atlas — API

API **read-only** sobre los outputs versionados del pipeline. *El pipeline produce los artefactos;
la API los hace explorables.*

> **No** corre modelos · **no** descarga h5ad · **no** toca S3 · **no** escribe nada.
> `make all` sigue siendo la fuente de verdad. La API solo lee `docs/tables/*.csv` al arrancar.

## Correr

```bash
pip install -r requirements.txt          # incluye fastapi + uvicorn (sección API)
make all                                  # genera los CSV que la API sirve
uvicorn api.main:app --reload --port 8000
```

- **Swagger UI** (demo sin frontend): http://localhost:8000/docs
- Frontend mínimo: abrir `frontend/index.html` (apunta a `localhost:8000` por defecto).

Si faltan tablas opcionales (reproducibilidad, edges), los endpoints siguen funcionando y devuelven
`available: false` / campos `null`.

## Endpoints

| Método | Ruta | Qué devuelve |
|---|---|---|
| GET | `/health` | estado + qué tablas se cargaron |
| GET | `/summary` | conteos, top reguladores, disponibilidad de repro/edges |
| GET | `/regulators` | lista filtrable (`q`, `regulator_class`, `limit`, `sort_by`) |
| GET | `/regulators/{gene}` | **perfil completo** de un regulador |
| GET | `/regulators/{gene}/edges` | edges robustos de ese regulador (si existen) |
| GET | `/classes/{global\|context-specific}` | reguladores por clase |
| GET | `/audit/reproducibility` | tabla de auditoría + cobertura |
| GET | `/audit/kd-gate` | comparación raw / KD-gated / EB |
| GET | `/edges/summary` | edges por regulador |
| GET | `/edges/downstream` | genes downstream más convergentes |

`sort_by` ∈ `core_rank` · `stability_frequency` · `reweighted_score`.

## Demo rápida (Swagger o curl)

```bash
curl localhost:8000/summary
curl localhost:8000/regulators/SGF29          # global, sobrevive todas las auditorías
curl localhost:8000/regulators/ZAP70          # context-specific (TCR)
curl "localhost:8000/regulators?regulator_class=global&sort_by=stability_frequency&limit=10"
curl localhost:8000/audit/reproducibility
```

## Diseño

- `data_loader.py` carga los CSV en memoria al startup e indexa por gen; un gen faltante en una tabla
  opcional → campos `null`.
- `schemas.py` define los modelos Pydantic (Swagger tipado).
- `main.py` expone los endpoints con CORS habilitado para el frontend local.
- Perfil por condición: derivado de `DE_stats.suppl_table.csv` si está local (solo groupby, sin modelo).
