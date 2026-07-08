"""CD4 Perturb-seq Regulator Atlas — read-only API over the `make all` outputs.

The pipeline produces the artifacts; this API only makes them explorable. It runs no models,
downloads no h5ad, touches no S3. `make all` remains the source of truth.

    uvicorn api.main:app --reload --port 8000
    → http://localhost:8000/docs   (Swagger)
"""
from __future__ import annotations
from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse

from api.data_loader import get_store, ROOT
from api.schemas import Health, Summary, RegulatorSummary, RegulatorProfile

FRONTEND = ROOT / "frontend" / "index.html"

app = FastAPI(
    title="CD4 Perturb-seq Regulator Atlas",
    description="Search, rank, audit and explore robust CD4+ T cell regulators. "
                "Read-only API over versioned pipeline outputs.",
    version="1.0.0",
)
app.add_middleware(
    CORSMiddleware, allow_origins=["*"], allow_methods=["GET"], allow_headers=["*"],
)


@app.get("/", include_in_schema=False)
def index():
    """Serves the UI from the same origin as the API (no CORS/file:// issues)."""
    if FRONTEND.exists():
        return FileResponse(FRONTEND)
    return {"app": "CD4 Perturb-seq Regulator Atlas", "docs": "/docs"}


@app.get("/health", response_model=Health, tags=["meta"])
def health():
    s = get_store()
    return {"status": "ok", "n_regulators": len(s.ranking), "tables_loaded": s.loaded}


@app.get("/summary", response_model=Summary, tags=["meta"])
def summary():
    return get_store().summary()


@app.get("/genes", tags=["meta"])
def genes():
    """All ranked gene symbols (for autocomplete / fuzzy search in the UI)."""
    return get_store().ranking["gene"].tolist()


@app.get("/regulators", response_model=list[RegulatorSummary], tags=["regulators"])
def regulators(
    q: str | None = Query(None, description="search by gene symbol"),
    regulator_class: str | None = Query(None, enum=["global", "condition-specific"]),
    limit: int = Query(50, ge=1, le=2000),
    sort_by: str = Query("core_rank", enum=["core_rank", "stability_frequency", "reweighted_score"]),
):
    return get_store().list_regulators(q=q, regulator_class=regulator_class,
                                       limit=limit, sort_by=sort_by)


@app.get("/regulators/{gene}", response_model=RegulatorProfile, tags=["regulators"])
def regulator(gene: str):
    p = get_store().profile(gene.upper())
    if p is None:
        raise HTTPException(404, f"regulator '{gene}' not found (did it pass the KD gate?)")
    return p


@app.get("/regulators/{gene}/edges", tags=["regulators"])
def regulator_edges(gene: str, limit: int = Query(200, ge=1, le=5000)):
    s = get_store()
    if s.edges is None:
        return {"gene": gene.upper(), "edges_available": False, "edges": []}
    edges = s.edges_for(gene.upper(), top=limit)
    return {"gene": gene.upper(), "edges_available": True, "n": len(edges), "edges": edges}


@app.get("/audit/reproducibility", tags=["audit"])
def audit_reproducibility():
    s = get_store()
    from api.data_loader import _records
    if s.audit is None:
        return {"available": False, "coverage": None, "rows": []}
    cov = _records(s.coverage)[0] if s.coverage is not None else None
    return {"available": True, "coverage": cov, "rows": _records(s.audit)}


@app.get("/audit/kd-gate", tags=["audit"])
def audit_kd_gate(limit: int = Query(50, ge=1, le=2000)):
    s = get_store()
    from api.data_loader import _records
    if s.baseline is None:
        return {"available": False, "rows": []}
    df = s.baseline.sort_values("raw_rank").head(limit)
    return {"available": True, "rows": _records(df)}


@app.get("/classes/{regulator_class}", tags=["regulators"])
def regulators_by_class(regulator_class: str):
    if regulator_class not in ("global", "condition-specific"):
        raise HTTPException(400, "regulator_class must be 'global' or 'condition-specific'")
    return get_store().list_regulators(regulator_class=regulator_class, limit=2000)


@app.get("/programs/pca", tags=["programs"])
def programs_pca():
    """PCA scores of the transcriptional fingerprint panel (side analysis)."""
    return get_store().programs_pca()


@app.get("/programs/neighbors/{gene}", tags=["programs"])
def programs_neighbors(gene: str):
    """Transcriptomically similar regulators (nearest neighbors by fingerprint)."""
    return get_store().programs_neighbors(gene.upper())


@app.get("/programs/clusters", tags=["programs"])
def programs_clusters():
    """Agnostic perturbation clusters by fingerprint similarity (structural)."""
    return get_store().programs_clusters()


@app.get("/programs/summary", tags=["programs"])
def programs_summary():
    """Program-level headline: labeled programs (members + markers), permutation validation
    of the known complexes (SAGA/Mediator/TCR), and the promoted/demoted coherence result."""
    return get_store().programs_summary()


@app.get("/programs/findings", tags=["programs"])
def programs_findings():
    """Consolidated per-regulator findings: program label, nearest complex, neighbors, markers."""
    return get_store().programs_findings()


@app.get("/programs/classes", tags=["programs"])
def programs_classes():
    """Balanced 30-regulator panel: do different regulator classes converge on different
    downstream programs? Returns per-class target counts, interferon specificity, and the
    pairwise Jaccard of class target sets."""
    return get_store().programs_classes()


@app.get("/programs/class-targets", tags=["programs"])
def programs_class_targets(regulator_class: str = Query(..., alias="class",
                           description="regulator class name (e.g. 'SAGA/chromatin')")):
    """Convergent target list (ISG-flagged) for one regulator class. Uses a query parameter so
    class names containing '/' (e.g. SAGA/chromatin) work."""
    return get_store().class_targets(regulator_class)


@app.get("/edges/summary", tags=["edges"])
def edges_summary():
    s = get_store()
    from api.data_loader import _records
    if s.edge_summary is None:
        return {"available": False, "rows": []}
    return {"available": True, "rows": _records(s.edge_summary)}


@app.get("/edges/downstream", tags=["edges"])
def edges_downstream(limit: int = Query(50, ge=1, le=2000)):
    s = get_store()
    from api.data_loader import _records
    if s.downstream is None:
        return {"available": False, "rows": []}
    return {"available": True, "rows": _records(s.downstream.head(limit))}
