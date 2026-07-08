"""Pydantic response models — CD4 Perturb-seq Regulator Atlas API."""
from __future__ import annotations
from typing import Optional
from pydantic import BaseModel


class Health(BaseModel):
    status: str
    n_regulators: int
    tables_loaded: dict[str, bool]


class Summary(BaseModel):
    title: str = "CD4 Perturb-seq Regulator Atlas"
    n_regulators: int
    n_global: Optional[int] = None
    n_context_specific: Optional[int] = None
    top_regulators: list[str]
    conditions: list[str]
    reproducibility_available: bool
    edges_available: bool
    programs_available: bool = False
    class_programs_available: bool = False


class RegulatorSummary(BaseModel):
    gene: str
    core_rank: int
    regulator_class: Optional[str] = None
    regpower_eb_mean: Optional[float] = None
    stability_frequency: Optional[float] = None
    guide_corr: Optional[float] = None
    donor_corr: Optional[float] = None


class EdgeItem(BaseModel):
    measured_gene: str
    condition: str
    log_fc: float
    theta_post_mean: float
    p_abs_effect_gt_1p5x: float
    direction: str


class RegulatorProfile(BaseModel):
    gene: str
    core_rank: int
    regpower_eb_mean: Optional[float] = None
    regpower_eb_sd: Optional[float] = None
    p_top_1pct: Optional[float] = None
    expected_downstream: Optional[float] = None
    # class
    regulator_class: Optional[str] = None
    condition_specificity: Optional[float] = None
    kd_significant_conditions: Optional[int] = None
    peak_condition: Optional[str] = None
    # per-condition profile (n_downstream)
    condition_profile: Optional[dict[str, int]] = None
    # quality audit
    stability_frequency: Optional[float] = None
    median_rank: Optional[float] = None
    rank_iqr: Optional[float] = None
    dropped_by_kd_gate: Optional[bool] = None
    raw_rank: Optional[int] = None
    # reproducibility audit (optional)
    guide_corr: Optional[float] = None
    donor_corr: Optional[float] = None
    single_guide_frac: Optional[float] = None
    n_guides: Optional[float] = None
    repro_weight: Optional[float] = None
    reweighted_rank: Optional[int] = None
    reweighted_score: Optional[float] = None
    donor_metadata: Optional[str] = None
    audit_status: Optional[str] = None
    audit_reason: Optional[str] = None
    # edges (optional)
    n_downstream_edges: Optional[int] = None
    top_downstream_edges: Optional[list[EdgeItem]] = None
    # transcriptional program (only for panel regulators)
    in_program_panel: Optional[bool] = None
    program_label: Optional[str] = None
    nearest_complex: Optional[str] = None
    nearest_complex_cosine: Optional[float] = None
    transcriptomic_neighbors: Optional[list[dict]] = None
    program_markers: Optional[list[dict]] = None
    # text
    interpretation: str
