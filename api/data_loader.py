"""Read-only loading of the `make all` outputs into memory.

Golden rule: does NOT run models, does NOT read h5ad, does NOT touch S3. It only serves
versioned CSVs generated offline by the pipeline. If an optional table is missing, its fields
are returned as null and the API keeps working.
"""
from __future__ import annotations
import math
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
TABLES = ROOT / "docs" / "tables"
DATA = ROOT / "data" / "suppl_tables"


def _clean(v):
    if isinstance(v, float) and math.isnan(v):
        return None
    if pd.isna(v) if not isinstance(v, (list, dict)) else False:
        return None
    return v


def _records(df):
    return [{k: _clean(v) for k, v in row.items()} for row in df.to_dict("records")]


def _read(path):
    try:
        return pd.read_csv(path) if path.exists() else None
    except Exception:
        return None


class DataStore:
    """In-memory indexes built at startup."""

    def __init__(self):
        self.loaded: dict[str, bool] = {}
        self._load()

    # ---- load ----
    def _t(self, name, base=TABLES):
        df = _read(base / name)
        self.loaded[name] = df is not None
        return df

    def _read_json(self, name, base=TABLES):
        import json
        p = base / name
        try:
            return json.loads(p.read_text()) if p.exists() else None
        except Exception:
            return None

    def _load(self):
        self.ranking = self._t("hub_ranking_bayes.csv")          # master (all genes)
        self.review = self._t("top_regulators_for_review.csv")
        self.stability = self._t("hub_ranking_stability.csv")
        self.baseline = self._t("ranking_baseline_comparison.csv")
        self.classes = self._t("regulator_classes.csv")
        self.repro = self._t("hub_ranking_bayes_reproducibility_aware.csv")
        self.audit = self._t("reproducibility_audit.csv")
        self.coverage = self._t("reproducibility_coverage.csv")
        self.edges = self._t("robust_edges.csv")
        self.edge_summary = self._t("edge_summary_by_regulator.csv")
        self.downstream = self._t("top_downstream_genes.csv")
        # transcriptional programs / fingerprints (optional, `make fingerprints`)
        self.fp_pca = self._t("fingerprint_pca_scores.csv")
        self.fp_neighbors = self._t("fingerprint_neighbors.csv")
        self.fp_clusters = self._t("fingerprint_clusters.csv")
        self.fp_findings = self._t("fingerprint_findings.csv")
        self.fp_program_markers = self._t("fingerprint_program_markers.csv")
        self.fp_program_evidence = self._t("program_label_evidence.csv")
        self.fp_complex_validation = self._t("fingerprint_complex_validation.csv")
        self.fp_coherence = self._t("fingerprint_audit_coherence.csv")
        self.fp_summary = self._read_json("fingerprint_summary.json")
        self.de = self._t("DE_stats.suppl_table.csv", base=DATA)  # for the per-condition profile

        if self.ranking is None:
            raise RuntimeError(
                "Missing docs/tables/hub_ranking_bayes.csv — run `make all` before the API.")

        # per-gene indexes
        self.ranking = self.ranking.rename(columns={"target_contrast_gene_name": "gene"})
        self._idx = {name: self._index(df, key)
                     for name, df, key in [
                         ("review", self.review, "gene"),
                         ("stability", self.stability, "gene"),
                         ("baseline", self.baseline, "gene"),
                         ("classes", self.classes, "gene"),
                         ("repro", self.repro, "gene"),
                         ("audit", self.audit, "gene"),
                         ("findings", self.fp_findings, "gene"),
                     ]}
        # per-condition profile from DE_stats (n_downstream per gene×condition)
        self.cond_profile = {}
        if self.de is not None:
            g = self.de.groupby("target_contrast_gene_name")
            for gene, sub in g:
                self.cond_profile[gene] = {
                    r.culture_condition: int(r.n_downstream) for r in sub.itertuples()}
        # edges per regulator
        self.edges_by_reg = {}
        if self.edges is not None:
            for gene, sub in self.edges.groupby("perturbed_gene"):
                self.edges_by_reg[gene] = sub
        # nearest neighbors per gene (programs)
        self.fp_neighbors_by_gene = {}
        if self.fp_neighbors is not None:
            for gene, sub in self.fp_neighbors.groupby("gene"):
                self.fp_neighbors_by_gene[gene] = _records(sub.sort_values("rank"))
        # convergent downstream response genes per program
        self.program_markers = {}
        if self.fp_program_markers is not None:
            for prog, sub in self.fp_program_markers.groupby("program"):
                self.program_markers[prog] = _records(sub)

        self.conditions = (sorted(self.de["culture_condition"].unique().tolist())
                           if self.de is not None else ["Rest", "Stim8hr", "Stim48hr"])

    @staticmethod
    def _index(df, key):
        if df is None:
            return {}
        return {row[key]: {k: _clean(v) for k, v in row.items()}
                for row in df.to_dict("records")}

    # ---- queries ----
    def summary(self):
        n_glob = n_cs = None
        if self.classes is not None:
            vc = self.classes["regulator_class"].value_counts()
            n_glob = int(vc.get("global", 0))
            n_cs = int(vc.get("condition-specific", 0))
        top = self.ranking.sort_values("rank").head(10)["gene"].tolist()
        return {
            "n_regulators": int(len(self.ranking)),
            "n_global": n_glob,
            "n_context_specific": n_cs,
            "top_regulators": top,
            "conditions": self.conditions,
            "reproducibility_available": self.repro is not None,
            "edges_available": self.edges is not None,
            "programs_available": self.fp_pca is not None,
        }

    # ---- transcriptional programs / fingerprints ----
    def programs_pca(self):
        if self.fp_pca is None:
            return {"available": False, "points": []}
        cols = [c for c in ["gene", "condition", "source", "regulator_class", "cluster",
                            "program_label", "regpower_eb_mean", "n_downstream", "PC1", "PC2", "PC3"]
                if c in self.fp_pca.columns]
        return {"available": True, "points": _records(self.fp_pca[cols])}

    def programs_neighbors(self, gene):
        nb = self.fp_neighbors_by_gene.get(gene)
        if nb is None:
            return {"gene": gene, "available": gene in self.fp_neighbors_by_gene,
                    "in_panel": False, "neighbors": []}
        return {"gene": gene, "available": True, "in_panel": True, "neighbors": nb}

    def programs_clusters(self):
        if self.fp_clusters is None:
            return {"available": False, "clusters": []}
        return {"available": True, "clusters": _records(self.fp_clusters)}

    def programs_findings(self):
        if self.fp_findings is None:
            return {"available": False, "rows": []}
        return {"available": True, "rows": _records(self.fp_findings)}

    def programs_summary(self):
        """Program-level headline: labeled programs (members + markers), permutation
        validation of the known complexes, and the promoted/demoted coherence result."""
        if self.fp_pca is None:
            return {"available": False}
        evidence = _records(self.fp_program_evidence) if self.fp_program_evidence is not None else []
        for e in evidence:
            e["markers"] = self.program_markers.get(e.get("program_label"), [])
        return {
            "available": True,
            "n_panel": int(len(self.fp_pca)),
            "n_in_programs": int((self.fp_summary or {}).get("n_regulators_in_programs", 0)),
            "programs": (self.fp_summary or {}).get("programs", []),
            "program_evidence": evidence,
            "complex_validation": (_records(self.fp_complex_validation)
                                   if self.fp_complex_validation is not None else []),
            "pc1_abs_vs_ndownstream_spearman": (self.fp_summary or {}).get("pc1_abs_vs_ndownstream_spearman"),
            "audit_coherence": (self.fp_summary or {}).get("audit_coherence", {}),
        }

    def list_regulators(self, q=None, regulator_class=None, limit=50,
                        sort_by="core_rank"):
        df = self.ranking
        rows = []
        for r in df.to_dict("records"):
            gene = r["gene"]
            cls = self._idx["classes"].get(gene, {})
            rep = self._idx["repro"].get(gene, {})
            stab = self._idx["stability"].get(gene, {})
            rows.append({
                "gene": gene,
                "core_rank": int(r["rank"]),
                "regulator_class": cls.get("regulator_class"),
                "regpower_eb_mean": r.get("regpower_eb_mean"),
                "stability_frequency": stab.get("selection_frequency_top30"),
                "guide_corr": rep.get("guide_corr"),
                "donor_corr": rep.get("donor_corr"),
                "_reweighted_score": rep.get("reweighted_score"),
            })
        if q:
            ql = q.upper()
            rows = [x for x in rows if ql in x["gene"].upper()]
        if regulator_class:
            rows = [x for x in rows if x["regulator_class"] == regulator_class]
        key = {"core_rank": ("core_rank", False),
               "stability_frequency": ("stability_frequency", True),
               "reweighted_score": ("_reweighted_score", True)}.get(sort_by, ("core_rank", False))
        rows.sort(key=lambda x: (x[key[0]] is None, -(x[key[0]] or 0) if key[1] else (x[key[0]] or 0)))
        for x in rows:
            x.pop("_reweighted_score", None)
        return rows[:limit]

    def edges_for(self, gene, top=None):
        sub = self.edges_by_reg.get(gene)
        if sub is None:
            return []
        sub = sub.reindex(sub["theta_post_mean"].abs().sort_values(ascending=False).index)
        if top:
            sub = sub.head(top)
        out = []
        for r in sub.to_dict("records"):
            out.append({
                "measured_gene": r["measured_gene"],
                "condition": r["condition"],
                "log_fc": round(float(r["log_fc"]), 3),
                "theta_post_mean": round(float(r["theta_post_mean"]), 3),
                "p_abs_effect_gt_1p5x": round(float(r["p_abs_effect_gt_1p5x"]), 3),
                "direction": "up" if r["theta_post_mean"] > 0 else "down",
            })
        return out

    def profile(self, gene):
        row = self.ranking[self.ranking["gene"] == gene]
        if row.empty:
            return None
        r = {k: _clean(v) for k, v in row.iloc[0].to_dict().items()}
        rev = self._idx["review"].get(gene, {})
        cls = self._idx["classes"].get(gene, {})
        stab = self._idx["stability"].get(gene, {})
        base = self._idx["baseline"].get(gene, {})
        rep = self._idx["repro"].get(gene, {})
        aud = self._idx["audit"].get(gene, {})
        fnd = self._idx["findings"].get(gene, {})
        edges = self.edges_for(gene, top=10)
        prog = fnd.get("program_label")
        prog_nb = self.fp_neighbors_by_gene.get(gene, [])[:5]
        prog_markers = self.program_markers.get(prog, [])[:8] if prog and prog != "mixed" else []

        p = {
            "gene": gene,
            "core_rank": int(r["rank"]),
            "regpower_eb_mean": r.get("regpower_eb_mean"),
            "regpower_eb_sd": r.get("regpower_eb_sd"),
            "p_top_1pct": r.get("p_top_1pct"),
            "expected_downstream": r.get("expected_downstream"),
            "regulator_class": cls.get("regulator_class"),
            "condition_specificity": cls.get("condition_specificity"),
            "kd_significant_conditions": r.get("n_signif_conditions"),
            "peak_condition": cls.get("peak_condition") or rep.get("peak_condition"),
            "condition_profile": self.cond_profile.get(gene),
            "stability_frequency": stab.get("selection_frequency_top30"),
            "median_rank": stab.get("median_rank"),
            "rank_iqr": stab.get("rank_iqr"),
            "dropped_by_kd_gate": base.get("dropped_by_kd_gate"),
            "raw_rank": base.get("raw_rank"),
            "guide_corr": rep.get("guide_corr"),
            "donor_corr": rep.get("donor_corr"),
            "single_guide_frac": rep.get("single_guide_frac"),
            "n_guides": rep.get("n_guides_med"),
            "repro_weight": rep.get("repro_weight"),
            "reweighted_rank": rep.get("new_rank"),
            "reweighted_score": rep.get("reweighted_score"),
            "donor_metadata": rep.get("donor_metadata"),
            "audit_status": aud.get("status"),
            "audit_reason": aud.get("reason"),
            "n_downstream_edges": len(edges) or (0 if self.edges is not None else None),
            "top_downstream_edges": edges or None,
            # transcriptional program (only for panel regulators; None if not in the panel)
            "in_program_panel": bool(fnd),
            "program_label": prog,
            "nearest_complex": fnd.get("nearest_complex"),
            "nearest_complex_cosine": fnd.get("nearest_complex_cosine"),
            "transcriptomic_neighbors": prog_nb or None,
            "program_markers": prog_markers or None,
            "interpretation": self._interpret(gene, cls, stab, base, aud, rep, fnd),
        }
        return p

    @staticmethod
    def _interpret(gene, cls, stab, base, aud, rep, fnd=None):
        klass = cls.get("regulator_class")
        parts = []
        if klass == "global":
            parts.append("Global regulator (stable effect across conditions)")
        elif klass == "condition-specific":
            pk = cls.get("peak_condition")
            parts.append(f"Context-specific regulator (peaks in {pk})" if pk
                         else "Context-specific regulator")
        else:
            parts.append("Regulator")
        checks = []
        if base.get("dropped_by_kd_gate") is False:
            checks.append("validated knockdown")
        sf = stab.get("selection_frequency_top30")
        if sf is not None:
            checks.append(f"top-30 in {sf*100:.0f}% of bootstraps")
        st = aud.get("status")
        if st == "survives":
            checks.append("survives guide/donor reproducibility")
        elif st == "demoted":
            checks.append("demoted by reproducibility audit")
        elif st == "promoted":
            checks.append("promoted by reproducibility audit")
        if checks:
            parts.append("; ".join(checks))
        prog = (fnd or {}).get("program_label")
        if prog and prog != "mixed":
            nc = (fnd or {}).get("nearest_complex")
            cos = (fnd or {}).get("nearest_complex_cosine")
            parts.append(f"Maps to the {prog} program by transcriptional fingerprint"
                         + (f" (nearest {nc} centroid, cosine {cos})" if nc and cos is not None else ""))
        return ". ".join(parts) + "."


store: DataStore | None = None


def get_store() -> DataStore:
    global store
    if store is None:
        store = DataStore()
    return store
