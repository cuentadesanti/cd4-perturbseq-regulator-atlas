#!/usr/bin/env python3
"""Assemble the CD4+ T-cell operator manuscript into a single pandoc-ready body.

Reads the canonical section files (kept human-readable with inline
`[Author Year `DOI`]` citations and `Figure N (`path`)` refs) and emits
`_body.md` where:
  * inline citations become pandoc `[@citekey]` markers (resolved via refs.bib
    + Nature CSL by --citeproc),
  * figure/table path-refs collapse to plain `Figure N` / `Table 1` text, and
  * each figure float and the curated Table 1 are inserted after the paragraph
    that first mentions them, with author-numbered captions.

The source .md files are never modified; this only produces the build artifact.
"""
import csv
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]            # worktree root
MAN = REPO / "docs" / "manuscript"

# --- section order: abstract first, then §1..§7 -----------------------------
SECTION_FILES = [
    REPO / "docs" / "ABSTRACT.md",
    MAN / "01_introduction.md",
    MAN / "02_ranking.md",
    MAN / "03_operator.md",
    MAN / "04_programs.md",
    MAN / "05_k562.md",
    MAN / "06_discussion.md",
    MAN / "07_methods.md",
]

# --- DOI -> citekey, read straight from refs.bib so keys always match -------
def load_doi_keymap(bib: Path) -> dict:
    text = bib.read_text(encoding="utf-8")
    m = {}
    # split at each entry start; robust to one-line CrossRef entries
    for entry in re.split(r"(?=@\w+\{)", text):
        km = re.search(r"@\w+\{([^,]+),", entry)
        dm = re.search(r"DOI=\{([^}]+)\}", entry, flags=re.I)
        if km and dm:
            m[dm.group(1).strip().lower()] = km.group(1).strip()
    return m

DOI_KEY = load_doi_keymap(HERE / "refs.bib")

# --- figure floats (author numbering preserved) -----------------------------
FIGURES = {
    "1": ("docs/figures/19_reproducibility_aware_ranking_shift.png",
          "**Figure 1. Donor reproducibility is carried as an annotation, not used to "
          "re-sort the ranking.** Rank change induced when cross-donor reproducibility is "
          "folded into the breadth-based regulator ranking. Regulators that rank highly but "
          "fail cross-donor concordance — for example one carried by a single guide with no "
          "cross-guide check — stay visible at their rank, flagged rather than demoted."),
    "2": ("docs/figures/35_operator_completion_curve_3106.png",
          "**Figure 2. The operator's predictive subspace is low-dimensional and beats "
          "persistence at every rank.** Out-of-panel completion R² for the low-rank fit "
          "versus the persistence baseline across ranks, on 621 held-out out-of-panel "
          "regulators of the 3,106-regulator operator. The low-rank model beats persistence at "
          "every rank; the aggregate margin peaks at rank 7 (Δ R² = +0.379) and "
          "declines thereafter — a genuine turnover that scopes the predictive subspace to "
          "~7 dimensions."),
    "3": ("docs/figures/34_operator_cp_condition_factors_3106.png",
          "**Figure 3. The operator decomposes into gene programs gated by activation state.** "
          "CANDECOMP/PARAFAC condition-factor weights across Rest / Stim8hr / Stim48hr. Two "
          "factors are gated with bootstrap confidence intervals excluding a flat profile — "
          "one peaking in resting cells (factor 1) and one in early stimulation (factor 6) — "
          "while two further factors are constitutive, evidencing bidirectional gating rather "
          "than a single activation program."),
    "S1": ("docs/figures/32_operator_svd_scree_3106.png",
           "**Figure S1. The operator's own effective rank is far higher than its predictive "
           "subspace.** Singular-value scree of the 3,106-regulator operator. The spectrum does "
           "not collapse after seven components (effective rank ≈ 86); the ~7-dimensional "
           "structure identified in Figure 2 is the transferable part that generalizes to unseen "
           "regulators, not the operator's full rank."),
}

# --- Table 1: rank-7 (peak) slice read straight from the stratified CSV ------
TABLE1_CSV = REPO / "docs" / "tables" / "operator_completion_stratified_3106.csv"
# display order + human labels; keys are the CSV `stratum` values
TABLE1_STRATA = [
    ("all_novel", "All out-of-panel"),
    ("strong", "Strong (median split)"),
    ("weak", "Weak (median split)"),
    ("strong_panel_matched", "Strong, panel-matched"),
]


def num(x, signed=False):
    """3-decimal format; U+2212 minus; optional leading + ; never '-0.000'."""
    v = round(float(x), 3)
    if v == 0:
        v = 0.0
    s = f"{abs(v):.3f}"
    if v < 0:
        return "−" + s              # real minus sign
    return ("+" + s) if signed else s


def build_table1() -> str:
    rows = {r["stratum"]: r for r in csv.DictReader(TABLE1_CSV.open())
            if int(r["rank"]) == 7}
    missing = [k for k, _ in TABLE1_STRATA if k not in rows]
    if missing:
        sys.exit(f"Table 1: strata missing at rank 7 in CSV: {missing}")
    body = ["| Stratum | _n_ | Model R² | Persistence R² | Margin Δ R² |",
            "|:--|--:|--:|--:|--:|"]
    for key, label in TABLE1_STRATA:
        r = rows[key]
        margin = num(r["margin"], signed=True)
        if key == "all_novel":
            margin = f"**{margin}**"
        body.append(f"| {label} | {r['n_regulators']} | {num(r['r2_model'])} | "
                    f"{num(r['r2_persistence'])} | {margin} |")
    caption = (": **Table 1. The out-of-panel predictive advantage survives strength "
               "stratification.** Held-out late-stimulation completion at rank 7 (the peak "
               "rank), by knockdown-strength stratum, read from "
               "`operator_completion_stratified_3106.csv`. The margin is largest for weak "
               "regulators, where the persistence baseline collapses, but stays clearly "
               "positive even for the strongest, panel-matched regulators — for which "
               "persistence is a fair baseline — so the advantage is not an artifact of "
               "weak-regulator baseline collapse.")
    return "\n".join(body) + "\n\n" + caption


def convert_citations(text: str) -> str:
    """`[Author Year `DOI`; ...]` -> `[@k1; @k2]`; narrative Chevalley too."""
    # narrative: "Chevalley et al. (2025, `DOI`)" -> "Chevalley et al. [@key]"
    def narrative(m):
        doi = m.group(1).strip().lower()
        key = DOI_KEY.get(doi)
        return f"Chevalley et al. [@{key}]" if key else m.group(0)
    text = re.sub(r"Chevalley et al\.\s*\(20\d\d,\s*`(10\.[^`]+)`\)",
                  narrative, text, flags=re.S)

    # bracketed spans containing at least one backtick-DOI
    def bracket(m):
        inner = m.group(1)
        dois = re.findall(r"`(10\.[^`]+)`", inner)
        keys = [DOI_KEY.get(d.strip().lower()) for d in dois]
        if not keys or any(k is None for k in keys):
            return m.group(0)          # leave untouched if any unmapped
        return "[" + "; ".join(f"@{k}" for k in keys) + "]"
    text = re.sub(r"\[([^\[\]]*`10\.[^\[\]]*)\]", bracket, text)
    return text


def strip_float_refs(text: str) -> str:
    """`Figure N (`path`)` -> `Figure N`; same for Table 1."""
    text = re.sub(r"Figure(&nbsp;| )(S?\d+)\s*\(`docs/figures/[^`]+`\)",
                  r"Figure\g<1>\g<2>", text)
    text = re.sub(r"Table(&nbsp;| )(\d+)\s*\(`docs/tables/[^`]+`\)",
                  r"Table\g<1>\g<2>", text)
    return text


def insert_floats(text: str) -> str:
    """After the paragraph first mentioning each float, append the float block."""
    paras = re.split(r"\n\n+", text)
    fig_done, tbl_done = set(), False
    out = []
    for p in paras:
        out.append(p)
        for num in ("1", "2", "3", "S1"):
            if num in fig_done:
                continue
            if re.search(rf"Figure(&nbsp;| ){re.escape(num)}\b", p):
                path, cap = FIGURES[num]
                out.append(f"![{cap}]({path}){{width=92%}}")
                fig_done.add(num)
        if not tbl_done and re.search(r"Table(&nbsp;| )1\b", p):
            out.append(build_table1())
            tbl_done = True
    return "\n\n".join(out)


def main():
    parts = []
    for f in SECTION_FILES:
        raw = f.read_text(encoding="utf-8")
        if f.name == "ABSTRACT.md":
            raw = raw.replace("# Abstract (draft)", "# Abstract", 1)
        parts.append(raw.strip())
    body = "\n\n".join(parts) + "\n"

    body = convert_citations(body)
    body = strip_float_refs(body)
    body = insert_floats(body)

    # References heading so the auto-generated bibliography has a title
    body += "\n\n# References\n\n::: {#refs}\n:::\n"

    out = HERE / "_body.md"
    out.write_text(body, encoding="utf-8")

    # report leftovers a reviewer would want to know about
    leftover_doi = re.findall(r"`10\.\d{4,}/[^`]+`", body)
    leftover_path = re.findall(r"`docs/(?:figures|tables)/[^`]+`", body)
    print(f"wrote {out} ({len(body.split())} words)")
    print(f"citekeys in map: {sorted(set(DOI_KEY.values()))}")
    if leftover_doi:
        print(f"WARNING unresolved inline DOIs: {leftover_doi}", file=sys.stderr)
    if leftover_path:
        print(f"WARNING leftover asset paths: {leftover_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
