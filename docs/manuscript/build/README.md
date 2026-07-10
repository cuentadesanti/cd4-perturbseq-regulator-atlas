# Manuscript build — preprint PDF

Assembles the eight section files in `docs/manuscript/` (+ `docs/ABSTRACT.md`)
into a single submission-style **preprint PDF** via pandoc → LaTeX (tectonic).

```
bash docs/manuscript/build/build.sh
# -> docs/manuscript/build/manuscript.pdf
```

## What the pipeline does

1. **`preprocess.py`** reads the canonical `.md` sections (never modifies them) and emits `_body.md`:
   - inline `[Author Year \`DOI\`]` citations → pandoc `[@citekey]` markers (keys read from `refs.bib`);
   - the narrative `Chevalley et al. (2025, \`DOI\`)` citation likewise;
   - `Figure N (\`docs/figures/…png\`)` refs collapse to plain `Figure N`, and each figure
     is inserted as a float after the paragraph that first mentions it (author numbering 1, 2, 3, S1 kept);
   - **Table 1 is read straight from `docs/tables/operator_completion_stratified_3106.csv`**
     (`rank==7` rows) — never re-typed — so the stratum margins can't drift from source.
2. **pandoc `--citeproc`** resolves `refs.bib` against `nature.csl` → numbered superscript citations
   + a numbered bibliography.
3. **tectonic** renders the LaTeX to PDF.

## Files

| file | role |
|------|------|
| `build.sh` | entry point (re-runnable) |
| `preprocess.py` | concatenation + citation/figure/table transforms |
| `meta.yaml` | title/author/date + LaTeX preamble (unicode maps, caption styling) |
| `refs.bib` | 8 references, fetched per-DOI from CrossRef (only works cited in the body) |
| `nature.csl` | Nature numeric citation style |
| `_body.md`, `manuscript.pdf` | generated artifacts |

## Requirements

`pandoc`, `tectonic` (`brew install pandoc tectonic`). No system TeX install needed —
tectonic fetches packages on demand.

## Regenerating `refs.bib`

Each DOI is resolved via CrossRef content negotiation:

```
curl -sL -H "Accept: application/x-bibtex" "https://doi.org/<DOI>"
```

The body cites 8 DOIs (a subset of the ~15 verified in `../../literature_positioning.md`);
only cited works belong in the reference list.

## Switching venue

This is the bioRxiv/preprint build. To target a specific journal:
- swap `nature.csl` for the journal's CSL (from the citation-style-language/styles repo);
- trim the abstract (currently ~335 words) to the venue limit in `docs/ABSTRACT.md`;
- adjust `meta.yaml` (geometry, fontsize) as the template requires.
