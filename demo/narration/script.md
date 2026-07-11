# Narration script

Narration is now **beat-level**: each sentence lives in a `with self.beat("…")`
block in `demo.py`, wrapping only the animation it describes, so the voiceover
syncs to that specific moment (the "336" is spoken as the counter appears, "watch
this one" as the camera zooms). Each beat holds its last frame until its audio
ends, so nothing plays over black. Edit the text in `demo.py`; this file mirrors
it. Every number is from `docs/OPERATOR_ANALYSIS.md`.

**1 · Intro — one gene → fingerprint → matrix → operator**
- Meet the regulatory operator — a map of how one human immune cell is wired.
- Take a single gene, and switch it off with CRISPR interference.
- Then measure how every other gene in the cell responds. That readout — one number per gene — is the perturbation's fingerprint.
- Repeat for the next gene, and the next — thousands of perturbations — and the fingerprints stack into one regulator-by-gene map.
- Compare every pair of fingerprints, and the map folds into one square object — the operator. Everything that follows is asking this matrix questions.

**2 · Workflow**
- Before the findings, a word on method. Every claim in this project ran the same loop.
- A question becomes a hypothesis; the hypothesis becomes an executable test; and whatever it produces gets audited — adversarially.
- Only then does a claim get to survive, weaken, or fail.

**3 · Spectral denoising (336 → 92)**
- Inside the operator, which directions are real signal — and which are just noise?
- Random-matrix theory draws a clean line. The closed-form Marchenko-Pastur edge counts three hundred and thirty-six directions above the noise.
- But that edge is optimistic. A permutation null — shuffling the data to see what noise alone produces — pushes the threshold outward, and the honest count drops to ninety-two.
- Real, reproducible structure. And of those ninety-two, only about seven generalise to unseen conditions.

**4 · Complex I reveal**
- Cleaned and clustered, the operator falls into communities — eight of them; three are stable enough to trust.
- Watch this one.
- Its members name themselves — the NDUF genes, subunits of Mitochondrial Respiratory Chain Complex One.
- An external database, CORUM, confirms the identity at a false-discovery rate near one in ten million. The clustering was blind — the model was never told what these genes were.

**5 · Skepticism**
- Good science also reports what didn't hold. Claude was asked to argue against its own results.
- Take the SAGA chromatin module. Its subunits do cluster together — but the group's stability, 0.56, sits below our 0.80 bar. So we flag it, rather than claim it.
- And the hardest test: can a regulator's own features predict the response of a regulator we have never perturbed?
- The answer is no. Real features and randomly shuffled features both score essentially zero. A clean null — reported, not buried.

**6 · Reproducibility**
- None of this asks you to take our word for it.
- Every claim ships as an artifact — a manuscript, versioned tables, matched nulls, seventeen automated tests, an open repository — the whole core rebuilding in about eight seconds on a laptop.

**7 · Closing**
- Which leaves one honest lesson. Recoverable structure is not the same as inductive predictability.
- Claude helped us find the structure — and, just as importantly, where the claim stops.

---

## Voice

- Preview/default: macOS `say` (voice **Samantha**) — offline, zero cost.
- Final target: ElevenLabs **Elise** (`voice_id EST9Ui6982FZPSi7gCHi`), wired as
  the default `eleven` voice. **Blocked on a free plan**: the API returns
  HTTP 402 ("free users cannot use library voices") for every library voice,
  Elise included. After upgrading the ElevenLabs plan to any paid tier:
  `./render.sh final eleven` renders the final with Elise, no other change.
