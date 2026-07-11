# The regulatory operator — Manim demo

A ~2:10 continuous visual essay (Manim Community v0.20) that tells the operator
story from `docs/OPERATOR_ANALYSIS.md`: 1.8 TB → one correlation operator → the
336-vs-92 spectral denoising → the Complex I community naming itself → the honest
nulls → the closing thesis. Built as one `MovingCameraScene` so objects persist
and transform rather than cutting between slides.

## Quick start

```bash
# from repo root; uses the repo venv (.venv) which already has manim installed
cd demo
./render.sh preview          # fast 480p, silent — iterate on this
./render.sh final            # 1080p30, silent
./render.sh final say        # 1080p30 + offline macOS narration (no API key)
```

Output lands in `media/videos/demo/`.

## Structure

| file | role |
|------|------|
| `demo.py` | the `FullDemo` scene: 7 chapter methods + the `NAR` narration |
| `theme.py` | the design system — semantic colors, type scale, easing grammar |
| `components.py` | 5 reusable parts: `EvidenceCard`, `MetricGate`, `EvidenceChain`, `ClaimVerdict`, `spotlight` |
| `voice.py` | `SayService` (offline macOS TTS) and `ElevenLabsByID` (paid API) |
| `assets/` | real figures from `docs/figures/`, cropped for the reveal |
| `narration/script.md` | the spoken script (mirrors `NAR` in `demo.py`) |

Color is semantic, not decorative: blue = supported structure, orange = the
process / Claude, green = survived, yellow = limited evidence, red = rejected /
null, gray = noise & baselines.

## Narration

Each chapter is wrapped in a single `self.voiceover(...)`, so the visuals
auto-stretch to the spoken length — no hand-tuned `wait()` calibration. No
bookmarks are used, so no transcription model is needed.

- **`say`** (default, works now): macOS `say` → mp3, zero cost, offline.
  Pick a voice with `SAY_VOICE=Samantha ./render.sh final say`.
- **`eleven`** (higher quality, needs a paid key): put the key in
  `eleven_labs_raw_key.txt` and an allowed voice id in `ELEVEN_VOICE_ID`, then
  `./render.sh final eleven`. The provided free key returns **HTTP 402**
  ("free users cannot use library voices via the API"), so this path needs a
  plan upgrade or a personal voice id.

## Notes / constraints

- **No LaTeX on this machine**, so all text is Pango `Text` (not `Tex`/`MathTex`)
  and numeric counters are driven by a `ValueTracker`, not `Integer`. Keep it
  that way unless you install a TeX distribution.
- Font is `Helvetica Neue` (macOS built-in); change `FONT` in `theme.py` for Inter.
- The heatmap is raster (`assets/heatmap_panel.png`); everything overlaid on it —
  the highlight box, NDUF labels, identity card — is vector so it stays crisp
  under the camera zoom.
