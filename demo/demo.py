"""The regulatory operator — a ~3 minute continuous visual essay.

One MovingCameraScene, seven chapters. The camera is used as an argument
(zoom into a discovery, pull back for the epistemic turn), objects persist
and change meaning, and every number is animated from the real analysis in
docs/OPERATOR_ANALYSIS.md.

Render:
    manim -pql demo.py FullDemo      # fast preview
    manim -qh  demo.py FullDemo      # 1080p delivery
    manim -qh demo.py FullDemo -n 0,1  # single chapter by section index
"""

import os
from contextlib import contextmanager
from pathlib import Path
import numpy as np

from manim_voiceover import VoiceoverScene

from manim import (
    MovingCameraScene,
    VGroup,
    Group,
    ImageMobject,
    Text,
    Axes,
    Rectangle,
    RoundedRectangle,
    Line,
    DashedLine,
    Arrow,
    Dot,
    ValueTracker,
    always_redraw,
    FadeIn,
    FadeOut,
    Write,
    Create,
    GrowFromCenter,
    Transform,
    ReplacementTransform,
    FadeTransform,
    LaggedStart,
    Indicate,
    UP,
    DOWN,
    LEFT,
    RIGHT,
    ORIGIN,
    config,
)

import theme as T
from components import (
    txt,
    title,
    EvidenceCard,
    MetricGate,
    EvidenceChain,
    ClaimVerdict,
    spotlight,
)

ASSETS = Path(__file__).parent / "assets"
config.background_color = T.BG



class FullDemo(VoiceoverScene, MovingCameraScene):
    # ------------------------------------------------------------------ #
    def construct(self):
        self.base_width = self.camera.frame.get_width()
        self.base_center = self.camera.frame.get_center().copy()

        # NARRATE=say (default off) enables offline macOS voiceover;
        # NARRATE=eleven uses ElevenLabs (needs a paid key + voice id).
        self.narrate = os.environ.get("NARRATE", "").lower()
        if self.narrate == "say":
            from voice import SayService
            self.set_speech_service(
                SayService(voice=os.environ.get("SAY_VOICE", "Samantha"),
                           ffmpeg=os.environ.get("FFMPEG", "ffmpeg")))
        elif self.narrate == "eleven":
            from voice import ElevenLabsByID
            # Elise. Works once the ElevenLabs plan is upgraded off Free
            # (free plan 402s on all library voices via the API).
            self.set_speech_service(
                ElevenLabsByID(
                    voice_id=os.environ.get("ELEVEN_VOICE_ID",
                                            "EST9Ui6982FZPSi7gCHi"),
                    model="eleven_turbo_v2_5"))

        for name, fn in [
            ("intro", self.intro),
            ("workflow", self.claude_workflow),
            ("denoising", self.spectral_denoising),
            ("complex_i", self.complex_i_reveal),
            ("skepticism", self.scientific_skepticism),
            ("reproducibility", self.reproducibility),
            ("closing", self.closing),
        ]:
            self.next_section(name)
            self.run_chapter(name, fn)

    def run_chapter(self, name, fn):
        # No blanket clear — chapters cross-fade (see crossfade_in) so a
        # protagonist (the operator) can persist and there is never a
        # fade-to-empty between chapters.
        self.camera.frame.move_to(self.base_center).set_width(self.base_width)
        fn()

    def crossfade_in(self, *news, run_time=0.7, keep=()):
        """Fade the current frame out while the next content fades in — no
        empty gap between chapters. Pass keep=[mob] so a protagonist stays."""
        keep = list(keep)
        outs = [FadeOut(m) for m in self.mobjects if m not in keep]
        ins = [FadeIn(n) for n in news]
        if outs or ins:
            self.play(*outs, *ins, run_time=run_time)

    @contextmanager
    def beat(self, text):
        """Narrate one segment while its animations play, syncing at the
        sentence level. The voiceover holds the last frame until the audio
        finishes, so a beat is never shorter than its narration and never
        plays over black. Silent when narration is off."""
        if self.narrate in ("say", "eleven"):
            with self.voiceover(text=text) as tr:
                yield tr
        else:
            yield None
            self.wait(0.4)

    # ---- camera helpers ---------------------------------------------- #
    def focus_on(self, mob, width=None, run_time=1.4):
        self.play(
            self.camera.frame.animate.move_to(mob).set(
                width=width or mob.width * 1.4
            ),
            run_time=run_time,
            rate_func=T.EASE_REVEAL,
        )

    def reset_camera(self, run_time=1.2):
        self.play(
            self.camera.frame.animate.move_to(self.base_center).set_width(
                self.base_width
            ),
            run_time=run_time,
            rate_func=T.EASE_MOVE,
        )

    def clear_all(self, run_time=0.6):
        if self.mobjects:
            self.play(
                *[FadeOut(m) for m in self.mobjects], run_time=run_time
            )

    # ================================================================== #
    # 1 — INTRO: 1.8 TB collapses into one operator
    # ================================================================== #
    def intro(self):
        rng = np.random.default_rng(3)
        cell = 0.26

        with self.beat("We begin with a genome-scale screen of primary human "
                       "T-cells: thousands of gene regulators, each switched "
                       "off in turn, while single-cell sequencing reads the "
                       "response."):
            head = title("The regulatory operator", color=T.SIGNAL)
            sub = txt("Genome-scale CRISPRi Perturb-seq · primary human "
                      "CD4⁺ T cells", size=T.SUB_SIZE, color=T.FG)
            VGroup(head, sub).arrange(DOWN, buff=0.4)
            scale = txt("3,106 regulators · 4 donors × 3 states · single-cell "
                        "readout", size=T.SMALL_SIZE, color=T.MUTED)
            scale.next_to(sub, DOWN, buff=0.45)
            self.play(Write(head), run_time=T.T_NORMAL)
            self.play(FadeIn(sub, shift=UP * 0.2), run_time=T.T_MICRO)
            self.play(FadeIn(scale), run_time=T.T_MICRO)
        self.play(FadeOut(sub), FadeOut(scale),
                  head.animate.scale(0.5).to_edge(UP), run_time=0.6)

        # --- where does one number come from? one regulator, one gene ---
        regA = VGroup(Dot(radius=0.16, color=T.NEGATIVE),
                      txt("Regulator A", size=T.SMALL_SIZE, color=T.FG))
        regA[1].next_to(regA[0], DOWN, buff=0.15)
        regA.move_to([-5.1, 1.5, 0])
        cross = txt("✕", size=24, color=T.NEGATIVE, weight="BOLD").move_to(regA[0])
        crispr = txt("CRISPRi", size=T.SMALL_SIZE, color=T.NEGATIVE).next_to(
            regA[0], UP, buff=0.15)
        geneB = VGroup(Dot(radius=0.13, color=T.SIGNAL),
                       txt("Gene B", size=T.SMALL_SIZE, color=T.FG))
        geneB[1].next_to(geneB[0], DOWN, buff=0.15)
        geneB.move_to([-5.1, -1.4, 0])
        arr0 = Arrow(regA[0].get_bottom(), geneB[0].get_top(), buff=0.25,
                     color=T.MUTED, stroke_width=T.STROKE_MUTED)
        with self.beat("Take one regulator. Inhibit it, and follow a single "
                       "downstream gene."):
            self.play(FadeIn(regA), FadeIn(crispr), FadeIn(cross),
                      run_time=T.T_NORMAL)
            self.play(Create(arr0), FadeIn(geneB), run_time=T.T_NORMAL)

        # --- many cells: control vs perturbed distributions ---
        axis = Line([-1.4, -1.7, 0], [4.4, -1.7, 0], color=T.MUTED, stroke_width=2)
        axlab = txt("expression of Gene B", size=T.SMALL_SIZE,
                    color=T.MUTED).next_to(axis, DOWN, buff=0.15)

        def cloud(cx, cy, color, n=16):
            return VGroup(*[
                Dot(radius=0.05, color=color).move_to(
                    [cx + rng.normal(0, 0.35), cy + rng.uniform(-0.3, 0.3), 0])
                for _ in range(n)])
        ctrl = cloud(0.6, 0.3, T.MUTED)
        pert = cloud(0.6, -0.7, T.SIGNAL)
        ctrl_lab = txt("control", size=T.SMALL_SIZE, color=T.MUTED).next_to(
            ctrl, LEFT, buff=0.4)
        pert_lab = txt("perturbed", size=T.SMALL_SIZE, color=T.SIGNAL).next_to(
            pert, LEFT, buff=0.4)
        with self.beat("Across many cells, the screen compares control to "
                       "perturbed — the whole distribution shifts."):
            self.play(Create(axis), FadeIn(axlab), run_time=T.T_MICRO)
            self.play(FadeIn(ctrl), FadeIn(ctrl_lab),
                      FadeIn(pert), FadeIn(pert_lab), run_time=T.T_NORMAL)
            self.play(pert.animate.shift(RIGHT * 1.4), run_time=T.T_NORMAL,
                      rate_func=T.EASE_MOVE)

        # --- collapse to one standardized effect ---
        zval = txt("z = +2.4", size=T.SUB_SIZE, color=T.POSITIVE,
                   weight="BOLD").move_to([1.5, -0.2, 0])
        znote = txt("standardized effect · depth-independent", size=T.SMALL_SIZE,
                    color=T.MUTED).next_to(zval, DOWN, buff=0.25)
        with self.beat("That shift collapses into one standardized effect — "
                       "direction and evidence in a single number."):
            self.play(FadeOut(axis), FadeOut(axlab), FadeOut(ctrl_lab),
                      FadeOut(pert_lab),
                      FadeTransform(VGroup(ctrl, pert), zval), run_time=T.T_REVEAL)
            self.play(FadeIn(znote), run_time=T.T_MICRO)

        # --- every gene's effect stacks into the fingerprint ---
        vals = rng.normal(0, 1, 16)
        vals[0] = 2.4  # the z we just built = the top cell
        vec = VGroup(*[
            Rectangle(width=cell, height=cell, stroke_width=0,
                      fill_color=(T.POSITIVE if i == 0 else
                                  (T.SIGNAL if v > 0 else T.NEGATIVE)),
                      fill_opacity=float(min(abs(v), 1)) * 0.9)
            for i, v in enumerate(vals)
        ]).arrange(DOWN, buff=0.02).move_to(RIGHT * 3.7 + DOWN * 0.3)
        veclab = txt("fingerprint\n2,000 genes", size=T.SMALL_SIZE,
                     color=T.MUTED).next_to(vec, RIGHT, buff=0.3)
        with self.beat("Do that for every gene — thousands of standardized "
                       "effects stack into one transcriptional fingerprint."):
            self.play(FadeOut(regA), FadeOut(cross), FadeOut(crispr),
                      FadeOut(arr0), FadeOut(geneB), FadeOut(znote), run_time=0.4)
            self.play(ReplacementTransform(zval, vec[0]), run_time=0.7)
            self.play(LaggedStart(*[GrowFromCenter(s) for s in vec[1:]],
                                  lag_ratio=0.05), FadeIn(veclab), run_time=1.4)

        # stack thousands of fingerprints into the regulator x gene map
        rows, cols = 16, 22
        colgroups = []
        for j in range(cols):
            cvals = rng.normal(0, 1, rows)
            cg = VGroup()
            for i, v in enumerate(cvals):
                sq = Rectangle(width=cell, height=cell, stroke_width=0,
                               fill_color=(T.SIGNAL if v > 0 else T.NEGATIVE),
                               fill_opacity=float(min(abs(v), 1)) * 0.9)
                sq.move_to([(j - cols / 2 + 0.5) * cell,
                            (rows / 2 - i - 0.5) * cell, 0])
                cg.add(sq)
            colgroups.append(cg)
        mat = VGroup(*colgroups).move_to(DOWN * 0.3)
        mlab = txt("3,106 regulators × 2,000 genes", size=T.LABEL_SIZE,
                   color=T.MUTED).next_to(mat, DOWN, buff=0.4)
        with self.beat("One fingerprint per regulator — stack thousands into "
                       "the regulator-by-gene map."):
            self.play(FadeOut(veclab), run_time=0.3)
            self.play(ReplacementTransform(vec, colgroups[0]), run_time=0.8)
            self.play(LaggedStart(*[FadeIn(cg) for cg in colgroups[1:]],
                                  lag_ratio=0.05), run_time=1.6)
            self.play(FadeIn(mlab), run_time=0.4)

        # fold into the square operator
        with self.beat("Correlate every pair, and the map folds into one "
                       "square object — the operator."):
            op = self.make_operator_matrix(n=20, size=5.0).move_to(DOWN * 0.3)
            self.play(FadeTransform(mat, op), FadeOut(mlab), run_time=1.4)
            cap = txt("the operator — regulator × regulator", size=T.LABEL_SIZE,
                      color=T.SIGNAL).next_to(op, DOWN, buff=0.4)
            self.play(FadeIn(cap), run_time=0.5)
        self.operator = op          # protagonist: persists into denoising
        self.op_cap = cap

    def make_operator_matrix(self, n=22, size=4.2):
        """A synthetic block-correlation grid — the abstract operator."""
        rng = np.random.default_rng(7)
        blocks = [(0, 6), (6, 12), (12, 16)]  # a few coherent communities
        M = rng.normal(0, 0.25, (n, n))
        for a, b in blocks:
            v = rng.normal(0, 1, (b - a, 1))
            M[a:b, a:b] += 0.9 * (v @ v.T) / (np.abs(v).max() ** 2)
        M = (M + M.T) / 2
        np.fill_diagonal(M, 1.0)
        M = np.clip(M / np.abs(M).max(), -1, 1)

        cell = size / n
        grid = VGroup()
        for i in range(n):
            for j in range(n):
                val = M[i, j]
                color = T.SIGNAL if val > 0 else T.NEGATIVE
                sq = Rectangle(
                    width=cell, height=cell, stroke_width=0,
                    fill_color=color, fill_opacity=float(min(abs(val), 1)) * 0.9,
                )
                sq.move_to([(j - n / 2 + 0.5) * cell,
                            (n / 2 - i - 0.5) * cell, 0])
                grid.add(sq)
        return grid

    # ================================================================== #
    # 2 — WORKFLOW: how every claim was earned
    # ================================================================== #
    def claude_workflow(self):
        head = txt("Every claim ran the same loop", size=T.SUB_SIZE,
                   color=T.CLAUDE, weight="BOLD").to_edge(UP, buff=0.7)
        steps = ["Question", "Hypothesis", "Executable test",
                 "Artifact", "Adversarial audit"]
        rows = VGroup(*[
            VGroup(
                RoundedRectangle(corner_radius=0.12, width=4.4, height=0.72,
                                 stroke_width=T.STROKE_MUTED,
                                 stroke_color=T.CLAUDE,
                                 fill_color=T.CLAUDE, fill_opacity=0.06),
                txt(s, size=T.LABEL_SIZE, color=T.FG),
            )
            for s in steps
        ])
        for r in rows:
            r[1].move_to(r[0])
        rows.arrange(DOWN, buff=0.35).next_to(head, DOWN, buff=0.5)
        arrows = VGroup(*[
            Arrow(rows[i].get_bottom(), rows[i + 1].get_top(), buff=0.06,
                  color=T.CLAUDE, stroke_width=T.STROKE_MUTED,
                  max_tip_length_to_length_ratio=0.2)
            for i in range(len(rows) - 1)
        ])
        with self.beat("A word on method. Every claim ran the same loop —"):
            # park the operator on the side; it returns in denoising
            self.crossfade_in(head, keep=[self.operator])
            self.play(self.operator.animate.scale(0.5).to_edge(LEFT, buff=0.4),
                      FadeOut(self.op_cap), run_time=0.8)
        with self.beat("question, hypothesis, executable test, artifact — each "
                       "one audited adversarially."):
            self.play(
                LaggedStart(
                    *[FadeIn(m) for pair in zip(rows, list(arrows) + [None])
                      for m in ([pair[0], pair[1]] if pair[1] else [pair[0]])],
                    lag_ratio=0.4,
                ),
                run_time=2.2,
            )

        # concrete cases — this is where Claude actually bit
        ex = VGroup(
            txt("“Is SAGA a stable community?”   →   0.56 < 0.80   →   weakened",
                size=T.SMALL_SIZE, color=T.WARNING),
            txt("“Predict an unseen regulator?”   →   real ≈ shuffled   →   failed",
                size=T.SMALL_SIZE, color=T.NEGATIVE),
        ).arrange(DOWN, buff=0.3).next_to(rows, DOWN, buff=0.5)
        with self.beat("Concretely: is SAGA a stable community? The gate "
                       "weakened it. Predict an unseen regulator? The control "
                       "failed it."):
            self.play(LaggedStart(*[FadeIn(e, shift=UP * 0.15) for e in ex],
                                  lag_ratio=0.4), run_time=T.T_REVEAL)

    # ================================================================== #
    # 3 — DENOISING: 336 signal directions become 92
    # ================================================================== #
    def spectral_denoising(self):
        head = txt("Separating signal from noise", size=T.SUB_SIZE,
                   color=T.SIGNAL, weight="BOLD").to_edge(UP, buff=0.7)

        ax = Axes(
            x_range=[0, 3.2, 0.5], y_range=[0, 1.9, 0.5],
            x_length=9.5, y_length=4.2,
            axis_config={"color": T.MUTED, "stroke_width": T.STROKE_MUTED,
                         "include_ticks": True, "include_numbers": False},
            tips=False,
        ).to_edge(DOWN, buff=1.0)
        xlab = txt("eigenvalue", size=T.SMALL_SIZE, color=T.MUTED).next_to(
            ax, DOWN, buff=0.15)

        # Marchenko-Pastur bulk + a signal tail poking above the edge
        lam_minus, lam_plus, sigma2, q = 0.039, 1.49, 0.5, 0.52

        def mp(x):
            if x <= lam_minus or x >= lam_plus:
                return 0.0
            return np.sqrt((lam_plus - x) * (x - lam_minus)) / (
                2 * np.pi * q * sigma2 * x)

        xs = np.linspace(0.03, 3.15, 70)
        bars = VGroup()
        bar_meta = []  # (bar, x)
        for x in xs:
            h = mp(x)
            if x >= lam_plus:  # signal tail: low, declining, above the edge
                h = 0.28 * np.exp(-(x - lam_plus) / 1.4) + 0.05
                is_signal = True
            else:
                is_signal = False
            top = ax.c2p(x, h)
            base = ax.c2p(x, 0)
            height = max(top[1] - base[1], 0.001)
            bar = Rectangle(
                width=(9.5 / 70) * 0.9, height=height, stroke_width=0,
                fill_color=(T.SIGNAL if is_signal else T.MUTED),
                fill_opacity=0.85,
            )
            bar.move_to([base[0], base[1] + height / 2, 0])
            bars.add(bar)
            bar_meta.append((bar, x))

        with self.beat("Now the questions. Inside the operator, which "
                       "directions are real — and which are noise?"):
            # the parked operator becomes its own eigenvalue spectrum
            self.play(*[FadeOut(m) for m in self.mobjects
                        if m is not self.operator],
                      Create(ax), FadeIn(head, shift=DOWN * 0.2),
                      run_time=T.T_NORMAL)
            self.play(FadeTransform(self.operator, bars), FadeIn(xlab),
                      run_time=1.4)

        # global mode annotation (off-scale)
        gmode = txt("global mode  λ₀ = 167  (off-scale)", size=T.SMALL_SIZE,
                    color=T.MUTED).next_to(head, DOWN, buff=0.35).to_edge(LEFT, buff=1.2)
        # theoretical edge 1.49  ->  336 signal eigenvalues
        edge_t = DashedLine(ax.c2p(lam_plus, 0), ax.c2p(lam_plus, 1.8),
                            color=T.WARNING, stroke_width=T.STROKE_MAIN)
        edge_t_lab = txt("closed-form edge  λ₊ = 1.49", size=T.SMALL_SIZE,
                         color=T.WARNING).next_to(edge_t, UP, buff=0.1)
        # LaTeX-free numeric counter driven by a ValueTracker
        count_val = ValueTracker(336)
        c_anchor = np.array([3.4, 2.35, 0])
        counter = always_redraw(
            lambda: txt(f"{int(round(count_val.get_value()))}",
                        size=T.TITLE_SIZE, weight="BOLD",
                        color=(T.SIGNAL if count_val.get_value() <= 100
                               else T.WARNING)).move_to(c_anchor))
        counter_static = txt("336", size=T.TITLE_SIZE, color=T.WARNING,
                             weight="BOLD").move_to(c_anchor)
        counter_lab = txt("admitted · closed-form edge", size=T.SMALL_SIZE,
                          color=T.WARNING).move_to(c_anchor + DOWN * 0.7)
        with self.beat("Random-matrix theory draws a line. The closed-form "
                       "Marchenko-Pastur edge admits three hundred and "
                       "thirty-six directions."):
            self.play(FadeIn(gmode), Create(edge_t), FadeIn(edge_t_lab),
                      run_time=T.T_NORMAL)
            self.play(FadeIn(counter_static, shift=UP * 0.2),
                      FadeIn(counter_lab, shift=UP * 0.2), run_time=T.T_NORMAL)

        # empirical permutation null: slide the edge to 2.95, demote the tail
        edge_e = DashedLine(ax.c2p(2.95, 0), ax.c2p(2.95, 1.8),
                            color=T.NEGATIVE, stroke_width=T.STROKE_MAIN)
        edge_e_lab = txt("empirical null  λ₊ = 2.95", size=T.SMALL_SIZE,
                         color=T.NEGATIVE).next_to(edge_e, UP, buff=0.1)
        demoted = [b for b, x in bar_meta if lam_plus <= x < 2.95]
        counter_lab2 = txt("supported · empirical null", size=T.SMALL_SIZE,
                           color=T.SIGNAL).move_to(counter_lab)
        self.remove(counter_static)
        self.add(counter)
        with self.beat("But that edge is optimistic. A permutation null — the "
                       "same data with the signal shuffled out — pushes it "
                       "outward, and the honest count drops to ninety-two."):
            self.play(
                LaggedStart(*[b.animate.set_fill(T.MUTED) for b in demoted],
                            lag_ratio=0.02),
                count_val.animate.set_value(92),
                Create(edge_e), FadeIn(edge_e_lab),
                run_time=T.T_REVEAL,
            )
            self.play(Transform(counter_lab, counter_lab2), run_time=T.T_NORMAL)

        # freeze the live counter so the chapter can cross-fade cleanly
        counter92 = txt("92", size=T.TITLE_SIZE, color=T.SIGNAL,
                        weight="BOLD").move_to(c_anchor)
        self.remove(counter)
        self.add(counter92)

        note = txt("real structure — of which only ~7 transfer across "
                   "held-out conditions", size=T.SMALL_SIZE, color=T.MUTED)
        note.next_to(ax, UP, buff=0.2).to_edge(LEFT, buff=1.2)
        with self.beat("Real, reproducible structure — of which only about "
                       "seven carry across held-out conditions."):
            self.play(FadeIn(note), run_time=T.T_MICRO)

    # ================================================================== #
    # 4 — COMPLEX I: camera enters a community; it names itself
    # ================================================================== #
    def complex_i_reveal(self):
        heat = ImageMobject(str(ASSETS / "heatmap_panel.png"))
        heat.height = 6.6
        heat.move_to(ORIGIN)
        head = txt("8 regulator communities · 3 clear the ≥ 0.80 stability gate",
                   size=T.LABEL_SIZE, color=T.MUTED).to_edge(UP, buff=0.5)
        with self.beat("Cleaned and clustered, the operator falls into "
                       "communities — eight; three stable enough to trust."):
            self.crossfade_in(heat, head, run_time=0.9)

        # target: the compact strong module in the extreme bottom-right corner
        # (bounds read off assets/heatmap_panel.png).
        w_img = heat.width
        h_img = heat.height
        bx = heat.get_center()[0] + 0.44 * w_img
        by = heat.get_center()[1] - 0.44 * h_img
        box = Rectangle(width=0.12 * w_img, height=0.12 * h_img,
                        stroke_color=T.POSITIVE, stroke_width=T.STROKE_MAIN)
        box.move_to([bx, by, 0])

        # labels go to the LEFT of the box — the corner leaves no room on the
        # right — and the camera frames box + labels together.
        genes = ["NDUFA9", "NDUFS2", "NDUFV1", "NDUFS3", "NDUFB8", "NDUFS7"]
        gene_col = VGroup(*[
            txt(g, size=13, color=T.POSITIVE) for g in genes
        ]).arrange(DOWN, buff=0.1, aligned_edge=RIGHT)
        gene_col.next_to(box, LEFT, buff=0.4)
        focus_grp = Group(box, gene_col)

        mask = spotlight(box, opacity=0.78)
        with self.beat("Watch this one."):
            self.play(FadeIn(mask), Create(box), FadeOut(head),
                      run_time=T.T_NORMAL)
            self.focus_on(focus_grp, width=focus_grp.width * 1.4, run_time=1.6)

        with self.beat("Its members name themselves — the N-D-U-F genes, "
                       "subunits of Mitochondrial Respiratory Chain Complex One."):
            self.play(LaggedStart(*[FadeIn(g, shift=LEFT * 0.1)
                                    for g in gene_col],
                                  lag_ratio=0.25), run_time=2.0)

        ident = VGroup(
            txt("Community 7 · n = 87 · 45% donor-robust",
                size=T.LABEL_SIZE, color=T.FG),
            txt("Mitochondrial Respiratory Chain Complex I",
                size=T.BODY_SIZE, color=T.POSITIVE, weight="BOLD"),
            txt("CORUM BH-FDR = 1.4 × 10⁻⁷   (8/13 NADH-dehydrogenase subunits)",
                size=T.LABEL_SIZE, color=T.FG),
        ).arrange(DOWN, buff=0.25)
        panel = RoundedRectangle(
            corner_radius=0.15, width=ident.width + 1.0, height=ident.height + 0.7,
            stroke_width=T.STROKE_MUTED, stroke_color=T.POSITIVE,
            fill_color=T.FAINT, fill_opacity=0.5,
        ).move_to(ident)
        stamp = txt("community fixed first · CORUM queried after",
                    size=T.SMALL_SIZE, color=T.MUTED)
        group = VGroup(panel, ident).move_to(UP * 0.3)
        stamp.next_to(group, DOWN, buff=0.6)
        with self.beat("CORUM, an external database, confirms it — a false-"
                       "discovery rate near one in ten million. The clustering "
                       "used no annotations; we asked what it contained only "
                       "afterward."):
            self.play(FadeOut(gene_col), run_time=0.3)
            self.reset_camera(run_time=1.2)
            self.play(FadeOut(mask), FadeOut(heat), FadeOut(box), run_time=0.5)
            self.play(FadeIn(panel), run_time=0.5)
            self.play(LaggedStart(*[FadeIn(line, shift=UP * 0.1)
                                    for line in ident], lag_ratio=0.3),
                      run_time=T.T_REVEAL)
            self.play(FadeIn(stamp), run_time=T.T_NORMAL)

    # ================================================================== #
    # 5 — SKEPTICISM: what did NOT survive
    # ================================================================== #
    def scientific_skepticism(self):
        head = txt("Claude argued against the result too", size=T.SUB_SIZE,
                   color=T.CLAUDE, weight="BOLD").to_edge(UP, buff=0.7)
        with self.beat("Good science also reports what fails — so Claude "
                       "argued against its own results."):
            self.crossfade_in(head)

        gate = MetricGate("SAGA core module stability", value=0.56, gate=0.80,
                          lo=0.0, hi=1.0, width=7.0, accent=T.SIGNAL)
        gate.move_to(UP * 0.9)
        start_pos = gate.marker_position(0.0)
        gate.marker.move_to(start_pos)
        gate.val_lab.next_to(gate.marker, DOWN, buff=0.2)
        end_pos = gate.marker_position(0.56)
        verdicts = VGroup(
            ClaimVerdict("SAGA subunits co-cluster", "observed", ok=True),
            ClaimVerdict("a stable Leiden community", "below 0.80 gate", ok=False),
            ClaimVerdict("a convergent module", "supported independently", ok=True),
        ).arrange(DOWN, buff=0.28, aligned_edge=LEFT).next_to(gate, DOWN, buff=0.8)
        with self.beat("Take the SAGA module. Its subunits co-cluster — but as "
                       "a stable community it misses the gate. It stands only "
                       "as a convergent module, supported independently."):
            self.play(FadeIn(gate), run_time=T.T_NORMAL)
            self.play(
                gate.marker.animate.move_to(end_pos).set_color(T.WARNING),
                gate.val_lab.animate.move_to(end_pos + DOWN * 0.45).set_color(T.WARNING),
                run_time=T.T_REVEAL, rate_func=T.EASE_MOVE,
            )
            self.play(LaggedStart(*[FadeIn(v, shift=UP * 0.15) for v in verdicts],
                                  lag_ratio=0.4), run_time=T.T_REVEAL)

        # inductive null: two rows land on one shared zero line (x ≈ R²)
        q = txt("Can a regulator's features predict an unseen regulator?",
                size=T.LABEL_SIZE, color=T.FG).next_to(head, DOWN, buff=0.55)
        zero = DashedLine(UP * 1.0, DOWN * 1.0, color=T.MUTED,
                          stroke_width=2).move_to(DOWN * 0.4)
        x0 = zero.get_center()[0]
        zlab = txt("held-out R² = 0", size=T.SMALL_SIZE, color=T.MUTED).next_to(
            zero, UP, buff=0.15)
        ay = zero.get_center()[1] + 0.4
        by = zero.get_center()[1] - 0.4
        real_lab = txt("real features", size=T.SMALL_SIZE,
                       color=T.SIGNAL).move_to([-4.3, ay, 0])
        shuf_lab = txt("shuffled features", size=T.SMALL_SIZE,
                       color=T.MUTED).move_to([-4.3, by, 0])
        real_dot = Dot(radius=0.11, color=T.SIGNAL).next_to(real_lab, RIGHT, buff=0.3)
        shuf_dot = Dot(radius=0.11, color=T.MUTED).next_to(shuf_lab, RIGHT, buff=0.3)
        with self.beat("The hardest test: can a regulator's own features "
                       "predict a regulator we have never perturbed?"):
            self.play(FadeOut(gate), FadeOut(verdicts), run_time=0.5)
            self.play(FadeIn(q), Create(zero), FadeIn(zlab), run_time=T.T_NORMAL)
            self.play(FadeIn(real_lab), FadeIn(shuf_lab),
                      FadeIn(real_dot), FadeIn(shuf_dot), run_time=T.T_MICRO)

        real_val = txt("−0.0005", size=T.SMALL_SIZE, color=T.SIGNAL)
        shuf_val = txt("−0.00003", size=T.SMALL_SIZE, color=T.MUTED)
        concl = txt("real ≈ shuffled ≈ 0   —   a clean, reported null",
                    size=T.LABEL_SIZE, color=T.NEGATIVE,
                    weight="BOLD").next_to(zero, DOWN, buff=1.0)
        with self.beat("No. Real features and shuffled features both land on "
                       "zero. A clean null — reported, not buried."):
            self.play(real_dot.animate.move_to([x0 - 0.08, ay, 0]),
                      shuf_dot.animate.move_to([x0 - 0.02, by, 0]),
                      run_time=T.T_REVEAL, rate_func=T.EASE_MOVE)
            real_val.next_to(real_dot, RIGHT, buff=0.25)
            shuf_val.next_to(shuf_dot, RIGHT, buff=0.25)
            self.play(FadeIn(real_val), FadeIn(shuf_val),
                      FadeIn(concl, shift=UP * 0.2), run_time=T.T_NORMAL)

    # ================================================================== #
    # 6 — REPRODUCIBILITY: the artifacts
    # ================================================================== #
    def reproducibility(self):
        head = txt("The code writes the paper", size=T.SUB_SIZE,
                   color=T.SIGNAL, weight="BOLD").to_edge(UP, buff=0.8)
        # one pipeline: code -> tables -> figures -> manuscript
        files = [("analysis.py", T.CLAUDE), ("results.csv", T.SIGNAL),
                 ("Figure 3", T.SIGNAL), ("manuscript.pdf", T.POSITIVE)]
        chain = VGroup(*[
            VGroup(RoundedRectangle(corner_radius=0.1, width=3.4, height=0.7,
                                    stroke_width=T.STROKE_MUTED, stroke_color=c,
                                    fill_color=c, fill_opacity=0.06),
                   txt(f, size=T.LABEL_SIZE, color=T.FG))
            for f, c in files
        ])
        for r in chain:
            r[1].move_to(r[0])
        chain.arrange(DOWN, buff=0.45).next_to(head, DOWN, buff=0.55).shift(LEFT * 2.6)
        arrows = VGroup(*[
            Arrow(chain[i].get_bottom(), chain[i + 1].get_top(), buff=0.06,
                  color=T.MUTED, stroke_width=T.STROKE_MUTED,
                  max_tip_length_to_length_ratio=0.28)
            for i in range(len(chain) - 1)
        ])
        chips = VGroup(
            EvidenceCard("tests", "17 automated", accent=T.POSITIVE, width=3.2, height=1.2),
            EvidenceCard("source", "open", accent=T.SIGNAL, width=3.2, height=1.2),
            EvidenceCard("make all", "≈ 8 s · laptop", accent=T.CLAUDE, width=3.2, height=1.2),
        ).arrange(DOWN, buff=0.4).to_edge(RIGHT, buff=1.1)
        with self.beat("None of this asks you to take our word for it."):
            self.crossfade_in(head)
        with self.beat("The same code that runs the analysis writes the "
                       "tables, the figures, and the manuscript itself."):
            self.play(LaggedStart(
                *[FadeIn(m) for pair in zip(chain, list(arrows) + [None])
                  for m in ([pair[0], pair[1]] if pair[1] else [pair[0]])],
                lag_ratio=0.35), run_time=2.2)
        with self.beat("Seventeen tests, open source, and the whole core "
                       "rebuilding in about eight seconds on a laptop."):
            self.play(LaggedStart(*[FadeIn(c, shift=LEFT * 0.3) for c in chips],
                                  lag_ratio=0.2), run_time=T.T_REVEAL)

    # ================================================================== #
    # 7 — CLOSING: the boundary of the claim
    # ================================================================== #
    def closing(self):
        left = txt("Recoverable structure", size=T.SUB_SIZE, color=T.SIGNAL,
                   weight="BOLD")
        neq = txt("≠", size=64, color=T.NEGATIVE, weight="BOLD")
        right = txt("Inductive predictability", size=T.SUB_SIZE, color=T.MUTED,
                    weight="BOLD")
        row = VGroup(left, neq, right).arrange(RIGHT, buff=0.6).move_to(UP * 0.6)
        thesis = VGroup(
            txt("Claude helped us find the structure —", size=T.BODY_SIZE,
                color=T.CLAUDE),
            txt("and the boundary of the claim.", size=T.BODY_SIZE,
                color=T.CLAUDE, weight="BOLD"),
        ).arrange(DOWN, buff=0.25).next_to(row, DOWN, buff=1.2)
        with self.beat("Which leaves one honest lesson. Recoverable structure "
                       "is not the same as inductive predictability."):
            self.crossfade_in(left)
            self.play(Write(neq), run_time=T.T_MICRO)
            self.play(FadeIn(right, shift=LEFT * 0.2), run_time=T.T_NORMAL)
            self.play(Indicate(neq, color=T.NEGATIVE, scale_factor=1.3),
                      run_time=T.T_NORMAL)
        with self.beat("Claude helped us find the structure — and, just as "
                       "importantly, where the claim stops."):
            self.play(FadeIn(thesis, shift=UP * 0.2), run_time=T.T_REVEAL)

        # end card
        endcard = VGroup(
            txt("The Empirical Regulatory Operator of CD4⁺ T Cells",
                size=T.BODY_SIZE, color=T.FG, weight="BOLD"),
            txt("Built with Claude · Life Sciences", size=T.LABEL_SIZE,
                color=T.CLAUDE),
            txt("github.com/cuentadesanti/cd4-perturbseq-regulator-atlas",
                size=T.SMALL_SIZE, color=T.MUTED),
        ).arrange(DOWN, buff=0.4)
        with self.beat("The empirical regulatory operator of CD4 T-cells — "
                       "built with Claude."):
            self.crossfade_in(endcard, run_time=1.0)
        self.wait(1.5)
