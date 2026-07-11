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
        # Reset between chapters; each chapter narrates itself in beats.
        if self.mobjects:
            self.play(*[FadeOut(m) for m in self.mobjects], run_time=0.4)
        self.camera.frame.move_to(self.base_center).set_width(self.base_width)
        fn()

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

        with self.beat("Meet the regulatory operator — a map of how one "
                       "human immune cell is wired."):
            head = title("The regulatory operator", color=T.SIGNAL)
            sub = txt("Genome-scale CRISPRi Perturb-seq · primary human "
                      "CD4⁺ T cells", size=T.SUB_SIZE, color=T.FG)
            VGroup(head, sub).arrange(DOWN, buff=0.4)
            self.play(Write(head), run_time=T.T_NORMAL)
            self.play(FadeIn(sub, shift=UP * 0.2), run_time=T.T_MICRO)
        self.play(FadeOut(sub), head.animate.scale(0.5).to_edge(UP),
                  run_time=0.6)

        # a row of genes; knock one down
        genes = VGroup(*[Dot(radius=0.13, color=T.SIGNAL, fill_opacity=0.45)
                         for _ in range(13)]).arrange(RIGHT, buff=0.4)
        genes.move_to(UP * 1.7)
        k = 6
        cross = txt("✕", size=30, color=T.NEGATIVE, weight="BOLD").move_to(genes[k])
        klab = txt("CRISPRi", size=T.SMALL_SIZE, color=T.NEGATIVE).next_to(
            genes[k], UP, buff=0.25)
        with self.beat("Take a single gene, and switch it off with CRISPR "
                       "interference."):
            self.play(LaggedStart(*[FadeIn(g) for g in genes], lag_ratio=0.05),
                      run_time=1.0)
            self.play(genes[k].animate.set_color(T.NEGATIVE).set_opacity(1.0),
                      FadeIn(cross), FadeIn(klab), run_time=0.8)

        # measure every other gene -> one response vector (the fingerprint)
        vals = rng.normal(0, 1, 16)
        vec = VGroup(*[
            Rectangle(width=cell, height=cell, stroke_width=0,
                      fill_color=(T.SIGNAL if v > 0 else T.NEGATIVE),
                      fill_opacity=float(min(abs(v), 1)) * 0.9)
            for v in vals
        ]).arrange(DOWN, buff=0.02).move_to(RIGHT * 3.6 + DOWN * 0.3)
        veclab = txt("response fingerprint\n(2,000 genes)", size=T.SMALL_SIZE,
                     color=T.MUTED).next_to(vec, RIGHT, buff=0.3)
        arr = Arrow(genes[k].get_bottom(), vec.get_top(), color=T.MUTED,
                    stroke_width=2, buff=0.2)
        with self.beat("Then measure how every other gene in the cell "
                       "responds. That readout — one number per gene — is the "
                       "perturbation's fingerprint."):
            self.play(Create(arr), run_time=0.6)
            self.play(LaggedStart(*[GrowFromCenter(s) for s in vec],
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
        with self.beat("Repeat for the next gene, and the next — thousands of "
                       "perturbations — and the fingerprints stack into one "
                       "regulator-by-gene map."):
            self.play(FadeOut(genes), FadeOut(cross), FadeOut(klab),
                      FadeOut(arr), FadeOut(veclab), run_time=0.4)
            self.play(ReplacementTransform(vec, colgroups[0]), run_time=0.8)
            self.play(LaggedStart(*[FadeIn(cg) for cg in colgroups[1:]],
                                  lag_ratio=0.05), run_time=1.6)
            self.play(FadeIn(mlab), run_time=0.4)

        # fold into the square operator
        with self.beat("Compare every pair of fingerprints, and the map folds "
                       "into one square object — the operator. Everything that "
                       "follows is asking this matrix questions."):
            op = self.make_operator_matrix(n=20, size=5.0).move_to(DOWN * 0.3)
            self.play(FadeTransform(mat, op), FadeOut(mlab), run_time=1.4)
            cap = txt("the operator — regulator × regulator", size=T.LABEL_SIZE,
                      color=T.SIGNAL).next_to(op, DOWN, buff=0.4)
            self.play(FadeIn(cap), run_time=0.5)

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
                   color=T.CLAUDE, weight="BOLD").to_edge(UP, buff=0.8)
        steps = ["Question", "Hypothesis", "Executable test",
                 "Artifact", "Adversarial audit"]
        rows = VGroup(*[
            VGroup(
                RoundedRectangle(corner_radius=0.12, width=4.6, height=0.8,
                                 stroke_width=T.STROKE_MUTED,
                                 stroke_color=T.CLAUDE,
                                 fill_color=T.CLAUDE, fill_opacity=0.06),
                txt(s, size=T.LABEL_SIZE, color=T.FG),
            )
            for s in steps
        ])
        for r in rows:
            r[1].move_to(r[0])
        rows.arrange(DOWN, buff=0.45).next_to(head, DOWN, buff=0.6)
        arrows = VGroup(*[
            Arrow(rows[i].get_bottom(), rows[i + 1].get_top(), buff=0.08,
                  color=T.CLAUDE, stroke_width=T.STROKE_MUTED,
                  max_tip_length_to_length_ratio=0.2)
            for i in range(len(rows) - 1)
        ])
        with self.beat("Before the findings, a word on method. Every claim in "
                       "this project ran the same loop."):
            self.play(FadeIn(head, shift=DOWN * 0.2), run_time=T.T_MICRO)
        with self.beat("A question becomes a hypothesis; the hypothesis "
                       "becomes an executable test; and whatever it produces "
                       "gets audited — adversarially."):
            self.play(
                LaggedStart(
                    *[FadeIn(m) for pair in zip(rows, list(arrows) + [None])
                      for m in ([pair[0], pair[1]] if pair[1] else [pair[0]])],
                    lag_ratio=0.4,
                ),
                run_time=2.6,
            )

        verdict = VGroup(
            txt("survives", size=T.LABEL_SIZE, color=T.POSITIVE, weight="BOLD"),
            txt("/", size=T.LABEL_SIZE, color=T.MUTED),
            txt("weakens", size=T.LABEL_SIZE, color=T.WARNING, weight="BOLD"),
            txt("/", size=T.LABEL_SIZE, color=T.MUTED),
            txt("fails", size=T.LABEL_SIZE, color=T.NEGATIVE, weight="BOLD"),
        ).arrange(RIGHT, buff=0.25).next_to(rows, DOWN, buff=0.5)
        with self.beat("Only then does a claim get to survive, weaken, or fail."):
            self.play(FadeIn(verdict, shift=UP * 0.2), run_time=T.T_NORMAL)

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

        with self.beat("Inside the operator, which directions are real signal "
                       "— and which are just noise?"):
            self.play(FadeIn(head, shift=DOWN * 0.2), run_time=T.T_MICRO)
            self.play(Create(ax), FadeIn(xlab), run_time=T.T_NORMAL)
            self.play(LaggedStart(*[GrowFromCenter(b) for b in bars],
                                  lag_ratio=0.01), run_time=1.6)

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
                        size=T.TITLE_SIZE, color=T.SIGNAL, weight="BOLD"
                        ).move_to(c_anchor))
        counter_static = txt("336", size=T.TITLE_SIZE, color=T.SIGNAL,
                             weight="BOLD").move_to(c_anchor)
        counter_lab = txt("signal eigenvalues", size=T.LABEL_SIZE,
                          color=T.SIGNAL).move_to(c_anchor + DOWN * 0.75)
        with self.beat("Random-matrix theory draws a clean line. The "
                       "closed-form Marchenko-Pastur edge counts three hundred "
                       "and thirty-six directions above the noise."):
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
        self.remove(counter_static)
        self.add(counter)
        with self.beat("But that edge is optimistic. A permutation null — "
                       "shuffling the data to see what noise alone produces — "
                       "pushes the threshold outward, and the honest count "
                       "drops to ninety-two."):
            self.play(
                LaggedStart(*[b.animate.set_fill(T.MUTED) for b in demoted],
                            lag_ratio=0.02),
                count_val.animate.set_value(92),
                Create(edge_e), FadeIn(edge_e_lab),
                run_time=T.T_REVEAL,
            )
            self.play(Indicate(counter_lab, color=T.SIGNAL, scale_factor=1.15),
                      run_time=T.T_NORMAL)

        note = txt("real structure — but only ~7 of these generalise across "
                   "conditions", size=T.SMALL_SIZE, color=T.MUTED)
        note.next_to(ax, UP, buff=0.2).to_edge(LEFT, buff=1.2)
        with self.beat("Real, reproducible structure. And of those ninety-two, "
                       "only about seven generalise to unseen conditions."):
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
                       "communities — eight of them; three are stable enough "
                       "to trust."):
            self.play(FadeIn(heat), FadeIn(head), run_time=T.T_NORMAL)

        # target: the compact strong module in the lower-right corner
        w_img = heat.width
        h_img = heat.height
        bx = heat.get_center()[0] + 0.42 * w_img
        by = heat.get_center()[1] - 0.36 * h_img
        block = [bx, by, 0]
        box = Rectangle(width=0.16 * w_img, height=0.16 * h_img,
                        stroke_color=T.POSITIVE, stroke_width=T.STROKE_MAIN)
        box.move_to(block)

        mask = spotlight(box, opacity=0.78)
        with self.beat("Watch this one."):
            self.play(FadeIn(mask), Create(box), FadeOut(head),
                      run_time=T.T_NORMAL)
            self.focus_on(box, width=box.width * 3.2, run_time=1.6)

        # the genes name themselves — vector text, crisp under zoom.
        # font sizes are small here because the camera is zoomed ~3x in.
        genes = ["NDUFA9", "NDUFS2", "NDUFV1", "NDUFS3", "NDUFB8", "NDUFS7"]
        gene_col = VGroup(*[
            txt(g, size=13, color=T.POSITIVE) for g in genes
        ]).arrange(DOWN, buff=0.11, aligned_edge=LEFT)
        gene_col.next_to(box, RIGHT, buff=0.3)
        with self.beat("Its members name themselves — the N-D-U-F genes, "
                       "subunits of Mitochondrial Respiratory Chain Complex One."):
            self.play(LaggedStart(*[FadeIn(g, shift=RIGHT * 0.1)
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
        stamp = txt("annotations added after clustering — the model was never told",
                    size=T.SMALL_SIZE, color=T.MUTED)
        group = VGroup(panel, ident).move_to(UP * 0.3)
        stamp.next_to(group, DOWN, buff=0.6)
        with self.beat("An external database, CORUM, confirms the identity at "
                       "a false-discovery rate near one in ten million. The "
                       "clustering was blind — the model was never told what "
                       "these genes were."):
            self.play(FadeOut(gene_col), run_time=0.3)
            self.reset_camera(run_time=1.2)
            self.play(FadeOut(mask), FadeOut(heat), FadeOut(box), run_time=0.5)
            self.play(FadeIn(panel), Write(ident), run_time=T.T_REVEAL)
            self.play(FadeIn(stamp), run_time=T.T_NORMAL)

    # ================================================================== #
    # 5 — SKEPTICISM: what did NOT survive
    # ================================================================== #
    def scientific_skepticism(self):
        head = txt("Claude argued against the result too", size=T.SUB_SIZE,
                   color=T.CLAUDE, weight="BOLD").to_edge(UP, buff=0.7)
        with self.beat("Good science also reports what didn't hold. Claude was "
                       "asked to argue against its own results."):
            self.play(FadeIn(head, shift=DOWN * 0.2), run_time=T.T_MICRO)

        gate = MetricGate("SAGA core module stability", value=0.56, gate=0.80,
                          lo=0.0, hi=1.0, width=7.0, accent=T.SIGNAL)
        gate.move_to(UP * 0.7)
        # build with marker at 0, then walk it to 0.56 (stops short of the gate)
        start_pos = gate.marker_position(0.0)
        gate.marker.move_to(start_pos)
        gate.val_lab.next_to(gate.marker, DOWN, buff=0.2)
        end_pos = gate.marker_position(0.56)
        verdicts = VGroup(
            ClaimVerdict("SAGA subunits co-cluster into one group", "recovered", ok=True),
            ClaimVerdict("a tight, gate-passing stable module", "below 0.80 gate", ok=False),
        ).arrange(DOWN, buff=0.35, aligned_edge=LEFT).next_to(gate, DOWN, buff=1.0)
        with self.beat("Take the SAGA chromatin module. Its subunits do cluster "
                       "together — but the group's stability, zero point five "
                       "six, sits below our zero point eight bar. So we flag "
                       "it, rather than claim it."):
            self.play(FadeIn(gate), run_time=T.T_NORMAL)
            self.play(
                gate.marker.animate.move_to(end_pos).set_color(T.WARNING),
                gate.val_lab.animate.move_to(end_pos + DOWN * 0.45).set_color(T.WARNING),
                run_time=T.T_REVEAL, rate_func=T.EASE_MOVE,
            )
            self.play(LaggedStart(*[FadeIn(v, shift=UP * 0.15) for v in verdicts],
                                  lag_ratio=0.4), run_time=T.T_REVEAL)

        # the inductive null: real ≈ shuffled ≈ 0
        q = txt("Can regulator features predict an unseen regulator's response?",
                size=T.LABEL_SIZE, color=T.FG).next_to(head, DOWN, buff=0.7)
        zero_axis = Line(LEFT * 4.5, RIGHT * 4.5, color=T.MUTED,
                         stroke_width=T.STROKE_MUTED).move_to(DOWN * 0.6)
        zero_tick = txt("held-out R² = 0", size=T.SMALL_SIZE,
                        color=T.MUTED).next_to(zero_axis.get_center(), UP, buff=0.15)
        real = VGroup(Dot(radius=0.12, color=T.SIGNAL),
                      txt("real features   R² ≈ −0.0005", size=T.SMALL_SIZE,
                          color=T.SIGNAL))
        real[0].move_to(zero_axis.get_center() + LEFT * 3.2)
        real[1].next_to(real[0], DOWN, buff=0.3)
        shuf = VGroup(Dot(radius=0.12, color=T.MUTED),
                      txt("shuffled features   R² ≈ −0.00003", size=T.SMALL_SIZE,
                          color=T.MUTED))
        shuf[0].move_to(zero_axis.get_center() + RIGHT * 3.2)
        shuf[1].next_to(shuf[0], UP, buff=0.3)
        with self.beat("And the hardest test: can a regulator's own features "
                       "predict the response of a regulator we have never "
                       "perturbed?"):
            self.play(FadeOut(gate), FadeOut(verdicts), run_time=0.5)
            self.play(FadeIn(q), run_time=T.T_MICRO)
            self.play(Create(zero_axis), FadeIn(zero_tick), run_time=T.T_NORMAL)
            self.play(FadeIn(real), FadeIn(shuf), run_time=T.T_NORMAL)

        concl = txt("real ≈ shuffled ≈ 0   —   a clean, reported null "
                    "(linear and non-linear)", size=T.LABEL_SIZE,
                    color=T.NEGATIVE, weight="BOLD").next_to(zero_axis, DOWN, buff=1.2)
        with self.beat("The answer is no. Real features and randomly shuffled "
                       "features both score essentially zero. A clean null — "
                       "reported, not buried."):
            self.play(
                real[0].animate.move_to(zero_axis.get_center() + LEFT * 0.15),
                shuf[0].animate.move_to(zero_axis.get_center() + RIGHT * 0.15),
                run_time=T.T_REVEAL, rate_func=T.EASE_MOVE,
            )
            self.play(FadeIn(concl, shift=UP * 0.2), run_time=T.T_NORMAL)

    # ================================================================== #
    # 6 — REPRODUCIBILITY: the artifacts
    # ================================================================== #
    def reproducibility(self):
        head = txt("Every claim ships as an artifact", size=T.SUB_SIZE,
                   color=T.SIGNAL, weight="BOLD").to_edge(UP, buff=0.8)
        cards = [
            ("manuscript", "peer-review draft", T.SIGNAL),
            ("tables", "versioned CSVs", T.SIGNAL),
            ("nulls", "matched permutation", T.WARNING),
            ("tests", "17 automated", T.POSITIVE),
            ("repo", "open source", T.SIGNAL),
            ("make all", "≈ 8 s · laptop", T.CLAUDE),
        ]
        cgroup = VGroup(*[
            EvidenceCard(cap, val, accent=acc, width=3.4, height=1.5)
            for cap, val, acc in cards
        ]).arrange_in_grid(rows=2, cols=3, buff=0.5).next_to(head, DOWN, buff=0.7)
        dirs = [LEFT, UP, RIGHT, LEFT, DOWN, RIGHT]
        with self.beat("None of this asks you to take our word for it."):
            self.play(FadeIn(head, shift=DOWN * 0.2), run_time=T.T_MICRO)
        with self.beat("Every claim ships as an artifact — a manuscript, "
                       "versioned tables, matched nulls, seventeen automated "
                       "tests, an open repository — the whole core rebuilding "
                       "in about eight seconds on a laptop."):
            self.play(
                LaggedStart(*[FadeIn(c, shift=-d * 0.4)
                              for c, d in zip(cgroup, dirs)], lag_ratio=0.15),
                run_time=2.4,
            )

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
            self.play(FadeIn(left, shift=RIGHT * 0.2), run_time=T.T_NORMAL)
            self.play(Write(neq), run_time=T.T_MICRO)
            self.play(FadeIn(right, shift=LEFT * 0.2), run_time=T.T_NORMAL)
            self.play(Indicate(neq, color=T.NEGATIVE, scale_factor=1.3),
                      run_time=T.T_NORMAL)
        with self.beat("Claude helped us find the structure — and, just as "
                       "importantly, where the claim stops."):
            self.play(FadeIn(thesis, shift=UP * 0.2), run_time=T.T_REVEAL)
        self.play(FadeOut(row), FadeOut(thesis), run_time=1.0)
