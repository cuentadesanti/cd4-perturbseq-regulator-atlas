"""The regulatory operator — a biology-led ~2:50 visual essay.

One MovingCameraScene. The story is the discoveries; rigor is felt through a
few decisive controls, not narrated as separate chapters. A single protagonist
(fingerprint → operator → denoised operator → community matrix → Complex I)
stays alive across the opening scenes instead of fading to black between them.

Render:
    manim -pql demo.py FullDemo                 # fast preview
    NARRATE=say ./render.sh final say           # 1080p + narration
    BURN_SUBS=1 NARRATE=say ./render.sh final say  # + burned-in subtitles
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
    Rectangle,
    RoundedRectangle,
    Circle,
    Line,
    DashedLine,
    Arrow,
    DoubleArrow,
    ArcBetweenPoints,
    Dot,
    FadeIn,
    FadeOut,
    Write,
    Create,
    GrowFromCenter,
    GrowArrow,
    Wiggle,
    Transform,
    ReplacementTransform,
    TransformFromCopy,
    FadeTransform,
    LaggedStart,
    Indicate,
    rate_functions,
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

        self.narrate = os.environ.get("NARRATE", "").lower()
        if self.narrate == "say":
            from voice import SayService
            self.set_speech_service(
                SayService(voice=os.environ.get("SAY_VOICE", "Samantha"),
                           rate=int(os.environ.get("SAY_RATE", "138")),
                           ffmpeg=os.environ.get("FFMPEG", "ffmpeg")))
        elif self.narrate == "deepgram":
            from voice import DeepgramService
            self.set_speech_service(DeepgramService(
                voice=os.environ.get("DG_VOICE", "aura-asteria-en")))
        elif self.narrate == "eleven":
            from voice import ElevenLabsByID
            self.set_speech_service(
                ElevenLabsByID(
                    voice_id=os.environ.get("ELEVEN_VOICE_ID",
                                            "EST9Ui6982FZPSi7gCHi"),
                    model="eleven_turbo_v2_5"))

        for name, fn in [
            ("operator", self.scene_operator),
            ("denoise", self.scene_denoise),
            ("complex_i", self.scene_complex_i),
            ("saga", self.scene_saga),
            ("controls", self.scene_controls),
            ("repro", self.scene_repro),
            ("close", self.scene_close),
        ]:
            self.next_section(name)  # render organization only — not a reset
            self.camera.frame.move_to(self.base_center).set_width(self.base_width)
            fn()

    # ---- helpers ----------------------------------------------------- #
    @contextmanager
    def beat(self, text):
        if self.narrate in ("say", "eleven", "deepgram"):
            with self.voiceover(text=text) as tr:
                yield tr
        else:
            yield None
            self.wait(0.4)

    def crossfade_in(self, *news, run_time=0.7, keep=()):
        keep = list(keep)
        outs = [FadeOut(m) for m in self.mobjects if m not in keep]
        ins = [FadeIn(n) for n in news]
        if outs or ins:
            self.play(*outs, *ins, run_time=run_time)

    def focus_on(self, mob, width=None, run_time=1.4):
        self.play(
            self.camera.frame.animate.move_to(mob).set(
                width=width or mob.width * 1.4),
            run_time=run_time, rate_func=T.EASE_REVEAL)

    def reset_camera(self, run_time=1.2):
        self.play(
            self.camera.frame.animate.move_to(self.base_center).set_width(
                self.base_width),
            run_time=run_time, rate_func=T.EASE_MOVE)

    def make_operator_matrix(self, n=20, size=5.0, seed=7):
        rng = np.random.default_rng(seed)
        blocks = [(0, 6), (6, 12), (12, 16)]
        M = rng.normal(0, 0.25, (n, n))
        for a, b in blocks:
            b = min(b, n)
            if b - a < 2:
                continue
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
                sq = Rectangle(
                    width=cell, height=cell, stroke_width=0,
                    fill_color=(T.SIGNAL if val > 0 else T.NEGATIVE),
                    fill_opacity=float(min(abs(val), 1)) * 0.9)
                sq.move_to([(j - n / 2 + 0.5) * cell,
                            (n / 2 - i - 0.5) * cell, 0])
                grid.add(sq)
        return grid

    def crispri_intro(self, rng):
        """~8 s schematic CRISPRi diagram with two distinct transcription
        cycles: (1) one successful pass, then (2) a second polymerase that is
        blocked by a dCas9–KRAB complex parked on the promoter. Transcription
        drops qualitatively (not a literal %), then the regulator hands off to
        three downstream genes (↑ / ↓ / ≈). Returns the end-state nodes."""
        S, C, M = T.SIGNAL, T.CLAUDE, T.MUTED
        dna_y, pol_y = 0.0, 0.40  # Pol travels above the upper DNA strand
        gene_c = np.array([0.0, dna_y - 0.09, 0])
        gene_w, gene_h = 3.0, 0.42
        gene_left = gene_c[0] - gene_w / 2
        gene_bot = gene_c[1] - gene_h / 2  # -0.30

        def make_polymerase():
            body = RoundedRectangle(width=0.66, height=0.44, corner_radius=0.22,
                                    fill_opacity=0.92, stroke_width=0, fill_color=S)
            return VGroup(body, txt("Pol", size=11, color=T.BG).move_to(body))

        # genomic layout: promoter → gap → TSS arrow → gap → gene body
        dna_top = Line([-4.3, dna_y, 0], [4.3, dna_y, 0], stroke_width=3, color=M)
        dna_bot = dna_top.copy().shift(DOWN * 0.18)
        promoter = Line([-2.95, dna_y, 0], [-2.25, dna_y, 0], stroke_width=9, color=S)
        prom_lab = txt("promoter", size=11, color=M).next_to(promoter, DOWN, buff=0.24)
        tss = Arrow([-2.02, dna_y, 0], [-1.68, dna_y, 0], buff=0, stroke_width=2.2,
                    color=M, max_tip_length_to_length_ratio=0.35)
        gene = RoundedRectangle(width=gene_w, height=gene_h, corner_radius=0.12,
                                stroke_width=2, stroke_color=S, fill_color=S,
                                fill_opacity=0.12).move_to(gene_c)
        gene_lab = txt("target regulator", size=17, color=M).set_opacity(0.8)
        gene_lab.next_to(gene, UP, buff=0.80)  # clears Pol II with margin

        self.play(Create(dna_top), Create(dna_bot), run_time=0.5)
        self.play(FadeIn(promoter), FadeIn(prom_lab), Create(tss),
                  FadeIn(gene), FadeIn(gene_lab), run_time=0.7)

        # (1) one successful pass — a convoy of polymerases making several RNAs
        pol_a = make_polymerase().move_to([-4.7, pol_y, 0])
        pol_b = make_polymerase().move_to([-6.0, pol_y, 0])
        transcripts = VGroup()
        for frac, L in zip((0.20, 0.36, 0.52, 0.68, 0.84),
                           (1.0, 0.9, 0.8, 0.72, 0.64)):
            x = gene_left + frac * gene_w
            stem = Line([x, gene_bot, 0], [x, gene_bot - 0.16, 0],
                        stroke_width=3, color=S)
            arc = ArcBetweenPoints([x, gene_bot - 0.16, 0],
                                   [x + 0.42 * L, gene_bot - 0.16 - 0.40 * L, 0],
                                   angle=-0.6, stroke_width=3, color=S)
            transcripts.add(VGroup(stem, arc))
        self.play(FadeIn(pol_a, shift=RIGHT * 0.2), run_time=0.35)
        self.play(pol_a.animate.move_to([1.6, pol_y, 0]),
                  FadeIn(pol_b, shift=RIGHT * 0.2),
                  LaggedStart(*[Create(t) for t in transcripts], lag_ratio=0.35),
                  run_time=2.3, rate_func=rate_functions.linear)
        self.play(pol_a.animate.shift(RIGHT * 1.1).set_opacity(0.0),
                  pol_b.animate.move_to([1.6, pol_y, 0]),
                  run_time=0.9, rate_func=rate_functions.linear)
        self.play(pol_b.animate.shift(RIGHT * 1.1).set_opacity(0.0), run_time=0.5)
        self.wait(0.4)  # hold after successful transcription
        self.remove(pol_a, pol_b)

        # (2) recruit dCas9–KRAB to the promoter (Claude orange; no DNA cut)
        guide = ArcBetweenPoints([-2.6, 1.8, 0], promoter.get_center() + UP * 0.04,
                                 angle=-0.9, stroke_width=2.5, color=C)
        dcas9 = RoundedRectangle(width=0.86, height=0.58, corner_radius=0.24,
                                 fill_opacity=0.95, stroke_width=0, fill_color=C)
        krab = RoundedRectangle(width=0.4, height=0.32, corner_radius=0.13,
                                fill_opacity=0.95, stroke_width=0, fill_color=C)
        krab.next_to(dcas9, UP + RIGHT, buff=-0.14)
        crispri = VGroup(dcas9, krab,
                         txt("CRISPRi", size=10, color=T.BG).move_to(dcas9))
        crispri.move_to([-2.6, 1.6, 0])
        self.play(Create(guide), FadeIn(crispri, shift=DOWN * 0.2), run_time=0.70)
        self.play(crispri.animate.move_to([-2.6, pol_y, 0]),
                  promoter.animate.set_color(C), run_time=0.85,
                  rate_func=T.EASE_MOVE)
        self.play(FadeOut(guide), run_time=0.25)
        self.wait(0.35)  # hold after binding

        # (2b) a second polymerase enters and stalls just before the complex
        pol_blocked = make_polymerase().move_to([-4.5, pol_y, 0])
        self.play(FadeIn(pol_blocked, shift=RIGHT * 0.2), run_time=0.35)
        self.play(pol_blocked.animate.move_to([-3.5, pol_y, 0]), run_time=0.90,
                  rate_func=T.EASE_MOVE)
        self.play(pol_blocked.animate.shift(RIGHT * 0.15), run_time=0.40,
                  rate_func=rate_functions.there_and_back)
        self.play(pol_blocked.animate.set_opacity(0.0), run_time=0.40)
        self.remove(pol_blocked)

        # (2c) reduced transcription, shown qualitatively
        red_lab = txt("transcription reduced", size=T.SMALL_SIZE,
                      color=C).to_edge(UP, buff=0.9)
        self.play(LaggedStart(*[t.animate.set_opacity(0.1) for t in transcripts],
                              lag_ratio=0.3), FadeIn(red_lab),
                  gene.animate.set_fill(M, opacity=0.1), run_time=1.25)
        self.wait(0.25)  # hold reduced state

        # hand off: regulator node → downstream genes (↑ / ↓ / ≈)
        reg_node = VGroup(Circle(radius=0.33, fill_opacity=0.9, stroke_width=0,
                                 fill_color=M), txt("R", size=T.SMALL_SIZE, color=T.BG))
        reg_node[1].move_to(reg_node[0])
        reg_node.move_to([-4.2, -0.1, 0])
        molecular = VGroup(dna_top, dna_bot, promoter, prom_lab, tss, gene,
                           gene_lab, crispri, transcripts, red_lab)
        self.play(FadeOut(molecular), FadeIn(reg_node), run_time=0.6)

        down_nodes = VGroup()
        down_marks = VGroup()
        for y, col, mk in zip((1.3, 0.0, -1.3),
                              (T.POSITIVE, T.NEGATIVE, M), ("↑", "↓", "≈")):
            node = Circle(radius=0.26, fill_opacity=0.9, stroke_width=0,
                          fill_color=col).move_to([2.6, y, 0])
            down_nodes.add(node)
            down_marks.add(txt(mk, size=T.SUB_SIZE, color=col,
                               weight="BOLD").next_to(node, RIGHT, buff=0.25))
        arrows = VGroup(*[
            Arrow(reg_node.get_right(), n.get_left(), buff=0.12, stroke_width=2,
                  color=M) for n in down_nodes])
        self.play(LaggedStart(*[GrowArrow(a) for a in arrows], lag_ratio=0.15),
                  FadeIn(down_nodes), run_time=0.7)
        self.play(LaggedStart(*[FadeIn(m, shift=RIGHT * 0.1) for m in down_marks],
                              lag_ratio=0.2), run_time=0.55)
        return reg_node, down_nodes, down_marks, arrows

    # ================================================================== #
    # 1 — BIOLOGICAL QUESTION + OPERATOR CONSTRUCTION
    # ================================================================== #
    def scene_operator(self):
        rng = np.random.default_rng(3)
        cell = 0.26

        q = VGroup(
            txt("Do thousands of gene perturbations form", size=T.SUB_SIZE,
                color=T.FG),
            txt("a hidden regulatory architecture?", size=T.TITLE_SIZE,
                color=T.SIGNAL, weight="BOLD"),
            txt("— or just a noisy list of hits?", size=T.BODY_SIZE,
                color=T.MUTED),
        ).arrange(DOWN, buff=0.35)
        with self.beat("Do thousands of gene perturbations form a hidden "
                       "regulatory architecture — or are they just a noisy "
                       "list of hits?"):
            self.play(FadeIn(q[0], shift=UP * 0.1), run_time=T.T_MICRO)
            self.play(Write(q[1]), run_time=T.T_NORMAL)
            self.play(FadeIn(q[2], shift=UP * 0.1), run_time=T.T_MICRO)
            self.wait(0.8)
        self.play(FadeOut(q), run_time=0.5)

        # --- CRISPRi molecular animation → downstream genes ---
        with self.beat("CRISPR interference reduces the expression of one "
                       "regulator without cutting the DNA. The screen then "
                       "measures how the rest of the transcriptome responds."):
            reg_node, down_nodes, down_marks, arrows = self.crispri_intro(rng)

        # --- one downstream gene's cell observations → a standardized effect ---
        # the z-score is shown as the horizontal shift between the two means.
        def cloud(cx, cy, color, n=16):
            return VGroup(*[
                Dot(radius=0.05, color=color).move_to(
                    [cx + rng.normal(0, 0.32), cy + rng.uniform(-0.28, 0.28), 0])
                for _ in range(n)])
        gb = down_nodes[0]  # the up-regulated downstream gene
        axis_y = -1.7
        ctrl_cx, shift = 0.1, 2.2
        pert_cx = ctrl_cx + shift
        axis = Line([-1.6, axis_y, 0], [4.2, axis_y, 0], color=T.MUTED, stroke_width=2)
        axlab = txt("expression of one downstream gene", size=T.SMALL_SIZE,
                    color=T.MUTED).next_to(axis, DOWN, buff=0.15)
        ctrl = cloud(ctrl_cx, 0.45, T.MUTED)
        pert = cloud(ctrl_cx, -0.5, T.SIGNAL)  # starts aligned with control
        ctrl_lab = txt("control", size=T.SMALL_SIZE, color=T.MUTED).next_to(
            ctrl, LEFT, buff=0.4)
        pert_lab = txt("perturbed", size=T.SMALL_SIZE, color=T.SIGNAL).next_to(
            pert, LEFT, buff=0.4)
        mean_ctrl = DashedLine([ctrl_cx, 0.78, 0], [ctrl_cx, axis_y, 0],
                               color=T.MUTED, stroke_width=2)
        mean_pert = DashedLine([pert_cx, -0.2, 0], [pert_cx, axis_y, 0],
                               color=T.SIGNAL, stroke_width=2)
        delta = DoubleArrow([ctrl_cx, axis_y + 0.32, 0], [pert_cx, axis_y + 0.32, 0],
                            color=T.POSITIVE, stroke_width=3, buff=0, tip_length=0.16)
        dlab = txt("horizontal shift  (Δ)", size=T.SMALL_SIZE,
                   color=T.POSITIVE).next_to(delta, UP, buff=0.08)
        zval = txt("z = +2.4", size=T.SUB_SIZE, color=T.POSITIVE, weight="BOLD")
        with self.beat("Those single-cell observations become one signed, "
                       "standardized effect — the horizontal shift between the "
                       "means, measured in standard deviations."):
            self.play(FadeOut(reg_node), FadeOut(down_nodes[1]),
                      FadeOut(down_nodes[2]), FadeOut(down_marks),
                      FadeOut(arrows), run_time=0.4)
            self.play(FadeTransform(gb, pert), FadeIn(ctrl), Create(axis),
                      FadeIn(axlab), FadeIn(ctrl_lab), FadeIn(pert_lab),
                      run_time=T.T_NORMAL)
            self.wait(0.3)
            # the perturbed distribution shifts clearly to the right
            self.play(pert.animate.shift(RIGHT * shift), run_time=0.9,
                      rate_func=T.EASE_MOVE)
            # mark each mean and the horizontal difference on the axis
            self.play(Create(mean_ctrl), Create(mean_pert), run_time=0.55)
            self.play(GrowFromCenter(delta), FadeIn(dlab), run_time=0.55)
            self.wait(0.35)
            # that horizontal shift IS the standardized effect
            zval.next_to(delta, UP, buff=0.12)
            self.play(ReplacementTransform(dlab, zval), run_time=0.7)
            self.wait(0.45)

        # fingerprint
        vals = rng.normal(0, 1, 16)
        vals[0] = 2.4
        vec = VGroup(*[
            Rectangle(width=cell, height=cell, stroke_width=0,
                      fill_color=(T.POSITIVE if i == 0 else
                                  (T.SIGNAL if v > 0 else T.NEGATIVE)),
                      fill_opacity=float(min(abs(v), 1)) * 0.9)
            for i, v in enumerate(vals)
        ]).arrange(DOWN, buff=0.02).move_to(RIGHT * 3.7 + DOWN * 0.3)
        veclab = txt("transcriptional\nfingerprint", size=T.SMALL_SIZE,
                     color=T.MUTED).next_to(vec, RIGHT, buff=0.3)
        with self.beat("Every gene gives one; together they form a "
                       "transcriptional fingerprint."):
            self.play(*[FadeOut(m) for m in self.mobjects if m is not zval],
                      run_time=0.4)
            self.play(ReplacementTransform(zval, vec[0]), run_time=0.6)
            self.play(LaggedStart(*[GrowFromCenter(s) for s in vec[1:]],
                                  lag_ratio=0.05), FadeIn(veclab), run_time=1.2)

        # many fingerprints -> matrix -> operator
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
        with self.beat("Thousands of fingerprints form the regulatory "
                       "operator."):
            self.play(FadeOut(veclab), run_time=0.3)
            self.play(ReplacementTransform(vec, colgroups[0]), run_time=0.7)
            self.play(LaggedStart(*[FadeIn(cg) for cg in colgroups[1:]],
                                  lag_ratio=0.05), run_time=1.4)
            self.play(FadeIn(mlab), run_time=0.3)
            op = self.make_operator_matrix(n=20, size=5.0).move_to(DOWN * 0.3)
            self.play(FadeTransform(mat, op), FadeOut(mlab), run_time=1.2)
            cap = txt("the operator — regulator × regulator", size=T.LABEL_SIZE,
                      color=T.SIGNAL).next_to(op, DOWN, buff=0.4)
            self.play(FadeIn(cap), run_time=0.4)
        self.operator = op
        self.op_cap = cap

    # ================================================================== #
    # 2 — DENOISING, ONE IDEA (protagonist morphs into the community matrix)
    # ================================================================== #
    def scene_denoise(self):
        op = self.operator
        rng = np.random.default_rng(5)
        cells = list(op)
        noise = [c for c in cells if rng.random() < 0.7]

        cap = txt("92 empirically supported signal directions",
                  size=T.LABEL_SIZE, color=T.SIGNAL, weight="BOLD")
        sub = txt("empirical permutation null", size=T.SMALL_SIZE, color=T.MUTED)
        VGroup(cap, sub).arrange(DOWN, buff=0.2).to_edge(UP, buff=0.9)

        with self.beat("Most of that apparent structure was noise. An "
                       "empirical permutation control kept only the "
                       "reproducible signal —"):
            self.play(FadeOut(self.op_cap), run_time=0.3)
            self.play(LaggedStart(*[c.animate.set_opacity(0.05) for c in noise],
                                  lag_ratio=0.002), run_time=1.6)
            self.play(FadeIn(cap, shift=DOWN * 0.15), FadeIn(sub), run_time=0.6)

        heat = ImageMobject(str(ASSETS / "heatmap_panel.png"))
        heat.height = 6.2
        heat.move_to(DOWN * 0.2)
        with self.beat("the surviving directions that reorganized into the "
                       "communities we analyzed next."):
            self.play(FadeOut(cap), FadeOut(sub),
                      FadeTransform(op, heat), run_time=1.4)
        self.community_heat = heat

    # ================================================================== #
    # 3 — COMPLEX I  (the climax — longest scene)
    # ================================================================== #
    def scene_complex_i(self):
        heat = self.community_heat
        head = txt("community structure   z = 259 vs matched null",
                   size=22, color=T.MUTED).to_edge(UP, buff=0.7)
        with self.beat("Block-ordered, the operator's community structure sits "
                       "far above a matched null."):
            self.play(FadeIn(head, shift=DOWN * 0.15), run_time=T.T_NORMAL)

        w_img, h_img = heat.width, heat.height
        bx = heat.get_center()[0] + 0.44 * w_img
        by = heat.get_center()[1] - 0.44 * h_img
        box = Rectangle(width=0.12 * w_img, height=0.12 * h_img,
                        stroke_color=T.POSITIVE, stroke_width=T.STROKE_MAIN)
        box.move_to([bx, by, 0])
        genes = ["NDUFA9", "NDUFS2", "NDUFV1", "NDUFS3", "NDUFB8", "NDUFS7"]
        gene_col = VGroup(*[txt(g, size=13, color=T.POSITIVE) for g in genes]
                          ).arrange(DOWN, buff=0.1, aligned_edge=RIGHT)
        gene_col.next_to(box, LEFT, buff=0.4)
        focus_grp = Group(box, gene_col)

        mask = spotlight(box, opacity=0.6)
        with self.beat("Watch one community."):
            self.play(FadeIn(mask), Create(box), FadeOut(head), run_time=0.7)
            self.focus_on(focus_grp, width=focus_grp.width * 1.4, run_time=1.8)
            self.wait(0.5)

        with self.beat("Its members name themselves — the NDUF genes of "
                       "Mitochondrial Complex One."):
            self.play(LaggedStart(*[FadeIn(g, shift=LEFT * 0.1)
                                    for g in gene_col],
                                  lag_ratio=0.3), run_time=2.2)
            self.wait(0.6)

        ident = VGroup(
            txt("recovered stable community · n = 87 · 45% donor-robust",
                size=T.LABEL_SIZE, color=T.FG),
            txt("Mitochondrial Respiratory Chain Complex I",
                size=T.BODY_SIZE, color=T.POSITIVE, weight="BOLD"),
            txt("BH-FDR 1.4 × 10⁻⁷ · 8/13 CORUM holoenzyme members · "
                "10 NDUF-family genes", size=T.LABEL_SIZE, color=T.FG),
        ).arrange(DOWN, buff=0.25)
        panel = RoundedRectangle(
            corner_radius=0.15, width=ident.width + 1.0, height=ident.height + 0.7,
            stroke_width=T.STROKE_MUTED, stroke_color=T.POSITIVE,
            fill_color=T.FAINT, fill_opacity=0.5).move_to(ident)
        stamp = txt("BIOLOGICAL ANNOTATIONS ADDED AFTER CLUSTERING",
                    size=T.SMALL_SIZE, color=T.WARNING)
        group = VGroup(panel, ident).move_to(UP * 0.3)
        stamp.next_to(group, DOWN, buff=0.6)
        with self.beat("The clustering used no biological annotations. Only "
                       "after the community was fixed did we ask CORUM what it "
                       "contained — an eight-in-thirteen match at a false-"
                       "discovery rate near one in ten million."):
            self.play(
                FadeOut(gene_col), FadeOut(mask), FadeOut(heat), FadeOut(box),
                self.camera.frame.animate.move_to(self.base_center).set_width(
                    self.base_width),
                run_time=1.1, rate_func=T.EASE_MOVE)
            self.play(FadeIn(panel), run_time=0.5)
            self.play(LaggedStart(*[FadeIn(line, shift=UP * 0.1)
                                    for line in ident], lag_ratio=0.3),
                      run_time=T.T_REVEAL)
            self.play(FadeIn(stamp), run_time=T.T_NORMAL)
        self.wait(2.3)  # let the main discovery land

    # ================================================================== #
    # 4 — SAGA: a second, convergent biological result
    # ================================================================== #
    def scene_saga(self):
        chain = EvidenceChain(
            ["CORUM SAGA identity", "CP-factor concentration",
             "K562 module cohesion (z = 16)", "cross-cell-type concordance"],
            "convergent SAGA-centered module", accent=T.SIGNAL)
        chain.scale(0.9).move_to(ORIGIN)
        with self.beat("A second module, centered on SAGA, was backed by "
                       "independent evidence — complex identity, tensor "
                       "factors, and an external K562 screen."):
            self.crossfade_in(chain, run_time=0.7)

        # cell-type structure: SAGA universal, immune cytoskeleton T-specific
        head = txt("conserved vs. cell-type-specific", size=T.LABEL_SIZE,
                   color=T.MUTED).to_edge(UP, buff=0.8)
        colx = {"CD4⁺ T": -1.2, "K562": 1.2, "RPE1": 3.6}
        col_labels = VGroup(*[
            txt(k, size=T.SMALL_SIZE, color=T.MUTED).move_to([x, 1.7, 0])
            for k, x in colx.items()])
        saga_lab = txt("SAGA  (SUPT20H)", size=T.LABEL_SIZE, color=T.SIGNAL).move_to([-4.6, 0.7, 0])
        rac_lab = txt("Rac→WAVE→actin  (DOCK2)", size=T.LABEL_SIZE,
                      color=T.WARNING).move_to([-4.6, -0.9, 0])
        saga_dots = VGroup(*[Dot(radius=0.14, color=T.SIGNAL).move_to([x, 0.7, 0])
                             for x in colx.values()])
        # T-specific: present only in CD4
        rac_dots = VGroup(
            Dot(radius=0.14, color=T.WARNING).move_to([colx["CD4⁺ T"], -0.9, 0]),
            Dot(radius=0.09, color=T.MUTED, fill_opacity=0.3).move_to([colx["K562"], -0.9, 0]),
            Dot(radius=0.09, color=T.MUTED, fill_opacity=0.3).move_to([colx["RPE1"], -0.9, 0]))
        tags = VGroup(
            txt("conserved", size=T.SMALL_SIZE, color=T.SIGNAL).next_to(saga_dots, RIGHT, buff=0.5),
            txt("T-cell-specific", size=T.SMALL_SIZE, color=T.WARNING).next_to(rac_dots, RIGHT, buff=0.5))
        with self.beat("Across CD4, K562 and RPE1 cells, SAGA behaved as "
                       "broadly conserved machinery, while an immune "
                       "cytoskeletal program stayed T-cell-specific."):
            self.crossfade_in(head, col_labels, saga_lab, rac_lab, run_time=0.7)
            self.play(LaggedStart(*[GrowFromCenter(d) for d in saga_dots],
                                  lag_ratio=0.15), FadeIn(tags[0]), run_time=1.0)
            self.play(LaggedStart(*[GrowFromCenter(d) for d in rac_dots],
                                  lag_ratio=0.15), FadeIn(tags[1]), run_time=1.0)

        # disease nomination vs matched null (Complex I as control)
        d1 = ClaimVerdict("SAGA module — autoimmune-risk enriched",
                          "nominated", ok=True)
        d2 = ClaimVerdict("Complex I — not enriched", "control ✓", ok=False)
        dis = VGroup(d1, d2).arrange(DOWN, buff=0.35, aligned_edge=LEFT)
        disnote = txt("length-matched, MHC-excluded null · enrichment, not a "
                      "causal target", size=T.SMALL_SIZE, color=T.MUTED)
        block = VGroup(dis, disnote).arrange(DOWN, buff=0.5).move_to(ORIGIN)
        with self.beat("The SAGA module also carried autoimmune-risk signal "
                       "under a matched null — a nomination, not a cause — "
                       "while Complex One did not, a clean control."):
            self.crossfade_in(block, run_time=0.7)
        self.wait(1.0)

    # ================================================================== #
    # 5 — RIGOR THROUGH TWO DECISIVE CONTROLS
    # ================================================================== #
    def scene_controls(self):
        head = txt("Claude helped decide where the claims had to stop",
                   size=T.SUB_SIZE, color=T.CLAUDE, weight="BOLD").to_edge(
            UP, buff=0.7)

        # Control A — SAGA stability gate
        gate = MetricGate("SAGA partition stability", value=0.56, gate=0.80,
                          lo=0.0, hi=1.0, width=7.0, accent=T.SIGNAL)
        gate.move_to(UP * 0.9)
        gate.marker.move_to(gate.marker_position(0.0))
        gate.val_lab.next_to(gate.marker, DOWN, buff=0.2)
        end_pos = gate.marker_position(0.56)
        verdicts = VGroup(
            ClaimVerdict("SAGA subunits co-cluster", "observed", ok=True),
            ClaimVerdict("Stable recovered Leiden community", "rejected", ok=False),
            ClaimVerdict("Convergent regulatory module", "supported", ok=True),
        ).arrange(DOWN, buff=0.28, aligned_edge=LEFT).next_to(gate, DOWN, buff=0.8)
        with self.beat("Claude also challenged the appealing readings. SAGA's "
                       "convergent support was strong — but its partition "
                       "stability, 0.56, fell below our 0.8 gate. So: not a "
                       "recovered stable community."):
            self.crossfade_in(head, gate, run_time=0.7)
            self.play(
                gate.marker.animate.move_to(end_pos).set_color(T.WARNING),
                gate.val_lab.animate.move_to(end_pos + DOWN * 0.45).set_color(T.WARNING),
                run_time=T.T_REVEAL, rate_func=T.EASE_MOVE)
            self.play(LaggedStart(*[FadeIn(v, shift=UP * 0.15) for v in verdicts],
                                  lag_ratio=0.4), run_time=T.T_REVEAL)

        # Control B — unseen-regulator prediction (R² axis centred on zero)
        q = txt("predict a completely unseen regulator from side-features?",
                size=T.LABEL_SIZE, color=T.FG).next_to(head, DOWN, buff=0.55)
        zero = DashedLine(UP * 1.0, DOWN * 1.0, color=T.MUTED,
                          stroke_width=2).move_to(DOWN * 0.4)
        x0 = zero.get_center()[0]
        zlab = txt("held-out R² = 0", size=T.SMALL_SIZE, color=T.MUTED).next_to(
            zero, UP, buff=0.15)
        ay, by_ = zero.get_center()[1] + 0.4, zero.get_center()[1] - 0.4
        real_lab = txt("real features", size=T.SMALL_SIZE, color=T.SIGNAL).move_to([-4.3, ay, 0])
        shuf_lab = txt("shuffled features", size=T.SMALL_SIZE, color=T.MUTED).move_to([-4.3, by_, 0])
        real_dot = Dot(radius=0.11, color=T.SIGNAL).next_to(real_lab, RIGHT, buff=0.3)
        shuf_dot = Dot(radius=0.11, color=T.MUTED).next_to(shuf_lab, RIGHT, buff=0.3)
        concl = VGroup(
            txt("real ≈ shuffled ≈ 0", size=T.LABEL_SIZE, color=T.NEGATIVE, weight="BOLD"),
            txt("no detectable inductive signal", size=T.SMALL_SIZE, color=T.MUTED),
        ).arrange(DOWN, buff=0.15).next_to(zero, DOWN, buff=0.9)
        with self.beat("When we asked whether side-features could predict "
                       "completely unseen regulators, real and shuffled "
                       "features scored identically. Claude helped us stop "
                       "there and report the null."):
            self.play(FadeOut(gate), FadeOut(verdicts), run_time=0.5)
            self.play(FadeIn(q), Create(zero), FadeIn(zlab),
                      FadeIn(real_lab), FadeIn(shuf_lab), run_time=T.T_NORMAL)
            self.play(real_dot.animate.move_to([x0 - 0.08, ay, 0]),
                      shuf_dot.animate.move_to([x0 - 0.02, by_, 0]),
                      run_time=T.T_REVEAL, rate_func=T.EASE_MOVE)
            self.add(real_dot, shuf_dot)
            self.play(FadeIn(concl, shift=UP * 0.15), run_time=T.T_NORMAL)
            self.wait(0.4)

    # ================================================================== #
    # 6 — REPRODUCIBILITY, a short seal
    # ================================================================== #
    def scene_repro(self):
        files = [("analysis code", T.CLAUDE), ("versioned table", T.SIGNAL),
                 ("Figure 3", T.SIGNAL), ("manuscript", T.POSITIVE)]
        chain = VGroup(*[
            VGroup(RoundedRectangle(corner_radius=0.1, width=3.4, height=0.7,
                                    stroke_width=T.STROKE_MUTED, stroke_color=c,
                                    fill_color=c, fill_opacity=0.06),
                   txt(f, size=T.LABEL_SIZE, color=T.FG))
            for f, c in files])
        for r in chain:
            r[1].move_to(r[0])
        chain.arrange(DOWN, buff=0.4).shift(LEFT * 2.6)
        arrows = VGroup(*[
            Arrow(chain[i].get_bottom(), chain[i + 1].get_top(), buff=0.06,
                  color=T.MUTED, stroke_width=T.STROKE_MUTED,
                  max_tip_length_to_length_ratio=0.28)
            for i in range(len(chain) - 1)])
        chips = VGroup(
            EvidenceCard("data", "public", accent=T.SIGNAL, width=3.2, height=1.1),
            EvidenceCard("pipeline", "laptop · no 1.8 TB", accent=T.CLAUDE, width=3.2, height=1.1),
            EvidenceCard("controls", "matched nulls · tests", accent=T.POSITIVE, width=3.2, height=1.1),
        ).arrange(DOWN, buff=0.35).to_edge(RIGHT, buff=1.0)
        with self.beat("Every result here is backed by an executable analysis "
                       "and a versioned artifact in an open repository — and "
                       "the whole workflow runs on public data, without the "
                       "one-point-eight-terabyte cell-level atlas."):
            self.play(
                *[FadeOut(m) for m in self.mobjects],
                LaggedStart(
                    *[FadeIn(m) for pair in zip(chain, list(arrows) + [None])
                      for m in ([pair[0], pair[1]] if pair[1] else [pair[0]])],
                    lag_ratio=0.25),
                run_time=1.6)
            self.play(LaggedStart(*[FadeIn(c, shift=LEFT * 0.3) for c in chips],
                                  lag_ratio=0.2), run_time=1.0)

    # ================================================================== #
    # 7 — CLOSING: structure vs. predictability, side by side
    # ================================================================== #
    def scene_close(self):
        # left: recoverable structure (Complex I) · right: the prediction null
        op = self.make_operator_matrix(n=14, size=2.6).move_to([-3.3, 0.6, 0])
        opbox = Rectangle(width=2.6 * 4 / 14, height=2.6 * 4 / 14,
                          stroke_color=T.POSITIVE, stroke_width=2).move_to(
            op.get_corner(DOWN + RIGHT) + np.array([-0.37, 0.37, 0]))
        left_lab = VGroup(
            txt("recoverable structure", size=T.LABEL_SIZE, color=T.SIGNAL),
            txt("Complex I ✓", size=T.SMALL_SIZE, color=T.POSITIVE)
        ).arrange(DOWN, buff=0.15).next_to(op, DOWN, buff=0.4)

        zero = DashedLine(UP * 0.9, DOWN * 0.9, color=T.MUTED, stroke_width=2
                          ).move_to([3.3, 0.6, 0])
        rd = Dot(radius=0.1, color=T.SIGNAL).move_to(zero.get_center() + UP * 0.3 + LEFT * 0.06)
        sd = Dot(radius=0.1, color=T.MUTED).move_to(zero.get_center() + DOWN * 0.3 + LEFT * 0.02)
        right_lab = VGroup(
            txt("inductive prediction", size=T.LABEL_SIZE, color=T.MUTED),
            txt("unseen regulators ✗", size=T.SMALL_SIZE, color=T.NEGATIVE)
        ).arrange(DOWN, buff=0.15).next_to(zero, DOWN, buff=0.5)

        with self.beat("Recoverable structure is not the same as inductive "
                       "predictability."):
            self.crossfade_in(op, opbox, left_lab, zero, rd, sd, right_lab,
                              run_time=0.7)

        thesis = txt("Recoverable structure  ≠  inductive predictability",
                     size=T.SUB_SIZE, color=T.FG, weight="BOLD").to_edge(DOWN, buff=1.4)
        claude = txt("Claude helped us find the architecture — and the "
                     "boundary of the claim.", size=T.BODY_SIZE, color=T.CLAUDE)
        claude.next_to(thesis, DOWN, buff=0.35)
        with self.beat("Claude helped us find the architecture — and the "
                       "boundary of the claim."):
            self.play(FadeIn(thesis, shift=UP * 0.15), run_time=T.T_NORMAL)
            self.play(FadeIn(claude), run_time=T.T_NORMAL)
        self.wait(1.2)

        card = VGroup(
            txt("The Empirical Regulatory Operator of CD4⁺ T Cells",
                size=T.BODY_SIZE, color=T.FG, weight="BOLD"),
            txt("Built with Claude: Life Sciences", size=T.LABEL_SIZE, color=T.CLAUDE),
            txt("Researcher Track", size=T.SMALL_SIZE, color=T.MUTED),
            txt("github.com/cuentadesanti/cd4-perturbseq-regulator-atlas",
                size=T.SMALL_SIZE, color=T.MUTED),
        ).arrange(DOWN, buff=0.32)
        with self.beat("The empirical regulatory operator of CD4 T-cells — "
                       "built with Claude."):
            self.crossfade_in(card, run_time=0.8)
        self.wait(2.4)
