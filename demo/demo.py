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
    Line,
    DashedLine,
    Arrow,
    Dot,
    FadeIn,
    FadeOut,
    Write,
    Create,
    GrowFromCenter,
    Transform,
    ReplacementTransform,
    TransformFromCopy,
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

        self.narrate = os.environ.get("NARRATE", "").lower()
        if self.narrate == "say":
            from voice import SayService
            self.set_speech_service(
                SayService(voice=os.environ.get("SAY_VOICE", "Samantha"),
                           rate=int(os.environ.get("SAY_RATE", "138")),
                           ffmpeg=os.environ.get("FFMPEG", "ffmpeg")))
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
        if self.narrate in ("say", "eleven"):
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

        # one regulator, one gene
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
        with self.beat("For each inhibited regulator, the screen measures how "
                       "thousands of downstream genes respond across many "
                       "cells."):
            self.play(FadeIn(regA), FadeIn(crispr), FadeIn(cross),
                      run_time=T.T_NORMAL)
            self.play(Create(arr0), FadeIn(geneB), run_time=T.T_MICRO)

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
        zval = txt("z = +2.4", size=T.SUB_SIZE, color=T.POSITIVE,
                   weight="BOLD").move_to([1.5, -0.2, 0])
        znote = txt("signed, standardized effect", size=T.SMALL_SIZE,
                    color=T.MUTED).next_to(zval, DOWN, buff=0.25)
        with self.beat("Those observations become one signed, standardized "
                       "effect."):
            self.play(Create(axis), FadeIn(axlab), FadeIn(ctrl), FadeIn(ctrl_lab),
                      FadeIn(pert), FadeIn(pert_lab), run_time=T.T_NORMAL)
            self.play(pert.animate.shift(RIGHT * 1.4), run_time=0.7,
                      rate_func=T.EASE_MOVE)
            self.play(FadeOut(axis), FadeOut(axlab), FadeOut(ctrl_lab),
                      FadeOut(pert_lab),
                      FadeTransform(VGroup(ctrl, pert), zval), run_time=T.T_NORMAL)
            self.play(FadeIn(znote), run_time=T.T_MICRO)

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
            self.play(FadeOut(regA), FadeOut(cross), FadeOut(crispr),
                      FadeOut(arr0), FadeOut(geneB), FadeOut(znote), run_time=0.4)
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
        self.wait(1.8)  # let the main discovery land

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
        with self.beat("And when we asked whether side-features could predict "
                       "completely unseen regulators, real and shuffled "
                       "features scored identically. We stopped, and reported "
                       "the null."):
            self.play(FadeOut(gate), FadeOut(verdicts), run_time=0.5)
            self.play(FadeIn(q), Create(zero), FadeIn(zlab),
                      FadeIn(real_lab), FadeIn(shuf_lab), run_time=T.T_NORMAL)
            self.play(real_dot.animate.move_to([x0 - 0.08, ay, 0]),
                      shuf_dot.animate.move_to([x0 - 0.02, by_, 0]),
                      run_time=T.T_REVEAL, rate_func=T.EASE_MOVE)
            self.add(real_dot, shuf_dot)
            self.play(FadeIn(concl, shift=UP * 0.15), run_time=T.T_NORMAL)

        with self.beat("Claude helped generate the analysis — but its most "
                       "valuable work was deciding which conclusions had to "
                       "weaken, or fail."):
            self.play(Indicate(head, color=T.CLAUDE, scale_factor=1.05),
                      run_time=T.T_NORMAL)
            self.wait(0.6)

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
