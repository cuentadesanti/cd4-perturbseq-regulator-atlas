"""Five reusable components — the whole video is composed from these.

    txt / title      typographic helpers (consistent font + size)
    EvidenceCard     a labelled value chip (dataset facts, artifacts)
    MetricGate       a value approaching a pass/fail threshold on a track
    EvidenceChain    stacked evidence lines converging on a conclusion
    ClaimVerdict     a claim split into survives / weakens / fails
    spotlight        4-rectangle mask that darkens all but a target region

Keeping these small and generic is what buys the "one visual world" feel:
the same MetricGate renders the SAGA stability gate and the spectral edge.
"""

from manim import (
    VGroup,
    Group,
    Text,
    RoundedRectangle,
    Rectangle,
    Line,
    Dot,
    Arrow,
    DOWN,
    UP,
    LEFT,
    RIGHT,
    ORIGIN,
)
import numpy as np

import theme as T


# ---- typography ----------------------------------------------------------
def txt(s, size=T.BODY_SIZE, color=T.FG, weight="NORMAL", **kw):
    return Text(s, font=T.FONT, font_size=size, color=color, weight=weight, **kw)


def title(s, color=T.FG):
    return txt(s, size=T.TITLE_SIZE, color=color, weight="BOLD")


# ---- EvidenceCard --------------------------------------------------------
class EvidenceCard(VGroup):
    """A labelled value chip: a caption over a big value in an accent box."""

    def __init__(self, caption, value, accent=T.SIGNAL, width=3.2, height=1.5):
        super().__init__()
        box = RoundedRectangle(
            corner_radius=0.15,
            width=width,
            height=height,
            stroke_width=T.STROKE_MUTED,
            stroke_color=accent,
            fill_color=accent,
            fill_opacity=0.07,
        )
        cap = txt(caption, size=T.SMALL_SIZE, color=T.MUTED)
        val = txt(value, size=T.SUB_SIZE, color=T.FG, weight="BOLD")
        # shrink the value (and caption) to fit inside the card with padding
        max_w = width - 0.5
        if val.width > max_w:
            val.scale_to_fit_width(max_w)
        if cap.width > max_w:
            cap.scale_to_fit_width(max_w)
        stack = VGroup(cap, val).arrange(DOWN, buff=0.12).move_to(box)
        self.add(box, stack)
        self.box, self.caption_mob, self.value_mob = box, cap, val


# ---- MetricGate ----------------------------------------------------------
class MetricGate(VGroup):
    """A horizontal track with a pass/fail gate and a value marker.

    Build it at value=lo, then animate the marker to the true value with
    ``gate.marker_to(value)`` as an .animate target. If the value lands
    short of the gate it recolors to WARNING, communicating "recovered but
    below threshold" without a word.
    """

    def __init__(self, label, value, gate, lo=0.0, hi=1.0, width=7.0,
                 accent=T.SIGNAL):
        super().__init__()
        self.lo, self.hi, self.width_ = lo, hi, width
        self.gate_val, self.target_val = gate, value

        track = Line(LEFT * width / 2, RIGHT * width / 2,
                     stroke_width=T.STROKE_MUTED, color=T.MUTED)

        gx = self._x(gate)
        gate_line = Line(UP * 0.35, DOWN * 0.35, stroke_width=T.STROKE_MAIN,
                         color=T.WARNING).move_to(track.get_start() + RIGHT * (gx + width / 2))
        gate_lab = txt(f"gate {gate:.2f}", size=T.SMALL_SIZE, color=T.WARNING)
        gate_lab.next_to(gate_line, UP, buff=0.18)

        passed = value >= gate
        mark_color = accent if passed else T.WARNING
        marker = Dot(radius=0.11, color=mark_color)
        marker.move_to(track.get_start() + RIGHT * (self._x(value) + width / 2))
        val_lab = txt(f"{value:.2f}", size=T.LABEL_SIZE, color=mark_color, weight="BOLD")
        val_lab.next_to(marker, DOWN, buff=0.2)

        name = txt(label, size=T.LABEL_SIZE, color=T.FG)
        name.next_to(track, UP, buff=0.9).align_to(track, LEFT)

        self.add(name, track, gate_line, gate_lab, marker, val_lab)
        self.track = track
        self.marker = marker
        self.val_lab = val_lab
        self.gate_line = gate_line

    def _x(self, v):
        frac = (v - self.lo) / (self.hi - self.lo)
        return (frac - 0.5) * self.width_

    def marker_position(self, v):
        return self.track.get_start() + RIGHT * (self._x(v) + self.width_ / 2)


# ---- EvidenceChain -------------------------------------------------------
class EvidenceChain(VGroup):
    """Independent evidence lines stacked with '+', converging (↓) on a claim."""

    def __init__(self, items, conclusion, accent=T.SIGNAL):
        super().__init__()
        rows = VGroup()
        for i, it in enumerate(items):
            rows.add(txt(it, size=T.LABEL_SIZE, color=T.FG))
        rows.arrange(DOWN, buff=0.3)
        plus = VGroup(*[txt("+", size=T.LABEL_SIZE, color=T.MUTED)
                        for _ in range(len(items) - 1)])
        # interleave + between rows
        for j, p in enumerate(plus):
            p.move_to((rows[j].get_bottom() + rows[j + 1].get_top()) / 2)
        arrow = Arrow(UP * 0.25, DOWN * 0.25, buff=0, color=accent,
                      stroke_width=T.STROKE_MAIN)
        concl = txt(conclusion, size=T.BODY_SIZE, color=accent, weight="BOLD")
        block = VGroup(VGroup(rows, plus), arrow, concl).arrange(DOWN, buff=0.35)
        self.add(block)
        self.rows, self.arrow, self.conclusion = rows, arrow, concl


# ---- ClaimVerdict --------------------------------------------------------
class ClaimVerdict(VGroup):
    """A statement with a colored verdict token (✓ survives / ✗ fails)."""

    def __init__(self, claim, verdict, ok=True):
        super().__init__()
        mark = "✓" if ok else "✗"
        color = T.POSITIVE if ok else T.NEGATIVE
        c = txt(claim, size=T.LABEL_SIZE, color=T.FG)
        v = txt(f"{mark} {verdict}", size=T.LABEL_SIZE, color=color, weight="BOLD")
        VGroup(c, v).arrange(RIGHT, buff=0.5)
        self.add(c, v)
        self.claim_mob, self.verdict_mob = c, v


# ---- spotlight -----------------------------------------------------------
def spotlight(target, opacity=0.82, span=60.0):
    """Four dark rectangles that cover everything except ``target``'s bbox.

    Big ``span`` so it still darkens the scene under camera zoom. Returns a
    Group; FadeIn it to mask, FadeOut to release.
    """
    c = target.get_center()
    w, h = target.width, target.height
    x0, x1 = c[0] - w / 2, c[0] + w / 2
    y0, y1 = c[1] - h / 2, c[1] + h / 2

    def panel(cx, cy, pw, ph):
        r = Rectangle(width=pw, height=ph, stroke_width=0,
                      fill_color=T.BG, fill_opacity=opacity)
        r.move_to([cx, cy, 0])
        return r

    top = panel(c[0], y1 + span / 2, span * 2, span)
    bot = panel(c[0], y0 - span / 2, span * 2, span)
    left = panel(x0 - span / 2, c[1], span, h)
    right = panel(x1 + span / 2, c[1], span, h)
    return Group(top, bot, left, right)
