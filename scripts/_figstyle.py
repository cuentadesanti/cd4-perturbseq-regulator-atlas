"""Shared figure style — one entity, one colour, everywhere.

Import this into every analysis-figure script so a reviewer flipping through figs 26–31 is never
silently retrained (SAGA is always purple, interferon/ISG always red, generic/random always grey,
conditions always the grey→amber→red stimulation ramp).
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- canonical entity colours ---------------------------------------------------------------
SAGA = "#8e44ad"      # SAGA/chromatin (the "class purple")
ISG = "#c0392b"       # interferon / ISG / de-repression
GENERIC = "#9aa3b2"   # random / null / generic background / "other targets"
MEDIATOR = "#2f6fb0"
TCR = "#2e9e6b"
OTHER_CLASS = "#e0a441"
PROMOTED = "#d98c3f"
DEMOTED = "#7f8896"
MUTED = "#5d6779"      # reference lines / annotations
INK = "#222831"

CLASS_COLORS = {
    "SAGA/chromatin": SAGA, "Mediator": MEDIATOR, "TCR (context-specific)": TCR,
    "Other robust": OTHER_CLASS, "Repro-promoted": PROMOTED, "Demoted control": DEMOTED,
}
# grey (Rest) -> amber (Stim8hr) -> red (Stim48hr): cool→hot = more stimulation
CONDITION_RAMP = {"Rest": "#9aa3b2", "Stim8hr": "#e0a441", "Stim48hr": "#c0392b"}

# shared one-liners so jargon is defined on-figure exactly once per deck
ISG_DEF = "ISG = interferon-stimulated gene · de-repress = rises on knockdown"
# pointer used to harmonise the older 24×/19.2× figures with the PR-#8 specificity control
REQUALIFIED = ("magnitude is largely a general strong-perturbation effect — see the specificity "
               "control (fig 29); SAGA's distinction is the consistent de-repressive direction")


def apply_rc():
    plt.rcParams.update({
        "font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
        "axes.titlesize": 10.5, "axes.titleweight": "bold", "axes.edgecolor": "#c7ccd4",
        "axes.labelcolor": INK, "text.color": INK, "xtick.color": INK, "ytick.color": INK,
        "figure.facecolor": "white", "savefig.facecolor": "white",
    })


def callout(ax, text, xy, xytext, color=INK, fs=9, weight="bold", ha="center"):
    """Arrow annotation pointing at the decisive value."""
    ax.annotate(text, xy=xy, xytext=xytext, ha=ha, va="center", fontsize=fs, fontweight=weight,
                color=color, zorder=6,
                arrowprops=dict(arrowstyle="-|>", color=color, lw=1.4, shrinkA=2, shrinkB=3))


def footnote(fig, text, y=0.005):
    fig.text(0.005, y, text, fontsize=7.3, color=MUTED, ha="left", va="bottom", style="italic")
