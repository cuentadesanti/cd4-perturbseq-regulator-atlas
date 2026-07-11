"""Design system for the operator demo.

One place for every color, size, font and easing choice so the whole video
reads as a single visual world instead of a stack of slides. Colors are
*semantic*, not decorative:

    SIGNAL   (blue)   supported structure / the operator
    CLAUDE   (orange) the process, Claude-in-the-loop
    POSITIVE (green)  a claim that survived
    WARNING  (yellow) limited / partial evidence
    NEGATIVE (red)    a claim that was rejected / a null
    MUTED    (gray)   noise, baselines, shuffled controls
"""

from manim import rate_functions

# ---- palette -------------------------------------------------------------
BG = "#0E1117"
FG = "#F3F4F6"
MUTED = "#8B95A5"
SIGNAL = "#58C4DD"
CLAUDE = "#D99A6C"
POSITIVE = "#7CC576"
WARNING = "#F2C14E"
NEGATIVE = "#E06C75"
FAINT = "#252B36"  # panel fills, grid lines

# ---- type ----------------------------------------------------------------
# Helvetica Neue ships with macOS; swap here if you install Inter.
FONT = "Helvetica Neue"
TITLE_SIZE = 52
SUB_SIZE = 34
BODY_SIZE = 30
LABEL_SIZE = 24
SMALL_SIZE = 20

# ---- strokes -------------------------------------------------------------
STROKE_MAIN = 3.0
STROKE_MUTED = 1.5

# ---- easing grammar ------------------------------------------------------
# appearance -> ease out, movement -> smooth, big reveal -> ease in-out.
EASE_APPEAR = rate_functions.ease_out_cubic
EASE_MOVE = rate_functions.smooth
EASE_REVEAL = rate_functions.ease_in_out_sine
EASE_SNAP = rate_functions.rush_from  # quick, dry decisions

# ---- timing (seconds) ----------------------------------------------------
T_MICRO = 0.4
T_NORMAL = 1.0
T_REVEAL = 1.8
T_HOLD = 1.0  # the stillness after a climax is part of the animation
