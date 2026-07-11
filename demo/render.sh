#!/usr/bin/env bash
# Render the operator demo.
#
#   ./render.sh preview          fast 480p, silent
#   ./render.sh final            1080p30, silent
#   ./render.sh final say        1080p30 + offline macOS narration (no key)
#   ./render.sh final eleven     1080p30 + ElevenLabs narration (needs paid key)
#
# ElevenLabs: put the key in eleven_labs_raw_key.txt (or set ELEVEN_API_KEY),
# and a personal/allowed voice id in ELEVEN_VOICE_ID. The current free key
# returns HTTP 402 for library voices, so `say` is the working default.
set -euo pipefail
cd "$(dirname "$0")"

MODE="${1:-preview}"      # preview | final
NARRATE="${2:-}"          # ''(silent) | say | eleven

# find a .venv/bin/python by walking up from here; fall back to python3
find_python() {
  local d="$PWD"
  while [[ "$d" != "/" ]]; do
    [[ -x "$d/.venv/bin/python" ]] && { echo "$d/.venv/bin/python"; return; }
    d="$(dirname "$d")"
  done
  command -v python3
}
PY="${PYTHON:-$(find_python)}"
export PATH="/opt/homebrew/bin:$PATH"       # ffmpeg
export FFMPEG="$(command -v ffmpeg)"

if [[ "$NARRATE" == "eleven" ]]; then
  if [[ -z "${ELEVEN_API_KEY:-}" && -f eleven_labs_raw_key.txt ]]; then
    ELEVEN_API_KEY="$(tr -d '[:space:]' < eleven_labs_raw_key.txt)"
  fi
  export ELEVEN_API_KEY ELEVENLABS_API_KEY="${ELEVEN_API_KEY:-}"
fi
if [[ "$NARRATE" == "deepgram" ]]; then
  if [[ -z "${DEEPGRAM_API_KEY:-}" && -f deepgram_api_key.txt ]]; then
    DEEPGRAM_API_KEY="$(tr -d '[:space:]' < deepgram_api_key.txt)"
  fi
  export DEEPGRAM_API_KEY
fi
export NARRATE="$NARRATE"

case "$MODE" in
  preview) FLAGS="-ql" ;;
  final)   FLAGS="-qh --fps 30" ;;
  *) echo "unknown mode: $MODE" >&2; exit 1 ;;
esac

"$PY" -m manim $FLAGS demo.py FullDemo

# Optional: burn the manim-voiceover subtitles into a *_subbed.mp4 for silent
# async viewing.  BURN_SUBS=1 ./render.sh final say
MP4="$(ls -t media/videos/demo/*/FullDemo.mp4 2>/dev/null | head -1)"
SRT="${MP4%.mp4}.srt"
if [[ "${BURN_SUBS:-}" == "1" && -n "$MP4" && -f "$SRT" ]]; then
  SUBBED="${MP4%.mp4}_subbed.mp4"
  ffmpeg -y -i "$MP4" -vf "subtitles=${SRT}:force_style='FontName=Helvetica Neue,FontSize=16,Outline=1,Shadow=0,MarginV=36'" -c:a copy "$SUBBED"
  echo "Burned subtitles → $SUBBED"
fi
echo "Output under media/videos/demo/"
