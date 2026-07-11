"""ElevenLabs speech service that skips the voices-list call.

The stock ``ElevenLabsService`` always calls the ``voices()`` endpoint to
resolve a name to an id, which needs the ``voices_read`` API permission. A
text-to-speech-only key does not have it, so we construct the ``Voice``
straight from a known ``voice_id`` instead.
"""

import subprocess
from pathlib import Path

from elevenlabs import Voice
from manim_voiceover.services.base import SpeechService
from manim_voiceover.services.elevenlabs import ElevenLabsService
from manim_voiceover.helper import remove_bookmarks

# Legacy public voices (id -> character). Rachel = calm documentary narrator.
VOICE_IDS = {
    "Rachel": "21m00Tcm4TlvDq8ikWAM",
    "Adam": "pNInz6obpgDQGcFmaJgB",
    "Antoni": "ErXwobaYiN019PkySvjV",
    "Arnold": "VR6AewLTigWG4xSOukaG",
    "Sam": "yoZ06aMxZJJ28mfd3POQ",
}


class ElevenLabsByID(ElevenLabsService):
    def _select_voice(self, voice_name, voice_id):
        vid = voice_id or VOICE_IDS.get(voice_name) or VOICE_IDS["Rachel"]
        return Voice(voice_id=vid)


class SayService(SpeechService):
    """Zero-cost, offline speech via the macOS ``say`` command.

    No API key, no network, no transcription model. Bookmarks are not
    supported (they need word-level timestamps), but wrapping each chapter
    in a single ``self.voiceover(...)`` still auto-syncs the visuals to the
    narration length. Swap for ``ElevenLabsByID`` once a paid key is set.
    """

    def __init__(self, voice="Daniel", rate=170, ffmpeg="ffmpeg", **kwargs):
        self.say_voice = voice
        self.rate = int(rate)
        self.ffmpeg = ffmpeg
        SpeechService.__init__(self, **kwargs)  # no transcription_model

    def generate_from_text(self, text, cache_dir=None, path=None, **kwargs):
        cache = Path(cache_dir) if cache_dir is not None else Path(self.cache_dir)
        spoken = remove_bookmarks(text)
        input_data = {
            "input_text": spoken,
            "service": "say",
            "config": {"voice": self.say_voice, "rate": self.rate},
        }
        cached = self.get_cached_result(input_data, cache)
        if cached is not None:
            return cached

        audio_path = (self.get_audio_basename(input_data) + ".mp3"
                      if path is None else str(path))
        aiff = cache / (audio_path + ".aiff")
        subprocess.run(
            ["say", "-v", self.say_voice, "-r", str(self.rate),
             "-o", str(aiff), spoken],
            check=True,
        )
        subprocess.run(
            [self.ffmpeg, "-y", "-i", str(aiff), str(cache / audio_path)],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        aiff.unlink(missing_ok=True)
        return {"input_text": text, "input_data": input_data,
                "original_audio": audio_path}
