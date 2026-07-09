import sys
from pathlib import Path

# Drivers run from scripts/ (bare imports, like the repo's `import _figstyle`);
# put scripts/ on sys.path so the tests can `import _opkernels` the same way.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
