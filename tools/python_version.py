from __future__ import annotations
import re
from pathlib import Path

def get_python_version(path: str | Path = ".python-version") -> str | None:
    """
    Read .python-version and return 'MAJOR.MINOR' (e.g. '3.13').
    Accepts values like '3.13', '3.13.1', or pyenv-style strings.
    """
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"Could not find .python-version file: {path}")
    txt = p.read_text(encoding="utf-8").strip().splitlines()[0]
    m = re.search(r"(\d+)\.(\d+)", txt)  # grab major.minor; ignore patch/suffix
    if not m:
        raise ValueError(f"Could not parse Python version from {path}")
    return f"{m.group(1)}.{m.group(2)}" if m else None


if __name__ == "__main__":
    print(get_python_version())
