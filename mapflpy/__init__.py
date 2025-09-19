from . import tracer, scripts, utils, _typing

from importlib.metadata import version, PackageNotFoundError

__version__ = "0.0.0"

try:
    __version__ = version("mapflpy")
except PackageNotFoundError:
    from pathlib import Path
    try:
        import tomllib
    except ModuleNotFoundError:
        import tomli as tomllib

    try:
        data = tomllib.loads(
            (Path(__file__).resolve().parents[1] / "pyproject.toml").read_text(encoding="utf-8")
        )
        __version__ = data["project"]["version"]
    except Exception:
        pass

