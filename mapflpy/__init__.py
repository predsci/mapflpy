from . import tracer, scripts, utils, _typing

from importlib.metadata import version, PackageNotFoundError


try:
    __version__ = version("mapflpy")
except PackageNotFoundError:
    # Dev fallback: not installed (e.g., running from a checkout)
    __version__ = "0.0.0"
