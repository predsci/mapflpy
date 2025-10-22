"""
A Python package for tracing magnetic fieldlines in spherical coordinates.

This package provides tools for tracing magnetic fieldlines in spherical coordinate systems,
using PSI's cross-compiled ``mapfl`` Fortran library. The following modules are intended to
allow users a high-level interface to the underlying Fortran routines, as well as utilities for
visualizing and analyzing the traced fieldlines.
"""


from __future__ import annotations

import os
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
