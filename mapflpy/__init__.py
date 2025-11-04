"""
A Python package for tracing magnetic fieldlines in spherical coordinates.

This package provides tools for tracing magnetic fieldlines in spherical coordinate systems,
using PSI's cross-compiled ``mapfl`` Fortran library. The following modules are intended to
allow users a high-level interface to the underlying Fortran routines, as well as utilities for
visualizing and analyzing the traced fieldlines.
"""

# mapflpy/__init__.py
try:
    # If Meson generated this file:
    from ._version import __version__  # type: ignore[attr-defined]
except Exception:
    try:
        # Fallback to installed metadata (wheel/sdist)
        from importlib.metadata import version as _pkg_version
        __version__ = _pkg_version("mapflpy")  # type: ignore[assignment]
    except Exception:  # dev/editable without metadata
        __version__ = "0+unknown"  # type: ignore[assignment]

