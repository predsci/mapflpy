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
except ModuleNotFoundError as e:
    try:
        from importlib.metadata import version as _pkg_version
        from importlib.metadata import PackageNotFoundError
        from pathlib import Path
        # Fallback to installed metadata (wheel/sdist)
        __version__ = _pkg_version("mapflpy")  # type: ignore[assignment]
    except PackageNotFoundError as e:  # dev/editable without metadata
        try:
            import tomllib  # Python 3.11+
        except ModuleNotFoundError:  # pragma: no cover
            import tomli as tomllib  # pip install tomli

        pyproject = Path(__file__).parents[1].resolve() / 'pyproject.toml'
        data = tomllib.loads(pyproject.read_text())

        project_version = data.get("project", {}).get("version", "0+unknown")
        project_version = project_version.replace('"', '').replace("'", '')
        __version__ = project_version
