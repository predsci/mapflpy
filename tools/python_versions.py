#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path
from typing import Any


def _load_pyproject_toml(path: str | Path = "pyproject.toml") -> dict:
    try:
        import tomllib  # Python 3.11+
    except ModuleNotFoundError:  # pragma: no cover
        import tomli as tomllib  # pip install tomli

    pyproject = Path(path)
    if not pyproject.is_file():
        raise FileNotFoundError(f"Could not find pyproject.toml file: {path}")
    return tomllib.loads(pyproject.read_text())


def python_versions(
    pyproject: dict[str, Any], *, max_version: str | None = None
) -> list[str]:
    """
    Read a list of supported Python versions. Without ``max_version``, this
    will read the trove classifiers (recommended). With a ``max_version``, it
    will read the requires-python setting for a lower bound, and will use the
    value of ``max_version`` as the upper bound. (Reminder: you should never
    set an upper bound in ``requires-python``).

    .. note::
        This function requires the ``packaging`` library if ``max_version``
        is provided. This optional dependency can be avoided by specifying
        Python versions via classifiers (in the pyproject.toml) instead.

    .. note::
        This function is directly adapted from the same function in Nox
        (https://nox.thea.codes/en/stable/), licensed under the Apache License
        Version 2.0. If you are using Nox in your project, consider using
        Nox's built-in project module instead.
    """
    if max_version is None:
        # Classifiers are a list of every Python version
        from_classifiers = [
            c.split()[-1]
            for c in pyproject.get("project", {}).get("classifiers", [])
            if c.startswith("Programming Language :: Python :: 3.")
        ]
        if from_classifiers:
            return from_classifiers
        msg = 'No Python version classifiers found in "project.classifiers"'
        raise ValueError(msg)

    requires_python_str = pyproject.get("project", {}).get("requires-python", "")
    if not requires_python_str:
        msg = 'No "project.requires-python" value set'
        raise ValueError(msg)

    try:
        from packaging.specifiers import SpecifierSet
    except ModuleNotFoundError:
        raise ImportError(
            'The "packaging" library is required to parse "requires-python". '
            "This optional dependency can be avoided by specifying Python versions "
            "via classifiers (in the pyproject.toml) instead."
        ) from None
    for spec in SpecifierSet(requires_python_str):
        if spec.operator in {">", ">=", "~="}:
            min_minor_version = int(spec.version.split(".")[1])
            break
    else:
        msg = 'No minimum version found in "project.requires-python"'
        raise ValueError(msg)

    max_minor_version = int(max_version.split(".")[1])

    return [f"3.{v}" for v in range(min_minor_version, max_minor_version + 1)]


def get_python_versions(*args, **kwargs) -> list[str]:
    """Get the supported Python versions from pyproject.toml."""
    pyproject = _load_pyproject_toml(*args, **kwargs)
    return python_versions(pyproject)


if __name__ == "__main__":
    print(get_python_versions())
