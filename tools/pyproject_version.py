#!/usr/bin/env python3
from pathlib import Path


def _load_pyproject_toml(path: str | Path = "pyproject.toml") -> dict:
    try:
        import tomllib  # Python 3.11+
    except ModuleNotFoundError:  # pragma: no cover
        import tomli as tomllib  # pip install tomli

    pyproject = Path(path)
    if not pyproject.is_file():
        raise FileNotFoundError(f"Could not find pyproject.toml file: {path}")
    return tomllib.loads(pyproject.read_text())


def get_project_version(*args, **kwargs) -> str:
    """Get the project version from pyproject.toml.

    Returns
    -------
    str
        The project version string.
    """
    pyproject = _load_pyproject_toml(*args, **kwargs)
    project_version = (pyproject.get("project", {}).get("version", "0+unknown")
                       .replace('"', '')
                       .replace("'", ''))
    return project_version


if __name__ == "__main__":
    print(get_project_version())