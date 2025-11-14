#!/usr/bin/env python3
import os
import textwrap
from pathlib import Path


def init_version():
    # Load TOML (tomllib on 3.11+, fallback to tomli)
    try:
        import tomllib  # Python 3.11+
    except ModuleNotFoundError:  # pragma: no cover
        import tomli as tomllib  # pip install tomli

    pyproject = Path(__file__).parents[1].resolve() / 'pyproject.toml'
    with open(pyproject, 'rb') as f:
        data = tomllib.load(f)

    version = data.get("project", {}).get("version", '0.0.0')
    return version.replace('"', '').replace("'", '')

if __name__ == "__main__":
    version = init_version()
    print(version)