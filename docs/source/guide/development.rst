.. _development:

Development
===========

This document explains how to set up a dev environment, run tests, and build docs.

Requirements
------------
- Python 3.10+ (recommended)
- Conda (miniforge/miniconda) or pyenv + venv
- Git, C/C++ compilers; Fortran toolchain (if building Fortran parts)
- (macOS) Xcode CLT

Quick Start (Conda)
-------------------
.. code-block:: bash

   git clone git@github.org:predsci/mapflpy.git
   cd mapflpy
   conda env create -f envs/dev.yml   # or: conda create -n mapflpy-dev python=3.11
   conda activate mapflpy-dev
   pip install -e .[dev]              # installs nox, lint/test deps, etc.
   pre-commit install                 # enable git hooks

Useful Nox Sessions
-------------------
.. code-block:: bash

   nox -s lint          # ruff/black/mypy
   nox -s tests         # pytest (CPU)
   nox -s docs          # Sphinx + gallery (serial, stable)
   nox -s build         # build wheel
   nox -s repair        # delocate/auditwheel as applicable

Documentation
-------------
Docs live in ``docs/`` and use Sphinx + sphinx-gallery.

- Build locally: ``nox -s docs``
- Gallery execution is serial to avoid native teardown crashes on macOS.

Coding Style & Type Hints
-------------------------
- Ruff + Black for style/formatting; Mypy for typing (strict-ish).
- NumPy/Google-style docstrings; public APIs documented.

Native/Fortran Notes
--------------------
If you touch Fortran extensions:
- Prefer Meson/Ninja for builds; keep Python wrappers thin.
- Memory layout: be explicit about Fortran-contiguous arrays.
- Add regression tests for array shape/ownership semantics.

