"""
This subpackage houses the cross-compiled Fortran extension module used for magnetic fieldline tracing.

Upon building and installing the ``mapflpy`` package (using the ``meson`` build system), the Fortran
source code is compiled into a Python extension module named ``mapflpy_fortran`` and placed within
this directory. This module provides the low-level functionality required for efficient magnetic
fieldline tracing in spherical coordinates.
"""