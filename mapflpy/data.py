"""
Module for fetching HDF5 assets used through examples.

This module uses the `pooch` library to manage the downloading and caching of
HDF5 files containing magnetic field data. It defines functions to fetch
coronal and heliospheric magnetic field files, returning their paths in a
named tuple for easy access.

Currently, the available Coronal and Heliospheric magnetic field files are hosted
at https://predsci.com/~rdavidson/assets/ *viz.* a Thermo 2 steady-state run for
Carrington Rotation 2143.
"""


from __future__ import annotations
import hashlib
from collections import namedtuple
from pathlib import Path

import pooch


BASE_URL = "https://predsci.com/~rdavidson/assets/"
"""Base URL hosting magnetic field file assets."""

REGISTRY = Path(__file__).parents[1] / "registry.txt"
"""Path to the registry file listing available magnetic field files and their checksums."""


FETCHER = pooch.create(
    path=pooch.os_cache("mapflpy"),
    base_url=BASE_URL,
    env="MAPFLPY_CACHE",
)
"""Pooch fetcher for downloading and caching magnetic field files.

.. note::
    The cache directory can be overridden by setting the ``MAPFLPY_CACHE``
    environment variable to a desired path. Otherwise, the default cache
    directory is platform-dependent, as determined by :func:`pooch.os_cache`.
"""

FETCHER.load_registry(REGISTRY)


MagneticFieldFiles = namedtuple("MagneticFieldFiles", ["br", "bt", "bp"])
MagneticFieldFiles.__doc__ = (
    """Container for magnetic field file paths.

    Attributes
    ----------
    br : str
        File path to radial magnetic field data.
    bt : str
        File path to theta magnetic field data.
    bp : str
        File path to phi magnetic field data.
    """
)


def fetch_cor_magfiles() -> MagneticFieldFiles:
    """Download sample coronal magnetic field files using pooch.

    Returns
    -------
    MagneticFieldFiles
        Named tuple containing file paths for br, bt, and bp magnetic field data.
    """
    filepaths = [FETCHER.fetch(filename) for filename in FETCHER.registry.keys() if 'cor' in filename]
    return MagneticFieldFiles(*filepaths)


def fetch_hel_magfiles() -> MagneticFieldFiles:
    """Download sample heliospheric magnetic field files using pooch.

    Returns
    -------
    MagneticFieldFiles
        Named tuple containing file paths for br, bt, and bp magnetic field data.
    """
    filepaths = [FETCHER.fetch(filename) for filename in FETCHER.registry.keys() if 'hel' in filename]
    return MagneticFieldFiles(*filepaths)


def magfiles() -> list[str]:
    """List all available magnetic field files in the registry.

    Returns
    -------
    list[str]
        List of available magnetic field file names.
    """
    return list(FETCHER.registry.keys())