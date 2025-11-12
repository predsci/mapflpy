"""
Module for fetching HDF5 assets used through examples.

This module uses the ``pooch`` library to manage the downloading and caching of
HDF5 files containing magnetic field data. It defines functions to fetch
coronal and heliospheric magnetic field files, returning their paths in a
named tuple for easy access.

Currently, the available Coronal and Heliospheric magnetic field files are hosted
at https://www.predsci.com/doc/assets/ *viz.* a Thermo 2 steady-state run for
Carrington Rotation 2143.
"""


from __future__ import annotations
from collections import namedtuple

try:
    import pooch
except ImportError as e:
    raise ImportError(
        "Missing the optional 'pooch' dependency required for data fetching. "
        "Please install it via pip or conda to access the necessary datasets."
    ) from e


REGISTRY = {
	"2143-mast2-cor/br002.h5": "sha256:2af0563abc56933a91089284c584f613075c9cede63b91f68bf4767a0a5563d8",
	"2143-mast2-cor/bt002.h5": "sha256:09f0041c7d16057e0fba3ef9e5ea09db9cbc057ac17e712968356c35edb6f6ac",
	"2143-mast2-cor/bp002.h5": "sha256:f53e82c96395ad2f72c42a94d3195ea21d22167b7b22fcc1f19273f00bec9c18",
	"2143-mast2-hel/br002.h5": "sha256:2cad9d9dc55b0d6e213f6dde7c45e3c7686340096dda5e164d69406b2a4e0117",
	"2143-mast2-hel/bt002.h5": "sha256:7a52c64255adaf55df85d00d3cb772c19ad22dd23d98d264e067bc701b712e7d",
	"2143-mast2-hel/bp002.h5": "sha256:f4937469c7e737dd872dc8d1731d7fc2040245eb08be432ee11bdd2fd4ec420c",
}
"""Registry of available magnetic field files with their SHA256 hashes. 

This registry is used by the pooch fetcher to verify the integrity of
downloaded files, and is primarily intended for building sphinx-gallery
examples that require MHD data files.

The files listed here correspond to a Thermo 2 steady-state run for
Carrington Rotation 2143 – both the coronal and heliospheric domains.
"""


BASE_URL = "https://www.predsci.com/doc/assets/"
"""Base URL hosting magnetic field file assets.
"""


FETCHER = pooch.create(
    path=pooch.os_cache("psi"),
    base_url=BASE_URL,
    registry=REGISTRY,
    env="MAPFLPY_CACHE",
)
"""Pooch fetcher for downloading and caching magnetic field files.

.. note::
    The cache directory can be overridden by setting the ``MAPFLPY_CACHE``
    environment variable to a desired path. Otherwise, the default cache
    directory is platform-dependent, as determined by :func:`pooch.os_cache`.
    
.. note::
    The default (os-dependent) cache directory stores assets under a
    subdirectory named ``psi``. The reason for this naming choice – as opposed
    to ``mapflpy`` – is to maintain consistency with other PredSci packages
    that utilize the same asset hosting and caching mechanism.
"""


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