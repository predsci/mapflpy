import hashlib
from collections import namedtuple

import pooch
# python -m pytest -q /Users/rdavidson/MHDweb/mapflpy/tests --rootdir .

URLS = {
    "2143-mast2-cor/br002.h5": "https://predsci.com/~rdavidson/assets/2143-mast2-cor/br002.h5",
    "2143-mast2-cor/bt002.h5": "https://predsci.com/~rdavidson/assets/2143-mast2-cor/bt002.h5",
    "2143-mast2-cor/bp002.h5": "https://predsci.com/~rdavidson/assets/2143-mast2-cor/bp002.h5",
    "2143-mast2-hel/br002.h5": "https://predsci.com/~rdavidson/assets/2143-mast2-hel/br002.h5",
    "2143-mast2-hel/bt002.h5": "https://predsci.com/~rdavidson/assets/2143-mast2-hel/bt002.h5",
    "2143-mast2-hel/bp002.h5": "https://predsci.com/~rdavidson/assets/2143-mast2-hel/bp002.h5",
}

REGISTRY = {
    "2143-mast2-cor/br002.h5": "sha256:2af0563abc56933a91089284c584f613075c9cede63b91f68bf4767a0a5563d8",
    "2143-mast2-cor/bt002.h5": "sha256:09f0041c7d16057e0fba3ef9e5ea09db9cbc057ac17e712968356c35edb6f6ac",
    "2143-mast2-cor/bp002.h5": "sha256:f53e82c96395ad2f72c42a94d3195ea21d22167b7b22fcc1f19273f00bec9c18",
    "2143-mast2-hel/br002.h5": "sha256:2cad9d9dc55b0d6e213f6dde7c45e3c7686340096dda5e164d69406b2a4e0117",
    "2143-mast2-hel/bt002.h5": "sha256:7a52c64255adaf55df85d00d3cb772c19ad22dd23d98d264e067bc701b712e7d",
    "2143-mast2-hel/bp002.h5": "sha256:f4937469c7e737dd872dc8d1731d7fc2040245eb08be432ee11bdd2fd4ec420c",
}

FETCHER = pooch.create(
    path=pooch.os_cache("mapflpy"),
    base_url="https://predsci.com/~rdavidson/assets/",
    registry=REGISTRY,
    urls=URLS,
    env="MAPFLPY_CACHE",
)

MagneticFieldFiles = namedtuple("MagneticFieldFiles", ["br", "bt", "bp"])

def fetch_cor_magfiles() -> MagneticFieldFiles:
    """Download a file using pooch."""
    filepaths = [FETCHER.fetch(filename) for filename in URLS.keys() if 'cor' in filename]
    return MagneticFieldFiles(*filepaths)

def fetch_hel_magfiles() -> MagneticFieldFiles:
    """Download a file using pooch."""
    filepaths = [FETCHER.fetch(filename) for filename in URLS.keys() if 'hel' in filename]
    return MagneticFieldFiles(*filepaths)

d=1
