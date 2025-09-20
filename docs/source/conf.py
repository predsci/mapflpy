import sys
from datetime import datetime
from pathlib import Path

sys.path.insert(0, Path(__file__).resolve().parents[2].as_posix())
import mapflpy

project = "mapflpy"
author = "Predictive Science Inc"
copyright = f"{datetime.now():%Y}, {author}"
version = mapflpy.__version__
release = mapflpy.__version__
language = "en"

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx_autodoc_typehints",
    "sphinx_gallery.gen_gallery",
]

# --- Theme ---
logo = "https://www.predsci.com/corona/apr2024eclipse/images/psi_logo.png"
numfig = False
html_theme = "sphinx_book_theme"
html_static_path = ["_static"]
html_favicon = logo
html_theme_options = {
    "show_prev_next": False,
    "navigation_with_keys": False,
    "show_nav_level": 2,
    "navigation_depth": 3,
    "logo": {
        "text": f"{project} v{version}",
        "image_light": logo,
        "image_dark": logo,
    },
    'icon_links': [
        {
            'name': 'PSI Home',
            'url': 'https://www.predsci.com/',
            'icon': 'fa fa-home fa-fw',
            "type": "fontawesome",
        },
        {
            'name': 'Github',
            'url': 'https://github.com/predsci/',
            "icon": "fa-brands fa-github fa-fw",
            "type": "fontawesome",
        },
        {
            'name': 'Bitbucket',
            'url': 'https://bitbucket.org/predsci/',
            "icon": "fa-brands fa-bitbucket fa-fw",
            "type": "fontawesome",
        },
        {
            'name': 'Contact',
            'url': 'https://www.predsci.com/portal/contact.php',
            'icon': 'fa fa-envelope fa-fw',
            "type": "fontawesome",
        },
    ],
}

# --- MyST Markdown ---
myst_enable_extensions = [
    "dollarmath",   # $math$
    "amsmath",      # equations
    "colon_fence",  # ::: fences
]

# --- Autodoc / API style ---
autosummary_generate = True
autodoc_typehints = "description"   # show type hints in the doc body
autoclass_content = "both"          # include __init__ docstrings
add_module_names = False            # shorten paths in headers

# --- Intersphinx cross-links ---
# intersphinx_mapping = {
#     "python": (
#         "https://docs.python.org/3/",
#         None,
#     ),
#     "numpy": (
#         "https://numpy.org/doc/stable/",
#         None,
#     ),
#     "scipy": (
#         "https://docs.scipy.org/doc/scipy/reference/",
#         None,
#     ),
# }

# --- Sphinx-Gallery (Examples) ---
from sphinx_gallery.sorting import FileNameSortKey

sphinx_gallery_conf = {
    "examples_dirs": ["examples"],
    "gallery_dirs": ["../generated"],
    "within_subsection_order": FileNameSortKey,
    "filename_pattern": r".*\.py",
    "download_all_examples": False,
    "remove_config_comments": True,
    # Run examples during build (True). For quick local builds, set False.
    "plot_gallery": True,
}
