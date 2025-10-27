from __future__ import annotations

import builtins
import os
import sys
import importlib
import pkgutil
import re
import ast
import inspect
from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from types import ModuleType
from typing import Iterable

from sphinx_gallery.sorting import ExplicitOrder

sys.path.insert(0, Path(__file__).resolve().parents[2].as_posix())
import mapflpy


# ------------------------------------------------------------------------------
# Helper Functions
# ------------------------------------------------------------------------------

def _check_exclusions(name, patterns):
    for pattern in patterns:
        if re.search(pattern, name):
            return True
    return False

def _get_base_exclusions(private: bool, tests: bool, dunder: bool) -> list[str]:
    patterns = []
    if private:
        patterns.append("(?:^|\\.)_[^.]+(?:\\.|$)")
    if tests:
        patterns.append("(?:^|\\.)(?:tests?|test_[^.]+|[^.]+_test)(?:\\.|$)")
    if dunder:
        patterns.append("(?:^|\\.)__[^.]+__(?:\\.|$)")
    return patterns


@dataclass(slots=True)
class Node:
    name: str
    kind: str
    children: list[Node] = field(default_factory=list)


def build_node_tree(root: str | ModuleType, *args, **kwargs) -> Node:
    """
    Build a JSON-ish API tree suitable for Jinja templating.

    Parameters
    ----------
    root
        A package or module object, or its importable dotted name.
    sort_children
        If True, children lists are sorted deterministically (kind, name).

    Returns
    -------
    dict
        {"name": ..., "kind": ..., "children": [...]}
    """
    mod = importlib.import_module(root) if isinstance(root, str) else root
    node = _module_to_node(mod, *args, **kwargs)
    return node


# -------- Internal helpers ----------------------------------------------------

_KIND_ORDER = {
    "package": 0,
    "module": 1,
    "class": 2,
    "exception": 3,
    "method": 4,
    "class_method": 5,
    "static_method": 6,
    "property": 7,
    "data": 8,
    "function": 9,
    "attribute": 10,
    "descriptor": 11,
}


def _is_package(module: ModuleType) -> bool:
    return hasattr(module, "__path__")


def _fqn_for_class_member(owner: type, name: str, obj: object | None = None) -> str:
    """
    Build a fully-qualified dotted name for an attribute on `owner`.
    Prefer the owner's module/qualname; fall back to the object's module.
    """
    mod = getattr(owner, "__module__", None) or getattr(obj, "__module__", None) or "builtins"
    qual = getattr(owner, "__qualname__", getattr(owner, "__name__", ""))
    return f"{mod}.{qual}.{name}"


def _fqn_for_module_member(module: ModuleType, name: str, value: object) -> str:
    """Fully qualified dotted name for a member of `module`."""
    if inspect.isclass(value) or inspect.isfunction(value):
        mod = getattr(value, "__module__", module.__name__) or module.__name__
        qual = getattr(value, "__qualname__", getattr(value, "__name__", name))
        return f"{mod}.{qual}"
    # everything else: data/consts/enums instances, etc.
    return f"{module.__name__}.{name}"



def _iter_immediate_children(module: ModuleType) -> Iterable[tuple[str, bool]]:
    """
    Yield (child_name, is_pkg) for *immediate* children only.

    Uses pkgutil.iter_modules over module.__path__ if package; otherwise yields none.
    """
    if not _is_package(module):
        return []
    for mi in pkgutil.iter_modules(module.__path__):  # immediate only
        yield mi.name, mi.ispkg


def _defined_and_imported_names_via_ast(module: ModuleType) -> tuple[set[str], set[str]]:
    """
    Return (defined_names, imported_names) at module top level using AST.

    - defined_names: names assigned/annotated/aug-assigned, as well as class/def names
    - imported_names: direct names bound by 'import' and 'from ... import ...'

    If source is unavailable (C extensions, builtins), both sets may be empty.
    """
    src = None
    defined: set[str] = set()
    imported: set[str] = set()

    try:
        src = inspect.getsource(module)
    except:
        filename = getattr(module, "__file__", None)
        if filename and filename.endswith(".py"):
            try:
                with open(filename, "r", encoding="utf-8") as f:
                    src = f.read()
            except:
                pass

    try:
        tree = ast.parse(src)
    except:
        return defined, imported

    for node in tree.body:  # only top-level
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            defined.add(node.name)
        elif isinstance(node, ast.Assign):
            for tgt in node.targets:
                defined |= _collect_assigned_names(tgt)
        elif isinstance(node, (ast.AnnAssign, ast.AugAssign)):
            defined |= _collect_assigned_names(node.target)
        elif isinstance(node, ast.Import):
            for alias in node.names:
                imported.add(alias.asname or alias.name.split(".")[0])
        elif isinstance(node, ast.ImportFrom):
            # from .x import a as b
            for alias in node.names:
                if alias.name == "*":
                    # Can't know; skip (treat as neither defined nor imported here)
                    continue
                imported.add(alias.asname or alias.name)

    return defined, imported


def _collect_assigned_names(target: ast.AST) -> set[str]:
    out: set[str] = set()
    if isinstance(target, ast.Name):
        out.add(target.id)
    elif isinstance(target, (ast.Tuple, ast.List)):
        for elt in target.elts:
            out |= _collect_assigned_names(elt)
    # ignore attributes/subscripts; they don't bind module names
    return out


def _class_to_children(cls: type, /, sort_children: bool, exclude: Iterable) -> list[Node]:
    """
    Classify class members using inspect.classify_class_attrs, following the user's
    original filters:
      - skip dunder names
      - only include members whose defining_class is the class itself

    Kinds are mapped to:
      'method', 'class_method', 'static_method', 'property', 'data'
    """
    children: list[Node] = []
    try:
        members = inspect.classify_class_attrs(cls)
    except:
        return children

    for m in members:
        if m.defining_class is cls:
            fqn = _fqn_for_class_member(cls, m.name, m.object)
            if _check_exclusions(fqn, exclude):
                continue
            children.append(Node(m.name, m.kind.replace(" ", "_")))

    if sort_children:
        children.sort(key=lambda d: (_KIND_ORDER.get(d.kind, 99), d.name))
    return children


def _module_to_node(module: ModuleType, /, sort_children: bool, exclude: Iterable) -> Node:
    """
    Convert a module/package object to a {name, kind, children} node.
    """
    if _check_exclusions(module.__name__, exclude):
        raise ValueError(f"Root module ``{module.__name__}`` is excluded")
    name = module.__name__.split(".")[-1]
    children: list[Node] = []
    kind = "package" if _is_package(module) else "module"

    # Subpackages / submodules (immediate)
    if kind == "package":
        for child_name, is_pkg in _iter_immediate_children(module):
            fqn = f"{module.__name__}.{child_name}"
            try:
                child_mod = importlib.import_module(fqn)
            except Exception:
                # Skip broken/optional imports
                continue
            if _check_exclusions(fqn, exclude):
                continue
            children.append(_module_to_node(child_mod, sort_children, exclusions))

    # Functions / Classes / Exceptions (module-defined only)
    for member_name, value in inspect.getmembers(module):
        vmod = inspect.getmodule(value)

        # CLASSES (incl. exceptions)
        if inspect.isclass(value) and vmod is module:
            fqn = _fqn_for_module_member(module, member_name, value)
            if _check_exclusions(fqn, exclude):
                continue
            cls_children = _class_to_children(value, sort_children, exclusions)
            if issubclass(value, Exception):
                children.append(Node(member_name, "exception", cls_children))
            else:
                children.append(Node(member_name, "class", cls_children))

        # FUNCTIONS
        elif inspect.isfunction(value) and vmod is module:
            fqn = _fqn_for_module_member(module, member_name, value)
            if _check_exclusions(fqn, exclude):
                continue
            children.append(Node(member_name, "function"))

    # ATTRIBUTES (data) — top-level assigned names that aren't functions/classes/modules
    # Use module dict to read current values, but gate by AST-defined names and not imported
    # Members defined in this module (exclude imports)
    defined_names, imported_names = _defined_and_imported_names_via_ast(module)
    for member_name in defined_names - imported_names:
        if hasattr(module, member_name):
            value = getattr(module, member_name)
            if not inspect.ismodule(value) and not inspect.isclass(value) and not inspect.isfunction(value):
                fqn = _fqn_for_module_member(module, member_name, value)
                if _check_exclusions(fqn, exclude):
                    continue
                children.append(Node(member_name, "attribute"))

    # Deterministic child ordering (optional)
    if sort_children:
        children.sort(key=lambda d: (_KIND_ORDER.get(d.kind, 99), d.name))

    return Node(name, kind, children)


# ------------------------------------------------------------------------------
# Project Information
# ------------------------------------------------------------------------------
project = "mapflpy"
author = "Predictive Science Inc"
copyright = f"{datetime.now():%Y}, {author}"
version = mapflpy.__version__
release = mapflpy.__version__

# ------------------------------------------------------------------------------
# General Configuration
# ------------------------------------------------------------------------------
extensions = []

# --- HTML Theme
_logo = "https://www.predsci.com/corona/apr2024eclipse/images/psi_logo.png"
html_favicon = _logo
html_logo = _logo
html_theme = "sphinx_book_theme"
html_static_path = ["_static"]
html_theme_options = {
    "show_prev_next": False,
    "navigation_with_keys": False,
    "show_nav_level": 3,
    "navigation_depth": 5,
    "logo": {
        "text": f"{project} v{version}",
        "image_light": _logo,
        "image_dark": _logo,
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

# --- Python Syntax
add_module_names = False
python_maximum_signature_line_length = 80

# --- Templating
templates_path = ['_templates', ]

# ------------------------------------------------------------------------------
# Viewcode Configuration
# ------------------------------------------------------------------------------
extensions.append("sphinx.ext.viewcode")

viewcode_line_numbers = True

# ------------------------------------------------------------------------------
# Autosummary Configuration
# ------------------------------------------------------------------------------
extensions.append("sphinx.ext.autosummary")

root_package = project
exclude_private = False
exclude_tests = True
exclude_dunder = True
sort_members = True
exclusions = ['\\._abc_impl',
              'mapflpy.fortran',
              '_load_array_to_shared_memory',
              'scripts._.*',
              'data\\.[A-Z_]*$'
              ]

exclusions.extend(_get_base_exclusions(exclude_private, exclude_tests, exclude_dunder))
node_tree = build_node_tree(root_package, sort_members, exclusions)

autosummary_context = dict(pkgtree=asdict(node_tree))

# ------------------------------------------------------------------------------
# Autodoc Configuration
# ------------------------------------------------------------------------------
extensions.append("sphinx.ext.autodoc")

autodoc_typehints = "description"
autodoc_member_order = 'bysource'
autodoc_default_options = {
    "show-inheritance": True,
}
autodoc_type_aliases = {
    "NumberType": "NumberType",
    "PathType": "PathType",
    "ArrayType": "ArrayType",
    "MagneticFieldArrayType": "MagneticFieldArrayType",
    "DirectionType": "DirectionType",
    "MagneticFieldLabelType": "MagneticFieldLabelType",
    "ContextType": "ContextType",
}

# ------------------------------------------------------------------------------
# Napoleon Configuration
# ------------------------------------------------------------------------------
extensions.append('sphinx.ext.napoleon')

napoleon_use_ivar = True
napoleon_preprocess_types = True
napoleon_type_aliases = {
    "NumberType": "~mapflpy.globals.NumberType",
    "PathType": "~mapflpy.globals.PathType",
    "ArrayType": "~mapflpy.globals.ArrayType",
    "MagneticFieldArrayType": "~mapflpy.globals.MagneticFieldArrayType",
    "DirectionType": "~mapflpy.globals.DirectionType",
    "MagneticFieldLabelType": "~mapflpy.globals.MagneticFieldLabelType",
    "ContextType": "~mapflpy.globals.ContextType",
    "Traces": "~mapflpy.globals.Traces",
    "Polarity": "~mapflpy.globals.Polarity",
    "MagneticFieldFiles": "~mapflpy.data.MagneticFieldFiles",
    "ndarray": "~numpy.ndarray",
    "ChainMap": "~collections.ChainMap",
    "Callable": "~collections.abc.Callable",
    "MutableMapping": "~collections.abc.MutableMapping",
}

# ------------------------------------------------------------------------------
# Intersphinx Configuration
# ------------------------------------------------------------------------------
extensions.append("sphinx.ext.intersphinx")

DOCS = Path(__file__).resolve().parents[1]
INV = DOCS / "_intersphinx"
intersphinx_cache_limit = 30
intersphinx_mapping = {
    "python": (
        "https://docs.python.org/3/",
        (INV / "python-objects.inv").as_posix(),
    ),
    "numpy": (
        "https://numpy.org/doc/stable/",
        (INV / "numpy-objects.inv").as_posix(),
    ),
    "scipy": (
        "https://docs.scipy.org/doc/scipy/reference/",
        (INV / "scipy-objects.inv").as_posix(),
    ),
    "matplotlib": (
        "https://matplotlib.org/stable/",
        (INV / "matplotlib-objects.inv").as_posix(),
    ),
    "pooch": (
        "https://www.fatiando.org/pooch/latest/",
        (INV / "pooch-objects.inv").as_posix(),
    ),
}

# ------------------------------------------------------------------------------
# Sphinx-Gallery Configuration
# ------------------------------------------------------------------------------
extensions.append("sphinx_gallery.gen_gallery")

import matplotlib
matplotlib.use("Agg")
os.environ.setdefault('SPHINX_GALLERY_BUILD', '1')

sphinx_gallery_conf = {
    "examples_dirs": ["../../examples"],
    "gallery_dirs": ["gallery"],
    "within_subsection_order": "FileNameSortKey",
    "download_all_examples": False,
    "remove_config_comments": True,
    "filename_pattern": r"\.py$",
    # Run examples during build (True). For quick local builds, set False.
    "plot_gallery": True,
    "matplotlib_animations": True,
}

# from collections.abc import Mapping, Sequence
# from enum import Enum
# from types import MappingProxyType
# from pathlib import Path
# import importlib, json, pprint
#
# def _resolve_target(name, fallback):
#     try:
#         modname, attr = name.rsplit(".", 1)
#         mod = importlib.import_module(modname)
#         return getattr(mod, attr, fallback)
#     except Exception:
#         return fallback
#
# def _to_jsonable(obj, _seen: set[int] | None = None):
#     """
#     Convert arbitrary Python objects into something json.dumps() can handle.
#     - Mappings -> dict with stringified keys
#     - Sequences -> list
#     - Sets -> sorted list of repr()
#     - Enums -> .value (then recurse)
#     - Path -> str
#     - NumPy scalars -> Python scalars (duck-typed)
#     - Fallback -> repr(obj)
#     """
#     if _seen is None:
#         _seen = set()
#     oid = id(obj)
#     if oid in _seen:
#         return ""
#     _seen.add(oid)
#
#     # Enum
#     if isinstance(obj, Enum):
#         return _to_jsonable(obj.value, _seen)
#
#     # Pathlike
#     if isinstance(obj, Path):
#         return str(obj)
#
#     # Mapping (incl. mappingproxy)
#     if isinstance(obj, Mapping):
#         out = {}
#         for k, v in obj.items():
#             # JSON keys must be strings
#             ks = k if isinstance(k, str) else repr(k)
#             out[ks] = _to_jsonable(v, _seen)
#         return out
#
#     # Sequence but not str/bytes
#     if isinstance(obj, Sequence) and not isinstance(obj, (str, bytes, bytearray)):
#         return [_to_jsonable(x, _seen) for x in obj]
#
#     # Sets
#     if isinstance(obj, (set, frozenset)):
#         return sorted((repr(x) for x in obj))
#
#     # NumPy scalar duck-typing (avoid importing numpy here)
#     try:
#         # numpy scalars have .item() to get Python value
#         item = getattr(obj, "item", None)
#         if callable(item):
#             return item()
#     except Exception:
#         pass
#
#     # Basic types pass through
#     if isinstance(obj, (type(None), bool, int, float, str)):
#         return obj
#
#     # Last resort: repr
#     return repr(obj)
#
# def _pp_lines(mapping: Mapping) -> list[str]:
#     # Try JSON first
#     try:
#         safe = _to_jsonable(mapping)
#         txt = json.dumps(safe, indent=2, sort_keys=True, ensure_ascii=False)
#         return [".. code-block:: json", "", *("  " + ln for ln in txt.splitlines())]
#     except Exception:
#         # Fallback to pretty-printed Python
#         txt = pprint.pformat(mapping, width=88, sort_dicts=True, compact=False)
#         return [".. code-block:: python", "", *("  " + ln for ln in txt.splitlines())]
#
# def pretty_print_data(app, what, name, obj, options, lines):
#     # Cover variables exposed as either 'data' (module var) or 'attribute' (re-export)
#     if what not in {"data", "attribute"}:
#         return
#
#     target = _resolve_target(name, obj)
#
#     # Accept any Mapping, including MappingProxyType
#     if isinstance(target, Mapping) or isinstance(target, MappingProxyType):
#         lines[:] = ["Rendered value:", "", *_pp_lines(target)]
#
# def setup(app):
#     app.connect("autodoc-process-docstring", pretty_print_data)
#     return {"version": "1.1", "parallel_read_safe": True}
