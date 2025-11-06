# noxfile.py
from __future__ import annotations
import os
import platform
from pathlib import Path
import nox

# Speed up reruns; reuse envs when deps unchanged
nox.options.reuse_existing_virtualenvs = True


PY_VERSIONS = ["3.10", "3.11", "3.12", "3.13"]
PY_DEFAULT = "3.12"

ROOT = Path(__file__).parent.resolve()
ARTIFACTS = ROOT / ".nox" / "_artifacts"
DIST = ARTIFACTS / "dist"
WHEELHOUSE = ARTIFACTS / "wheelhouse"

pyproject = nox.project.load_toml("pyproject.toml")
PROJECT = pyproject["project"]["name"]
CONDA_ENV_BUILD_COMPILERS = pyproject["tool"][PROJECT].get("conda", [])

REPAIR_TOOLS: dict[str, list[str]] = {
    "linux": ["auditwheel"],
    "darwin": ["delocate"],
    "windows": ["delvewheel"],
}

DIST.mkdir(parents=True, exist_ok=True)
WHEELHOUSE.mkdir(parents=True, exist_ok=True)


def _darwin_sdk_env() -> dict[str, str]:
    """macOS: provide SDK + baseline target so the Fortran probe succeeds."""
    if platform.system() != "Darwin":
        return {}
    # Prefer already-set values; otherwise best-effort defaults
    env = {}
    env.setdefault("MACOSX_DEPLOYMENT_TARGET", os.environ.get("MACOSX_DEPLOYMENT_TARGET", "11.0"))
    # SDKROOT may be needed for clang/gfortran during Meson sanity checks
    if "SDKROOT" not in os.environ:
        try:
            import subprocess
            sdk = subprocess.check_output(["xcrun", "--show-sdk-path"], text=True).strip()
            env["SDKROOT"] = sdk
        except Exception:
            pass
    return env


def _build_env(session: nox.Session) -> Path:
    """Build a wheel into ./dist and return its path."""
    session.conda_install(
        *CONDA_ENV_BUILD_COMPILERS,
        channel="conda-forge"
    )
    session.env.update(_darwin_sdk_env())


def _dist_env(session: nox.Session) -> Path:
    """Environment for installing from built wheels."""
    session.install(
        *pyproject["project"].get("dependencies", []),
    )

    session.run(
        "python", "-m", "pip", "install",
        "--no-index", f"--find-links={WHEELHOUSE}",
        "--only-binary=:all:",               # avoid accidentally picking an sdist
        PROJECT,
    )


@nox.session(venv_backend='conda', python=PY_VERSIONS)
def build(session: nox.Session) -> None:
    """Build the package wheel (with compilers)."""
    _build_env(session)
    session.conda_install(
        *pyproject["build-system"].get("requires", []),
        *pyproject["project"].get("optional-dependencies", {}).get("build", []),
        channel="conda-forge",
    )
    session.run(
        "python", "-m", "build",
        "--wheel", "--outdir", DIST.as_posix(),
        external=False
    )

@nox.session(python=PY_DEFAULT)
def repair(session: nox.Session) -> None:
    """Repair wheels in dist/ into wheelhouse/ using the OS-specific tool."""
    platform_id = platform.system().lower()
    wheels = sorted(DIST.glob("*.whl"))

    match platform_id:
        case "linux":
            session.install("auditwheel")
            for whl in wheels:
                session.run("auditwheel", "show", str(whl))
                session.run("auditwheel", "repair", "-w", str(WHEELHOUSE), str(whl))
        case "darwin":
            session.install("delocate")
            for whl in wheels:
                session.run("delocate-listdeps", str(whl))
                session.run("delocate-wheel", "-w", str(WHEELHOUSE), str(whl))
        case "windows":
            session.install("delvewheel")
            for whl in wheels:
                session.run("python", "-m", "delvewheel", "show", str(whl))
                session.run("python", "-m", "delvewheel", "repair", "-w", str(WHEELHOUSE), str(whl))


@nox.session(python=PY_VERSIONS)
def test(session: nox.Session) -> None:
    """Build the wheel (with compilers), install it, then run pytest from a temp dir."""
    # Build wheel
    _dist_env(session)

    # Runtime/test deps
    session.install(*pyproject["project"].get("optional-dependencies", {}).get("test", []))

    tmp = session.create_tmp()
    session.chdir(tmp)

    # Pytest
    session.run(
        "python", "-m", "pytest",
        str(ROOT),
        "--cov=mapflpy", "--cov-report=term-missing",
        "--maxfail=1", "-q",
        "--import-mode=importlib",
    )

@nox.session(python=PY_DEFAULT)
def types(session: nox.Session) -> None:
    """Mypy type checking (analyzes source tree)."""
    session.install(*pyproject["project"].get("optional-dependencies", {}).get("types", []))

    session.run("mypy")

@nox.session(python=PY_DEFAULT)
def lint(session: nox.Session) -> None:
    """Ruff lint + format check."""
    session.install(*pyproject["project"].get("optional-dependencies", {}).get("lint", []))

    session.run("ruff", "check", PROJECT)
    session.run("ruff", "format", "--check", PROJECT)

# @nox.session(venv_backend="conda", python=PY)
# def docs(session: nox.Session) -> None:
#     """Build Sphinx docs against the installed wheel."""
#     # Build wheel (same as tests) and install
#     whl = _ensure_wheel(session)
#     session.install(str(whl))
#
#     # Docs deps (conda for heavy libs; pip for Python pkgs)
#     session.conda_install("--channel", "conda-forge", "numpy", "matplotlib", "hdf5")
#     session.install("sphinx", "sphinx-book-theme", "sphinx-gallery", "numpydoc", "pooch")
#
#     session.chdir(ROOT / "docs")
#     session.run("sphinx-build", "-b", "html", "source", "_build/html", "-W")
