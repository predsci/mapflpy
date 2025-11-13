# noxfile.py
from __future__ import annotations

import os
import json
import platform
import subprocess
from pathlib import Path
import nox

# Speed up reruns; reuse envs when deps unchanged
nox.options.reuse_existing_virtualenvs = False


PY_VERSIONS = ["3.10", "3.11", "3.12", "3.13"]
PY_TARGET = PY_VERSIONS[-1]

PROJECT_NAME_PATH = Path(__file__).parent.resolve()
ARTIFACTS = PROJECT_NAME_PATH / ".nox" / "_artifacts"
WHEEL_DIR = ARTIFACTS / "wheels"
SDIST_DIR = ARTIFACTS / "sdist"
WHEELHOUSE_DIR = ARTIFACTS / "wheelhouse"
DOCDIST_DIR = ARTIFACTS / "docs"

pyproject = nox.project.load_toml("pyproject.toml")
PROJECT_NAME = pyproject["project"]["name"]
CONDA_ENV_BUILD_COMPILERS = pyproject["tool"][PROJECT_NAME].get("conda", [])

REPAIR_TOOLS: dict[str, list[str]] = {
    "linux": ["auditwheel"],
    "darwin": ["delocate"],
    "windows": ["delvewheel"],
}

WHEEL_DIR.mkdir(parents=True, exist_ok=True)
WHEELHOUSE_DIR.mkdir(parents=True, exist_ok=True)
DOCDIST_DIR.mkdir(parents=True, exist_ok=True)


def _darwin_sdk_env() -> dict[str, str]:
    """macOS: provide SDK + baseline target so the Fortran probe succeeds."""
    if platform.system() != "Darwin":
        return {}
    # Prefer already-set values; otherwise best-effort defaults
    env = {}
    env.setdefault("MACOSX_DEPLOYMENT_TARGET", os.environ.get("MACOSX_DEPLOYMENT_TARGET", "11.0"))
    # SDKPROJECT_NAME_PATH may be needed for clang/gfortran during Meson sanity checks
    # if "SDKPROJECT_NAME_PATH" not in os.environ:
    #     try:
    #         import subprocess
    #         sdk = subprocess.check_output(["xcrun", "--show-sdk-path"], text=True).strip()
    #         env["SDKPROJECT_NAME_PATH"] = sdk
    #     except Exception:
    #         pass
    # return env
    try:
        sdk = subprocess.check_output(
            ["xcrun", "--sdk", "macosx", "--show-sdk-path"],
            text=True
        ).strip()
        env["SDKROOT"] = os.environ.get("SDKROOT", sdk)
        # Help the linkers/compilers see the SDK explicitly:
        env.setdefault("CFLAGS",    f"-isysroot {sdk} -mmacosx-version-min={env['MACOSX_DEPLOYMENT_TARGET']}")
        env.setdefault("CXXFLAGS",  f"-isysroot {sdk} -mmacosx-version-min={env['MACOSX_DEPLOYMENT_TARGET']}")
        env.setdefault("FCFLAGS",   f"-isysroot {sdk} -mmacosx-version-min={env['MACOSX_DEPLOYMENT_TARGET']}")
        env.setdefault("LDFLAGS",   f"-Wl,-syslibroot,{sdk} -mmacosx-version-min={env['MACOSX_DEPLOYMENT_TARGET']}")
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
        "--no-index", f"--find-links={WHEELHOUSE_DIR}",
        "--only-binary=:all:", "--report", "-",
        PROJECT_NAME,
    )


@nox.session(venv_backend='conda|mamba|micromamba', python=PY_VERSIONS)
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
        "--wheel", "--outdir", WHEEL_DIR.as_posix(),
        external=False
    )

@nox.session(python=PY_TARGET)
def repair(session: nox.Session) -> None:
    """Repair wheels in dist/ into wheelhouse/ using the OS-specific tool."""
    platform_id = platform.system().lower()
    wheels = sorted(WHEEL_DIR.glob("*.whl"))

    match platform_id:
        case "linux":
            session.install("auditwheel")
            for whl in wheels:
                session.run("auditwheel", "show", str(whl))
                session.run("auditwheel", "repair", "-w", str(WHEELHOUSE_DIR), str(whl))
        case "darwin":
            session.install("delocate")
            for whl in wheels:
                session.run("delocate-listdeps", str(whl))
                session.run("delocate-wheel", "-w", str(WHEELHOUSE_DIR), str(whl))
        case "windows":
            session.install("delvewheel")
            for whl in wheels:
                session.run("python", "-m", "delvewheel", "show", str(whl))
                session.run("python", "-m", "delvewheel", "repair", "-w", str(WHEELHOUSE_DIR), str(whl))


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
    session.run("pytest", PROJECT_NAME_PATH.as_posix())

@nox.session(venv_backend='conda|mamba|micromamba', python=PY_TARGET)
def sdist(session: nox.Session) -> None:
    """Build the package wheel (with compilers)."""
    _build_env(session)
    session.conda_install(
        *pyproject["build-system"].get("requires", []),
        *pyproject["project"].get("optional-dependencies", {}).get("build", []),
        channel="conda-forge",
    )
    session.run(
        "python", "-m", "build",
        "--sdist", "--outdir", SDIST_DIR.as_posix(),
        external=False
    )

@nox.session(python=PY_TARGET)
def types(session: nox.Session) -> None:
    """Mypy type checking (analyzes source tree)."""
    session.install(*pyproject["project"].get("optional-dependencies", {}).get("types", []))

    session.run("mypy")

@nox.session(python=PY_TARGET)
def lint(session: nox.Session) -> None:
    """Ruff lint + format check."""
    session.install(*pyproject["project"].get("optional-dependencies", {}).get("lint", []))

    session.run("ruff", "check", PROJECT_NAME)
    session.run("ruff", "format", "--check", PROJECT_NAME)

@nox.session(python=PY_TARGET)
def docs(session: nox.Session) -> None:
    """Build Sphinx docs against the installed wheel."""
    _dist_env(session)
    session.install(*pyproject["project"].get("optional-dependencies", {}).get("docs", []))
    args = pyproject["tool"].get("sphinx_build", {}).get("addopts", [])

    stamp = session.env_dir / ".mapflpy_tag.json"
    session.run(
        "python", "-c",
        (
            "import json; "
            "from importlib import metadata as md; "
            "d=md.distribution('mapflpy'); "
            "txt=d.read_text('WHEEL'); "
            "tags=[]; "
            "tags=[ln.split(':',1)[1].strip() for ln in txt.splitlines() if ln.startswith('Tag: ')]; "
            "out={'tags': tags}; "
            f"open('{stamp}','w').write(json.dumps(out))"
        )
    )
    tag = json.loads(Path(stamp).read_text())

    out_dir = DOCDIST_DIR / f"html-{tag.get('tags', ['none'])[0]}"
    src_dir = PROJECT_NAME_PATH / "docs" / "source"
    session.run("sphinx-build", src_dir.as_posix(), out_dir.as_posix(), *args)
