#!/usr/bin/env python3
from __future__ import annotations
import argparse
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

# ---- small helpers

def run(cmd: list[str], cwd: Path | None = None, env: dict[str, str] | None = None) -> None:
    where = f" (cwd={cwd})" if cwd else ""
    print("+", " ".join(cmd) + where, flush=True)
    subprocess.run(cmd, cwd=str(cwd) if cwd else None, env=env, check=True)

def read_pyproject_name(pp: Path) -> str | None:
    try:
        try:
            import tomllib  # 3.11+
        except ModuleNotFoundError:  # pragma: no cover
            import tomli as tomllib  # type: ignore
        data = tomllib.loads(pp.read_text(encoding="utf-8"))
        proj = data.get("project", {})
        name = proj.get("name")
        return name if isinstance(name, str) else None
    except Exception:
        return None

def select_best_wheel(dist_dir: Path) -> Path:
    wheels = sorted(dist_dir.glob("*.whl"))
    if not wheels:
        raise SystemExit(f"No wheels found in {dist_dir}")
    # Prefer the wheel that best matches current interpreter & platform
    try:
        from packaging import tags
        supported = list(tags.sys_tags())
        def score(p: Path) -> int:
            n = p.name
            for i, t in enumerate(supported):
                tag = f"{t.interpreter}-{t.abi}-{t.platform}"
                # allow -any platform too
                if tag in n or f"{t.interpreter}-{t.abi}-any" in n:
                    return -i
            return 10**9
        return sorted(wheels, key=score)[0]
    except Exception:
        # Fallback: newest by mtime
        return max(wheels, key=lambda p: p.stat().st_mtime)

def extract_extension_from_wheel(wheel: Path, dest_pkg_dir: Path) -> Path:
    # Discouraged pattern—only if you truly need the .so in-tree (e.g. docs hack).
    with zipfile.ZipFile(wheel) as zf:
        members = [m for m in zf.namelist()
                   if m.replace("\\", "/").startswith(f"{dest_pkg_dir.name}/fortran/")
                   and m.endswith((".so", ".dylib", ".pyd"))]
        if not members:
            raise SystemExit("No compiled extension found inside wheel.")
        member = members[0]
        out = dest_pkg_dir / "fortran" / Path(member).name
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_bytes(zf.read(member))
        print(f"Extracted {member} -> {out}")
        return out

# ---- main script

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build wheel (and sdist), install into a fresh venv, run tests. Optionally extract the compiled extension back to source."
    )
    parser.add_argument("-c", "--clean", action="store_true", help="Remove build artifacts before building")
    parser.add_argument("-s", "--sdist", action="store_true", help="Also build an sdist")
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1],
                        help="Project root (default: the parent of this script)")
    parser.add_argument("--package-name", default=None,
                        help="Package name to install from wheel (default: read from pyproject.toml)")
    parser.add_argument("--tests", type=Path, default=None,
                        help="Path to test directory (default: <root>/tests if exists)")
    parser.add_argument("--extract", action="store_true",
                        help="Extract the compiled extension from the wheel back into the source tree")
    parser.add_argument("--pytest-args", nargs=argparse.REMAINDER,
                        help="Extra args to pass to pytest (use after `--`)")
    args = parser.parse_args()

    root = args.root.resolve()
    dist = root / "dist"
    build_dir = root / "build"
    egg_info = next(root.glob("*.egg-info"), None)
    pkg_name = args.package_name or read_pyproject_name(root / "pyproject.toml") or ""
    if not pkg_name:
        raise SystemExit("Cannot determine package name; use --package-name or set [project].name in pyproject.toml.")

    tests_dir = args.tests or (root / "tests" if (root / "tests").is_dir() else None)

    if args.clean:
        for p in (dist, build_dir, root / ".pytest_cache"):
            if p.exists():
                print(f"Removing {p}")
                shutil.rmtree(p, ignore_errors=True)
        if egg_info and egg_info.exists():
            print(f"Removing {egg_info}")
            shutil.rmtree(egg_info, ignore_errors=True)

    # Build artifacts
    build_cmd = [sys.executable, "-m", "build", "--wheel"]
    if args.sdist:
        build_cmd.append("--sdist")
    build_cmd.append(str(root))
    run(build_cmd)

    wheel = select_best_wheel(dist)
    print(f"Selected wheel: {wheel.name}")

    # Optional: extract compiled extension back into source (not recommended)
    if args.extract:
        extract_extension_from_wheel(wheel, root / pkg_name)

    # Run tests (if any)
    if tests_dir and tests_dir.exists():
        cmd = [sys.executable, "-m", "pytest", str(tests_dir)]
        if args.pytest_args:
            cmd = [*cmd, *args.pytest_args]
        run(cmd)
    else:
        print("No tests directory found; skipping pytest.")

    print("All done.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
