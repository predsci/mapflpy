#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PARENT_DIR="$(dirname "$SCRIPT_DIR")"

rm -rf "$PARENT_DIR/dist" \
        "$PARENT_DIR/build" \
        "$PARENT_DIR/"*.egg-info \
        "$PARENT_DIR/mapflpy_fortran-f2pywrappers2.f90" \
        "$PARENT_DIR/mapflpy_fortranmodule.c" \
        "$PARENT_DIR/mapflpy/fortran/"*.so \
        "$PARENT_DIR/.pytest_cache" \

python -m build "$PARENT_DIR"

# pick the wheel
wheel="$(ls "$PARENT_DIR"/dist/*.whl | head -n1)"

# list members and select the .so you want
member="$(unzip -Z1 "$wheel" 'mapflpy/fortran/*.so' | head -n1)"

# extract that single file to the desired path
unzip -p "$wheel" "$member" > "$PARENT_DIR/mapflpy/fortran/$(basename "$member")"

rm -rf "$PARENT_DIR/dist" \
        "$PARENT_DIR/build" \
        "$PARENT_DIR/"*.egg-info \
        "$PARENT_DIR/mapflpy_fortran-f2pywrappers2.f90" \
        "$PARENT_DIR/mapflpy_fortranmodule.c"

python -m pytest "$PARENT_DIR/tests"