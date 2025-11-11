#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

rm -rf "$SCRIPT_DIR/dist" \
        "$SCRIPT_DIR/build" \
        "$SCRIPT_DIR/"*.egg-info \
        "$SCRIPT_DIR/mapflpy_fortran-f2pywrappers2.f90" \
        "$SCRIPT_DIR/mapflpy_fortranmodule.c" \
        "$SCRIPT_DIR/mapflpy/fortran/"*.so \
        "$SCRIPT_DIR/.pytest_cache" \

python -m build "$SCRIPT_DIR"

# pick the wheel
wheel="$(ls "$SCRIPT_DIR"/dist/*.whl | head -n1)"

# list members and select the .so you want
member="$(unzip -Z1 "$wheel" 'mapflpy/fortran/*.so' | head -n1)"

# extract that single file to the desired path
unzip -p "$wheel" "$member" > "$SCRIPT_DIR/mapflpy/fortran/$(basename "$member")"

rm -rf "$SCRIPT_DIR/dist" \
        "$SCRIPT_DIR/build" \
        "$SCRIPT_DIR/"*.egg-info \
        "$SCRIPT_DIR/mapflpy_fortran-f2pywrappers2.f90" \
        "$SCRIPT_DIR/mapflpy_fortranmodule.c"

python -m pytest "$SCRIPT_DIR/tests"