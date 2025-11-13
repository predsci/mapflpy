#!/bin/bash

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd -P)"
PARENT_DIR="$(dirname "$SCRIPT_DIR")"
TARGET_DIR="$PARENT_DIR/docs/_intersphinx"
if [ ! -d "$TARGET_DIR" ]; then
    mkdir -p "$TARGET_DIR"
fi
# this script updates the intersphinx files here
# make sure to follow potential redirects
curl -L https://docs.python.org/3/objects.inv > "${TARGET_DIR}"/python-objects.inv
curl -L https://docs.scipy.org/doc/scipy/objects.inv > "${TARGET_DIR}"/scipy-objects.inv
curl -L https://numpy.org/doc/stable/objects.inv > "${TARGET_DIR}"/numpy-objects.inv
curl -L https://matplotlib.org/stable/objects.inv > "${TARGET_DIR}"/matplotlib-objects.inv
curl -L https://docs.pytest.org/en/stable/objects.inv > "${TARGET_DIR}"/pytest-objects.inv
curl -L https://www.fatiando.org/pooch/latest/objects.inv > "${TARGET_DIR}"/pooch-objects.inv