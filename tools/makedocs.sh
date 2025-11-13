#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PARENT_DIR="$(dirname "$SCRIPT_DIR")"

fetch=0
remove=0
while getopts ":fr" name
do
   case $name in
     f)    fetch=1;;
     r)    remove=1;;
     \?) printf 'Usage: %s [-f] [-r]\n' "$0"; exit 2 ;;
   esac
done

if (( remove )); then
  rm -rf "$PARENT_DIR/docs/_build" \
         "$PARENT_DIR/docs/source/autodoc" \
         "$PARENT_DIR/docs/source/gallery"
fi

if (( fetch )); then
  rm -rf "$PARENT_DIR/docs/_intersphinx"
  chmod +x "$PARENT_DIR/docs/fetch.sh"
  "$PARENT_DIR/docs/fetch.sh"
fi

make -C "$PARENT_DIR/docs" html
