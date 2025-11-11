#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
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
  rm -rf "$SCRIPT_DIR/docs/_build" \
         "$SCRIPT_DIR/docs/source/autodoc" \
         "$SCRIPT_DIR/docs/source/gallery"
fi

if (( fetch )); then
  rm -rf "$SCRIPT_DIR/docs/_intersphinx"
  chmod +x "$SCRIPT_DIR/docs/fetch.sh"
  "$SCRIPT_DIR/docs/fetch.sh"
fi

make -C "$SCRIPT_DIR/docs" html
