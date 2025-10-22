#!/usr/bin/env bash

rm -rf _build
rm -rf source/autodoc
rm -rf source/gallery

#rm -rf _intersphinx
#./fetch.sh

make clean html
