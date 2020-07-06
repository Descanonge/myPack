#! /usr/bin/env bash

file=$1

build="_build_ds"
# flags="-fbounds-check"
flags=""


f2py -c $file.f90 -m $file \
    --build-dir "$build" \
    --f90flags="$flags" \
    --verbose
