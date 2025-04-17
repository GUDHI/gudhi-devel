#!/bin/bash

set -eu

brew update
brew install boost eigen gmp mpfr cgal
# For python>=3.9, numpy>=2.0 for package build and ABI compatibility with numpy 1.X and 2.X
# cf. https://numpy.org/doc/stable/dev/depending_on_numpy.html#numpy-2-0-specific-advice
python -m pip install --user numpy>=2.0
python -m pip install --user -r ext/gudhi-deploy/build-requirements.txt
python -m pip install --user twine delocate
./scripts/build_osx_universal_gmpfr.sh
export   GMP_LIB_DIR=$PWD/deps-uni/lib
export GMPXX_LIB_DIR=$PWD/deps-uni/lib
export  MPFR_LIB_DIR=$PWD/deps-uni/lib
