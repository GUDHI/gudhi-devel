#!/bin/bash

set -eu

brew update
brew install boost eigen gmp mpfr cgal
python -m pip install --user numpy~=1.21.4
python -m pip install --user -r ext/gudhi-deploy/build-requirements.txt
python -m pip install --user twine delocate
./scripts/build_osx_universal_gmpfr.sh
export   GMP_LIB_DIR=$PWD/deps-uni/lib
export GMPXX_LIB_DIR=$PWD/deps-uni/lib
export  MPFR_LIB_DIR=$PWD/deps-uni/lib
