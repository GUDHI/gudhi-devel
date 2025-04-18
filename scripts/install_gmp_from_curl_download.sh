#!/bin/bash

set -ue

export GMP_VERSION="6.3.0"

curl -LO "https://gmplib.org/download/gmp/gmp-${GMP_VERSION}.tar.gz"
tar xf gmp-${GMP_VERSION}.tar.gz
cd gmp-${GMP_VERSION}
./configure
make
make install
cd ..
rm -rf gmp-${GMP_VERSION}
