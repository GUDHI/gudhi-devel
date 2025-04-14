#!/bin/bash

set -ue

export MPFR_VERSION="4.2.2"

curl -LO "https://ftp.gnu.org/gnu/mpfr/mpfr-${MPFR_VERSION}.tar.gz"
tar xf mpfr-${MPFR_VERSION}.tar.gz
cd mpfr-${MPFR_VERSION}
./configure
make
make install
cd ..
rm -rf mpfr-${MPFR_VERSION}
