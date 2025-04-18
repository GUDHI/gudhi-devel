#!/bin/bash

set -eu

# Install recent boost headers to take advantage of the latest improvements - no dll required
curl -LO "https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz"
tar xf eigen-3.4.0.tar.gz
cd eigen-3.4.0
mkdir build
cd build
cmake ..
make all install
cd ../..
rm -rf eigen-3.4.0
