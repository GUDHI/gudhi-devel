#!/bin/bash

set -eu

# Install recent boost headers to take advantage of the latest improvements - no dll required
curl -LO "https://archives.boost.io/release/1.87.0/source/boost_1_87_0.tar.gz"
tar xf boost_1_87_0.tar.gz
cd boost_1_87_0
./bootstrap.sh
./b2 install
cd ..
rm -rf boost_1_87_0
