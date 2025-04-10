#!/bin/bash

set -ue

export CMAKE_VERSION="3.31"
export CMAKE_BUILD="6"

curl -LO "https://cmake.org/files/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}.${CMAKE_BUILD}.tar.gz"
tar xf cmake-${CMAKE_VERSION}.${CMAKE_BUILD}.tar.xz
cd cmake-${CMAKE_VERSION}.${CMAKE_BUILD}
./bootstrap
make
make install
