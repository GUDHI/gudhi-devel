#!/bin/bash

set -ue

export CGAL_VERSION="6.0.1"

curl -LO "https://github.com/CGAL/cgal/releases/download/v${CGAL_VERSION}/CGAL-${CGAL_VERSION}.tar.xz"
tar xf CGAL-${CGAL_VERSION}.tar.xz
cd CGAL-${CGAL_VERSION}
cmake -DCMAKE_BUILD_TYPE=Release .
make install
cd ..
rm -rf CGAL-${CGAL_VERSION}
