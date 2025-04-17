#!/bin/bash

set -eu

# Requires execution rights
chmod +x scripts/install_*.sh

if command -v yum 2>&1 >/dev/null
then
  echo "yum found"
  yum -y update
  # manylinux_2_28 means:
  #  - cmake is 3.26.5
  #  - eigen3-devel (header only) 3.3.4-6.el8
  #  - gmp-devel 6.1.2-11.el8
  #  - mpfr-devel 3.1.6-1.el8
  yum -y install curl cmake eigen3-devel gmp-devel mpfr-devel
  # yum install boost-devel would be 1.66.0
else
  echo $(uname -a)
  echo "No yum available"
  exit 1
fi

scripts/install_boost_from_curl_download.sh
# Install CGAL manually to detect eigen boost gmp and mpfr
scripts/install_cgal_from_curl_download.sh
