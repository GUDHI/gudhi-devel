#!/bin/bash

set -eu

# Requires execution rights
chmod +x scripts/install_*.sh

if command -v yum 2>&1 >/dev/null
then
  echo "yum found"
  yum -y update
  # eigen (header only) 3.3.7-1.el7 is ok here - not available on manylinux2014_i686
  yum -y install curl openssl-devel eigen3-devel
  # yum install cmake would be 2.8.12
  scripts/install_cmake_from_curl_download.sh
  # yum install boost would be 1.53.0
elif command -v apk 2>&1 >/dev/null
then
  echo "apk found"
  apk update
  # cmake 3.31.6 is ok here
  # eigen (header only) 3.4.0-r10 is ok here
  apk add curl cmake eigen
  # boost-dev (header only for python package) 1.84 would be here, but let's go with 1.87
else
  echo $(uname -a)
  echo "No apk nor yum available"
  exit 1
fi

scripts/install_boost_from_curl_download.sh
# Need to compile gmp and mpfr for wheel repair purpose
scripts/install_gmp_from_curl_download.sh
scripts/install_mpfr_from_curl_download.sh
# Install CGAL manually to detect eigen boost gmp and mpfr
scripts/install_cgal_from_curl_download.sh
