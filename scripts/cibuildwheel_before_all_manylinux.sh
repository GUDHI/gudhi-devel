#!/bin/bash

set -eu

yum -y update
yum -y install curl eigen3-devel

# yum install cmake would be 2.8.12
./install_cmake_from_curl_download.sh
./install_boost_from_curl_download.sh
./install_gmp_from_curl_download.sh
./install_mpfr_from_curl_download.sh
./install_cgal_from_curl_download.sh
