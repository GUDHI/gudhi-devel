#!/bin/bash

set -eu

apk update
# cmake 3.31.6 is ok here
apk add curl cmake

./install_boost_from_curl_download.sh
./install_gmp_from_curl_download.sh
./install_mpfr_from_curl_download.sh
./install_cgal_from_curl_download.sh
