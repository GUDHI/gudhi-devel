#!/bin/bash

set -eu

brew update
brew install curl cmake eigen

./install_boost_from_curl_download.sh
./install_gmp_from_curl_download.sh
./install_mpfr_from_curl_download.sh
./install_cgal_from_curl_download.sh
