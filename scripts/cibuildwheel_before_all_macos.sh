#!/bin/bash

set -eu

brew update
brew install curl cmake

scripts/install_eigen_from_curl_download.sh
scripts/install_boost_from_curl_download.sh
scripts/install_gmp_from_curl_download.sh
scripts/install_mpfr_from_curl_download.sh
scripts/install_cgal_from_curl_download.sh
