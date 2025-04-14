#!/bin/bash

set -eu

brew update
brew install curl cmake eigen

scripts/install_boost_from_curl_download.sh
scripts/install_gmp_from_curl_download.sh
scripts/install_mpfr_from_curl_download.sh
scripts/install_cgal_from_curl_download.sh
