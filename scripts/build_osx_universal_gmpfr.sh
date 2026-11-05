#!/bin/bash
set -e

# In the working directory, creates deps-uni/lib/*

# Downloading deps-amd64 (x86_64 gmp and mpfr) and deps-arm64 (arm64 gmp and mpfr)
# 2nd July 2026 version from ventura brew bottles - broken from brew version 6.X
curl -OL "https://gudhi.inria.fr/gmp_mpfr_ventura_brew_bottles.tgz"
curl -OL "https://gudhi.inria.fr/gmp_mpfr_ventura_brew_bottles.tgz.sha512"
shasum -a 512 -c *.sha512
tar xf gmp_mpfr_ventura_brew_bottles.tgz

# Merging
mkdir -p deps-uni/lib
GMP1=deps-amd64/gmp/*/lib/libgmp.*.dylib
GMP=`basename $GMP1`
GMPXX1=deps-amd64/gmp/*/lib/libgmpxx.*.dylib
GMPXX=`basename $GMPXX1`
MPFR1=deps-amd64/mpfr/*/lib/libmpfr.*.dylib
MPFR=`basename $MPFR1`
lipo -create $GMP1 deps-arm64/gmp/*/lib/$GMP -output deps-uni/lib/$GMP
lipo -create $GMPXX1 deps-arm64/gmp/*/lib/$GMPXX -output deps-uni/lib/$GMPXX
lipo -create $MPFR1 deps-arm64/mpfr/*/lib/$MPFR -output deps-uni/lib/$MPFR

# Necessary even for libs created by lipo
install_name_tool -id $PWD/deps-uni/lib/$GMP deps-uni/lib/$GMP
install_name_tool -id $PWD/deps-uni/lib/$GMPXX deps-uni/lib/$GMPXX
install_name_tool -id $PWD/deps-uni/lib/$MPFR deps-uni/lib/$MPFR
# Also fix dependencies
# otool gives twice the same dependency, keep only one (a loop would be safer...)
BADGMP=`otool -L deps-uni/lib/$MPFR|sed -ne 's/[[:space:]]*\(.*libgmp\..*dylib\).*/\1/p'|uniq`
install_name_tool -change $BADGMP $PWD/deps-uni/lib/$GMP deps-uni/lib/$MPFR
BADGMP=`otool -L deps-uni/lib/$GMPXX|sed -ne 's/[[:space:]]*\(.*libgmp\..*dylib\).*/\1/p'|uniq`
install_name_tool -change $BADGMP $PWD/deps-uni/lib/$GMP deps-uni/lib/$GMPXX

ln -s $GMP deps-uni/lib/libgmp.dylib
ln -s $GMPXX deps-uni/lib/libgmpxx.dylib
ln -s $MPFR deps-uni/lib/libmpfr.dylib

# Debug
ls -l deps-uni/lib
otool -L deps-uni/lib/*.*.dylib
