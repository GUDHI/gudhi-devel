#!/bin/bash
set -e

# oldest OSx version that brew supports for gmp and mpfr 'bottles'
OSX_VERSION='ventura'

# In the working directory, creates deps-uni/lib/*
# Assumes that the user has enough rights to run brew fetch

# Downloading
SED_PGM='s/^Downloaded to: |^Already downloaded: //p'
mkdir deps-amd64
cd deps-amd64
tar xf "`brew fetch --bottle-tag=x86_64_$OSX_VERSION gmp  | grep -F bottle.tar.gz | sed -Ene "$SED_PGM"`"
tar xf "`brew fetch --bottle-tag=x86_64_$OSX_VERSION mpfr | grep -F bottle.tar.gz | sed -Ene "$SED_PGM"`"
cd ..
mkdir deps-arm64
cd deps-arm64
tar xf "`brew fetch --bottle-tag=arm64_$OSX_VERSION gmp   | grep -F bottle.tar.gz | sed -Ene "$SED_PGM"`"
tar xf "`brew fetch --bottle-tag=arm64_$OSX_VERSION mpfr  | grep -F bottle.tar.gz | sed -Ene "$SED_PGM"`"
cd ..

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
