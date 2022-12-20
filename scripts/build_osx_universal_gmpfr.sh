#!/bin/bash
set -e

# In the working directory, creates deps-uni/lib/*
# Assumes that the user has enough rights to run brew fetch

# Downloading
mkdir deps-amd64
cd deps-amd64
tar xf "`brew fetch --bottle-tag=big_sur gmp        | sed -ne 's/^Downloaded to: //p'`"
tar xf "`brew fetch --bottle-tag=big_sur mpfr       | sed -ne 's/^Downloaded to: //p'`"
cd ..
mkdir deps-arm64
cd deps-arm64
tar xf "`brew fetch --bottle-tag=arm64_big_sur gmp  | sed -ne 's/^Downloaded to: //p'`"
tar xf "`brew fetch --bottle-tag=arm64_big_sur mpfr | sed -ne 's/^Downloaded to: //p'`"
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
BADGMP=`otool -L deps-uni/lib/$MPFR|sed -ne 's/[[:space:]]*\(.*libgmp\..*dylib\).*/\1/p'`
install_name_tool -change $BADGMP $PWD/deps-uni/lib/$GMP deps-uni/lib/$MPFR
BADGMP=`otool -L deps-uni/lib/$GMPXX|sed -ne 's/[[:space:]]*\(.*libgmp\..*dylib\).*/\1/p'`
install_name_tool -change $BADGMP $PWD/deps-uni/lib/$GMP deps-uni/lib/$GMPXX

ln -s $GMP deps-uni/lib/libgmp.dylib
ln -s $GMPXX deps-uni/lib/libgmpxx.dylib
ln -s $MPFR deps-uni/lib/libmpfr.dylib

# Debug
ls -l deps-uni/lib
otool -L deps-uni/lib/*.*.dylib
