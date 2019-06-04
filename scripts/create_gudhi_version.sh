#!/bin/bash

login="vrouvrea"
version="2.3.0"
cgaldir="/home/vincent/workspace/CGAL-4.11-HO/build"
cpucount=7


# We start from scripts dir in the dev branch
cd ..
RELATIVEURL=`svn info . |grep -F "Relative URL:" | awk '{print $NF}'`

if [ "$RELATIVEURL" != "^/trunk" ]
then
echo "Script must be launched in trunk and not in $RELATIVEURL"
exit
fi

rm -rf build; mkdir build; cd build; cmake  -DCMAKE_BUILD_TYPE=Debug -DDEBUG_TRACES=ON -DCGAL_DIR=${cgaldir} -DWITH_GUDHI_EXAMPLE=ON -DWITH_GUDHI_BENCHMARK=ON -DPython_ADDITIONAL_VERSIONS=3 ..
cmake  -DCMAKE_BUILD_TYPE=Debug .

CURRENTDIRECTORY=`pwd`
export PYTHONPATH=$CURRENTDIRECTORY/src/cython:$PYTHONPATH

make -j ${cpucount} all test

cd ..
svn st | grep -v GUDHIVersion.cmake | grep "^\?" | awk "{print \$2}" | xargs rm -rf

svn copy svn+ssh://${login}@scm.gforge.inria.fr/svnroot/gudhi/trunk svn+ssh://${login}@scm.gforge.inria.fr/svnroot/gudhi/tags/gudhi-release-${version} \
  -m "Creating a tag of Gudhi release version ${version}."

cd build
make user_version

userversiondir=`find . -type d -name "*_GUDHI_${version}" | sed 's/\.\///g'`
echo "User version directory = ${userversiondir}"

tar -czvf ${userversiondir}.tar.gz ${userversiondir}

userdocdir=${userversiondir/GUDHI/GUDHI_DOC}
echo "User documentation directory = ${userdocdir}"
mkdir ${userdocdir}
make doxygen

cp -R ${userversiondir}/doc/html ${userdocdir}/cpp
cd ${userversiondir}
rm -rf build; mkdir build; cd build; cmake  -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=./installed -DCGAL_DIR=${cgaldir} -DWITH_GUDHI_EXAMPLE=ON -DPython_ADDITIONAL_VERSIONS=3 ..

CURRENTDIRECTORY=`pwd`
export PYTHONPATH=$CURRENTDIRECTORY/cython:$PYTHONPATH

make sphinx

cp -R cython/sphinx ../../${userdocdir}/python
cd ../..
tar -czvf ${userdocdir}.tar.gz ${userdocdir}

cd ${userversiondir}/build
make -j ${cpucount} all test install

cd ../..
actualdir=`pwd`
echo "Library is available at ${actualdir}/${userversiondir}.tar.gz"
sha256sum ${userversiondir}.tar.gz
echo "Documentation is available at ${actualdir}/${userdocdir}.tar.gz"
