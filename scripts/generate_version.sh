#!/bin/bash
#usage bash generate_version.sh : dont generate if svn st non empty
#usage bash generate_version.sh -f : generate even if svn is empty
#usage bash generate_version.sh -f DIR : generate even if svn is empty and save library in dir
#
# 23/06/2015 - Remove source, add biblio, and doc
# 06/10/2015 - Replace static Version.txt with generated GUDHIVersion.cmake file
# VERSION CHECK
ROOT_DIR=..
VERSION_FILE="$ROOT_DIR/GUDHIVersion.cmake"

if [ ! -f $VERSION_FILE ]; then
    echo "File not found! : $VERSION_FILE - Please launch cmake first to generate file"
    exit 1
fi

# SVN STATUS CHECK OR IF FORCED BY USER
if [ "$1" != "-f" ]
then
  SVN_STATUS=`svn status $ROOT_DIR | grep -v $VERSION_FILE`
  echo $SVN_STATUS
fi


TARGET_DIR=""
if [ "$2" != "-f" ]
then
  TARGET_DIR=$2
  echo "Install folder : $TARGET_DIR"
fi

if [ "$SVN_STATUS" != "" ]
then
  echo "Svn status must be empty to create a version!"
  exit 1
fi

# TEMPORARY FOLDER CREATION
VERSION_DATE=`date +"%Y-%m-%d-%H-%M-%S"`
GUDHI="_GUDHI_"
VERSION_REVISION=`cat $VERSION_FILE`
VERSION_DIR="$VERSION_DATE$GUDHI$VERSION_REVISION"
echo $VERSION_DIR
mkdir "$VERSION_DIR"

# TOP LEVEL FILE COPY
cp $ROOT_DIR/README $VERSION_DIR
cp $ROOT_DIR/Conventions.txt $VERSION_DIR
cp $ROOT_DIR/COPYING $VERSION_DIR
cp -R $ROOT_DIR/data $VERSION_DIR
cp $ROOT_DIR/src/CMakeLists.txt $VERSION_DIR
cp $ROOT_DIR/src/Doxyfile $VERSION_DIR
cp -R $ROOT_DIR/biblio $VERSION_DIR
cp $ROOT_DIR/src/GUDHIConfigVersion.cmake.in $VERSION_DIR
cp $ROOT_DIR/src/GUDHIConfig.cmake.in $VERSION_DIR
cp $ROOT_DIR/CMakeGUDHIVersion.txt $VERSION_DIR
cp $ROOT_DIR/GUDHIVersion.cmake.in $VERSION_DIR

# PACKAGE LEVEL COPY
PACKAGE_INC_DIR="/include"
PACKAGE_EX_DIR="/example"
PACKAGE_CONCEPT_DIR="/concept"
PACKAGE_DOC_DIR="/doc"
for package in `ls $ROOT_DIR/src/`
do
  if [ -d "$ROOT_DIR/src/$package" ] && [ $package != "Bottleneck" ]
  then
    echo $package
    if [ "$package" == "cmake" ] || [ "$package" == "debian" ]
    then
      # SPECIFIC FOR CMAKE MODULES
      cp -R $ROOT_DIR/src/$package $VERSION_DIR
    elif [ "$package" == "GudhUI" ]
    then
      # SPECIFIC FOR GUDHI USER INTERFACE
      cp -R $ROOT_DIR/src/$package $VERSION_DIR
    elif [ "$package" == "cython" ]
    then
      # SPECIFIC FOR CYTHON INTERFACE
      cp -R $ROOT_DIR/src/$package $VERSION_DIR
    else
      # PACKAGE COPY
      if [ -d "$ROOT_DIR/src/$package$PACKAGE_INC_DIR" ]
      then
        if [ ! -d "$VERSION_DIR$PACKAGE_INC_DIR" ]
        then
          # MUST CREATE DIRECTORY ON FIRST LOOP
          mkdir $VERSION_DIR$PACKAGE_INC_DIR
        fi
        cp -R $ROOT_DIR/src/$package$PACKAGE_INC_DIR/* $VERSION_DIR$PACKAGE_INC_DIR/
      fi
      if [ -d "$ROOT_DIR/src/$package$PACKAGE_EX_DIR" ]
      then
        mkdir -p $VERSION_DIR$PACKAGE_EX_DIR/$package
        cp -R $ROOT_DIR/src/$package$PACKAGE_EX_DIR/* $VERSION_DIR$PACKAGE_EX_DIR/$package
      fi
      if [ -d "$ROOT_DIR/src/$package$PACKAGE_CONCEPT_DIR" ]
      then
        mkdir -p $VERSION_DIR$PACKAGE_CONCEPT_DIR/$package
        cp -R $ROOT_DIR/src/$package$PACKAGE_CONCEPT_DIR/* $VERSION_DIR$PACKAGE_CONCEPT_DIR/$package
      fi
      if [ -d "$ROOT_DIR/src/$package$PACKAGE_DOC_DIR" ]
      then
        mkdir -p $VERSION_DIR$PACKAGE_DOC_DIR/$package
        cp -R $ROOT_DIR/src/$package$PACKAGE_DOC_DIR/* $VERSION_DIR$PACKAGE_DOC_DIR/$package
      fi
    fi
  fi
done


#INSTALL to some directory 
if [ "$TARGET_DIR" != "" ]; then
  echo "Install in dir $TARGET_DIR"	
  mv "$VERSION_DIR" "$TARGET_DIR"
else
  # ZIP DIR AND REMOVE IT
  tar -zcf "$VERSION_DIR.tar.gz" "$VERSION_DIR"
  rm -rf "$VERSION_DIR"
fi




