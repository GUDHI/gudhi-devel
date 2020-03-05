# Create a new GUDHI version

We will consider that all operations will be performed in a brand new clone of the main project:
```bash
git clone https://github.com/GUDHI/gudhi-devel.git
cd gudhi-devel
```

## Version file modification

**Edit the file CMakeGUDHIVersion.txt**, and increment major, minor, or patch version number, in function of the version new delivery.
```bash
# cf. .gitignore - ignore this if it is a fresh clone version
rm -rf data/points/COIL_database/lucky_cat.off_dist data/points/COIL_database/lucky_cat.off_sc.dot data/points/KleinBottle5D.off_dist data/points/KleinBottle5D.off_sc.dot data/points/human.off_dist data/points/human.off_sc.off data/points/human.off_sc.txt
```

Checkin the modifications, build and test the version:
```bash
mkdir build
cd build
cmake -DCGAL_DIR=/your/path/to/CGAL -DWITH_GUDHI_EXAMPLE=ON -DWITH_GUDHI_BENCHMARK=ON  -DUSER_VERSION_DIR=gudhi.@GUDHI_VERSION@ -DPython_ADDITIONAL_VERSIONS=3 ..
make user_version
date +"%d-%m-%Y-%T" > gudhi.@GUDHI_VERSION@/timestamp.txt
tar -czvf gudhi.@GUDHI_VERSION@.tar.gz gudhi.@GUDHI_VERSION@
md5sum gudhi.@GUDHI_VERSION@.tar.gz > md5sum.txt
sha256sum gudhi.@GUDHI_VERSION@.tar.gz > sha256sum.txt
sha512sum gudhi.@GUDHI_VERSION@.tar.gz > sha512sum.txt

make -j all test
```

***[Check there are no error]***

## Create the documentation
```bash
mkdir gudhi.doc.@GUDHI_VERSION@
make doxygen  2>&1 | tee dox.log && grep warning dox.log
```

***[Check there are no error and the warnings]***

```bash
cp -R gudhi.@GUDHI_VERSION@/doc/html gudhi.doc.@GUDHI_VERSION@/cpp
cd gudhi.@GUDHI_VERSION@
rm -rf build; mkdir build; cd build; cmake -DCGAL_DIR=/your/path/to/CGAL -DWITH_GUDHI_EXAMPLE=ON -DPython_ADDITIONAL_VERSIONS=3 ..
export LC_ALL=en_US.UTF-8  # cf. bug
make sphinx
```

***[Check there are no error]***

```bash
cp -R python/sphinx ../../gudhi.doc.@GUDHI_VERSION@/python
cd ../..
tar -czvf gudhi.doc.@GUDHI_VERSION@.tar.gz gudhi.doc.@GUDHI_VERSION@

cd gudhi.@GUDHI_VERSION@/build
make all test
```

***[Check there are no error]***

## Upload the documentation

Upload by ftp the content of the directory gudhi.doc.@GUDHI_VERSION@/cpp in a new directory on ForgeLogin@scm.gforge.inria.fr:/home/groups/gudhi/htdocs/doc/@GUDHI_VERSION@

Upload by ftp the content of the directory gudhi.doc.@GUDHI_VERSION@/python in a new directory on ForgeLogin@scm.gforge.inria.fr:/home/groups/gudhi/htdocs/python/@GUDHI_VERSION@

Through ssh, make the **latest** link to your new version of the documentation:
```bash
ssh ForgeLogin@scm.gforge.inria.fr
cd /home/groups/gudhi/htdocs/doc
rm latest
ln -s @GUDHI_VERSION@ latest
cd /home/groups/gudhi/htdocs/python
rm latest
ln -s @GUDHI_VERSION@ latest
```

## Put a version label on files

* Go on page https://github.com/GUDHI/gudhi-devel/releases/new
* Name the tag: tags/gudhi-release-@GUDHI_VERSION@
* Name the release GUDHI @GUDHI_VERSION@
* Write the release note
* Drag'n drop *gudhi.@GUDHI_VERSION@.tar.gz*, *md5sum.txt*, *sha256sum.txt*, *sha512sum.txt* files
* Tick the *This is a pre-release* check button if this is a release candidate (untick if this is an official version)
* Click the *Publish the release* button

## Mail sending
Send version mail to the following lists :
* gudhi-devel@lists.gforge.inria.fr
* gudhi-users@lists.gforge.inria.fr (not for release candidate)

