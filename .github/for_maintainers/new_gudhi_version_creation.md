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
git submodule update --init
rm -rf build; mkdir build; cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCGAL_DIR=/your/path/to/CGAL -DWITH_GUDHI_EXAMPLE=ON -DWITH_GUDHI_BENCHMARK=ON  -DUSER_VERSION_DIR=gudhi.@GUDHI_VERSION@ -DPython_ADDITIONAL_VERSIONS=3 ..
make user_version
date +"%d-%m-%Y-%T" > gudhi.@GUDHI_VERSION@/timestamp.txt
tar -czvf gudhi.@GUDHI_VERSION@.tar.gz gudhi.@GUDHI_VERSION@
md5sum gudhi.@GUDHI_VERSION@.tar.gz > md5sum.txt
sha256sum gudhi.@GUDHI_VERSION@.tar.gz > sha256sum.txt
sha512sum gudhi.@GUDHI_VERSION@.tar.gz > sha512sum.txt

make -j 4 all && ctest -j 4 --output-on-failure
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
rm -rf build; mkdir build; cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCGAL_DIR=/your/path/to/CGAL -DWITH_GUDHI_EXAMPLE=ON -DPython_ADDITIONAL_VERSIONS=3 ..
export LC_ALL=en_US.UTF-8  # cf. bug https://github.com/GUDHI/gudhi-devel/issues/111
make sphinx
```

***[Check there are no error]***

```bash
cp -R python/sphinx ../../gudhi.doc.@GUDHI_VERSION@/python
cd ../..
tar -czvf gudhi.doc.@GUDHI_VERSION@.tar.gz gudhi.doc.@GUDHI_VERSION@

cd gudhi.@GUDHI_VERSION@/build
make -j 4 all && ctest -j 4 --output-on-failure
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
* Name the release GUDHI @GUDHI_VERSION@ release
* Write the release note
* Drag'n drop *gudhi.@GUDHI_VERSION@.tar.gz*, *md5sum.txt*, *sha256sum.txt*, *sha512sum.txt* files
* Tick the *This is a pre-release* check button if this is a release candidate (untick if this is an official version)
* Click the *Publish the release* button

## Pip package

The pip package construction shall be started on release creation, you just have to check
[gudhi github actions](https://github.com/GUDHI/gudhi-devel/actions) results.
The version number must be conform to [pep440](https://www.python.org/dev/peps/pep-0440/#pre-releases)

## Conda package

You have to fork [conda-forge/gudhi-feedstock](https://github.com/conda-forge/gudhi-feedstock).
The main changes consist into changing in the `recipe/meta.yaml`:
* `{% set version = "@GUDHI_VERSION@" %}`
* The cgal-cpp version number with the last one (you can find it [here](https://anaconda.org/conda-forge/cgal-cpp)) in the `host:` and the `run:` sections

Create a Pull Request (PR) from this fork.
If you need to update conda tools (conda-build, conda-smithy, ...), add a comment in your PR saying `@conda-forge-admin, please rerender`, it will done automatically (do not forget to `git pull` the changes).

## Docker image

You have to modify the
[Dockerfile_gudhi_installation](https://github.com/GUDHI/gudhi-deploy/blob/main/Dockerfile_for_gudhi_installation)
in gudhi-deploy repository in order to use the last release, cf. lines:
```
...
ARG GUDHI_VERSION="3.X.X"
...
```

After pushing the changes the docker image build will be automatically performed for
[latest_gudhi_version](https://hub.docker.com/repository/docker/gudhi/latest_gudhi_version)
docker image on docker hub.

***[Check there are no error]***

## Mail sending
Send version mail to the following lists :
* gudhi-devel@lists.gforge.inria.fr
* gudhi-users@lists.gforge.inria.fr (not for release candidate)
