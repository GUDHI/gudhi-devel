# Create a new GUDHI version

We will consider that all operations will be performed in a brand new clone of the main project:
```bash
git clone https://github.com/GUDHI/gudhi-devel.git
cd gudhi-devel
```

## Release note

Verify that [.github/next_release.md](.github/next_release.md) is up-to-date with the development performed on the release.
Also check [.github/changelog.md](.github/changelog.md) is appended with the release note.

## Python documentation version switcher

Verify that [.github/for_maintainers/switcher.json](.github/for_maintainers/switcher.json) contains the new version (not necessary for release candidates versions).

## Unlabelled issues and pull requests

Check there are no [unlabelled closed issues](https://github.com/GUDHI/gudhi-devel/issues?q=is%3Aissue+is%3Aclosed+no%3Alabel)
nor [unlabelled closed PRs](https://github.com/GUDHI/gudhi-devel/pulls?q=is%3Apr+is%3Aclosed+no%3Alabel) related to the release.

## Version file modification

**Edit the file CMakeGUDHIVersion.txt**, and increment major, minor, or patch version number, in function of the version new delivery.
```bash
# cf. .gitignore - ignore this if it is a fresh clone version
rm -rf data/points/COIL_database/lucky_cat.off_dist data/points/COIL_database/lucky_cat.off_sc.dot data/points/KleinBottle5D.off_dist data/points/KleinBottle5D.off_sc.dot data/points/human.off_dist data/points/human.off_sc.off data/points/human.off_sc.txt src/python/test/__pycache__ src/python/gudhi/*.cpp
```

Checkin the modifications, build and test the version:
```bash
git submodule update --init
rm -rf build; mkdir build; cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=$CONDA_PREFIX -DWITH_GUDHI_REMOTE_TEST=ON -DWITH_GUDHI_EXAMPLE=ON -DWITH_GUDHI_BENCHMARK=ON  -DUSER_VERSION_DIR=gudhi.@GUDHI_VERSION@ ..
make user_version
date +"%d-%m-%Y-%T" > gudhi.@GUDHI_VERSION@/timestamp.txt
tar -czvf gudhi.@GUDHI_VERSION@.tar.gz gudhi.@GUDHI_VERSION@
sha256sum gudhi.@GUDHI_VERSION@.tar.gz > sha256sum.txt
sha512sum gudhi.@GUDHI_VERSION@.tar.gz > sha512sum.txt

make && ctest --output-on-failure
```

***[Check there are no error]***

## Create the documentation
```bash
mkdir gudhi.doc.@GUDHI_VERSION@
```

***[Check there are no error and the warnings]***

```bash
cd gudhi.@GUDHI_VERSION@
rm -rf build; mkdir build; cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=$CONDA_PREFIX -DWITH_GUDHI_EXAMPLE=ON ..
make doxygen && grep warning doxygen.log
```

***[Check there are no error and the warnings]***

```bash
cp -R html ../../gudhi.doc.@GUDHI_VERSION@/cpp
export LC_ALL=en_US.UTF-8  # cf. bug https://github.com/GUDHI/gudhi-devel/issues/111
make sphinx
```

***[Check there are no error]***

```bash
cp -R python/sphinx ../../gudhi.doc.@GUDHI_VERSION@/python
cd ../..
tar -czvf gudhi.doc.@GUDHI_VERSION@.tar.gz gudhi.doc.@GUDHI_VERSION@

cd gudhi.@GUDHI_VERSION@/build
make && ctest --output-on-failure
```

***[Check there are no error]***

## Upload the documentation

[GUDHI GitHub pages](https://gudhi.github.io/) is only used as a _"qualification"_ web hosting service.
The _"production"_ web hosting service is https://files.inria.fr (cf. [this doc](https://doc-si.inria.fr/display/SU/Espace+web)
or [this one](https://www.nextinpact.com/article/30325/109058-se-connecter-a-serveur-webdav-sous-linux-macos-ou-windows)).

Upload the content of the directory gudhi.doc.@GUDHI_VERSION@/cpp in a new directory on gudhi WebDAV in doc/@GUDHI_VERSION@
Delete the directory doc/latest on gudhi WebDAV.
Copy gudhi WebDAV doc/@GUDHI_VERSION@ as doc/latest (no symbolic link with WebDAV).

Upload the content of the directory gudhi.doc.@GUDHI_VERSION@/python in a new directory on gudhi WebDAV in python/@GUDHI_VERSION@
Delete the directory python/latest on gudhi WebDAV.
Copy gudhi WebDAV python/@GUDHI_VERSION@ as python/latest (no symbolic link with WebDAV).

Upload [.github/for_maintainers/switcher.json](.github/for_maintainers/switcher.json) at root of gudhi WebDAV server.

## Put a version label on files

* Go on page https://github.com/GUDHI/gudhi-devel/releases/new
* Name the tag: tags/gudhi-release-@GUDHI_VERSION@
* Name the release GUDHI @GUDHI_VERSION@ release
* Write the release note
* Drag'n drop *gudhi.@GUDHI_VERSION@.tar.gz*, *sha256sum.txt*, *sha512sum.txt* files
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

After pushing the changes in gudhi-deploy repo main branch, the docker image needs to be built and pushed at
[latest_gudhi_version](https://hub.docker.com/repository/docker/gudhi/latest_gudhi_version)
docker image on docker hub:

```bash
docker build -t gudhi/latest_gudhi_version:3.X.X -f Dockerfile_for_gudhi_installation .
docker tag gudhi/latest_gudhi_version:3.X.X gudhi/latest_gudhi_version:latest
docker login
docker push gudhi/latest_gudhi_version:3.X.X gudhi/latest_gudhi_version:latest
```

***[Check there are no error]***

## Mail sending
Send version mail to the following lists :
* gudhi-devel@inria.fr
* gudhi-users@inria.fr (not for release candidate)

## File to update after the release

**Edit the file CMakeGUDHIVersion.txt**, and increment major, minor, or patch version number, in function of the future version (something like `3.X+1.0a0`).

Reset [.github/next_release.md](.github/next_release.md) with [.github/for_maintainers/next_release_template.md](.github/for_maintainers/next_release_template.md).