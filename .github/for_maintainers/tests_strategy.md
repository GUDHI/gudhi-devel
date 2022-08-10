# Tests strategy

This document tries to sum up the tests strategy that has been put in place for gudhi continuous integration.

The aim is to help maintainers to anticipate third parties modifications, updates.

## CMake options

[CMake GUDHI options](../../src/cmake/modules/GUDHI_options.cmake) allows to activate/deactivate what should be built and tested.
Note the special option `WITH_GUDHI_THIRD_PARTY` that, when set to `OFF`, accelerates doxygen documentation generation or `user_version` for instance.

## Builds

### Linux

As all the third parties are already installed (thanks to docker), the compilations have been separated in categories to be parallelized:

* examples (C++)
* tests (C++)
* utils (C++)
* doxygen (C++ documentation that is available in the artefacts)
* python (including documentation and code coverage that are available in the artefacts; here the WITH_GUDHI_REMOTE_TEST option is enabled which adds datasets fetching test)

(cf. `.circleci/config.yml`)

These build categories are done with and without CGAL, and, with and without Eigen to be sure the users won't be annoyed if a third party is missing.

With CGAL and with Eigen builds are performed inside the docker image `gudhi/ci_for_gudhi` based on `Dockerfile_for_circleci_image` file.
Without CGAL, and, with or without Eigen builds are performed inside the docker image `gudhi/ci_for_gudhi_wo_cgal` based on `Dockerfile_for_circleci_image_without_cgal` file.

#### Update docker images

C++ third parties installation is done thanks to apt on Ubuntu latest LTS.

Docker images need to be rebuilt and pushed each time `.github/build-requirements`, `.github/test-requirements`, when a new third party is added, when a new CGAL version improves gudhi performances, ...

```bash
docker build -f Dockerfile_for_circleci_image -t gudhi/ci_for_gudhi:latest .
docker build -f Dockerfile_for_circleci_image_without_cgal -t gudhi/ci_for_gudhi_wo_cgal:latest .
docker login # requires some specific rights on https://hub.docker.com/u/gudhi/repository/docker/gudhi
docker push gudhi/ci_for_gudhi:latest
docker push gudhi/ci_for_gudhi_wo_cgal:latest
```

### Windows

The compilations are not parallelized, as installation time (about 30 minutes) is too much compared to
build and tests timings (about 30 minutes). Builds and tests include:

* examples (C++)
* tests (C++)
* utils (C++)
* python (here the WITH_GUDHI_REMOTE_TEST option is enabled which adds datasets fetching test)

Doxygen (C++) is not generated.
(cf. `azure-pipelines.yml`)

C++ third parties installation is done thanks to [vcpkg](https://github.com/microsoft/vcpkg/).
In case of an installation issue, check in [vcpkg issues](https://github.com/microsoft/vcpkg/issues).

### OSx

The compilations are not parallelized, but they should, as installation time (about 4 minutes) is
negligible compared to build and tests timings (about 30 minutes). Builds and tests include:

* examples (C++)
* tests (C++)
* utils (C++)
* python (here the WITH_GUDHI_REMOTE_TEST option is enabled which adds datasets fetching test)
* Doxygen (C++)

(cf. `azure-pipelines.yml`)

C++ third parties installation is done thanks to [brew](https://formulae.brew.sh/formula/).
In case of an installation issue, check in formula issues.

## Pip packaging

Pip packaging is done in 2 parts:

* on push and pull requests, the wheels are built (pip package dry-run)
* on releases, the wheels are built and sent to pypi.org (package)

Only the Linux pip package is based on a docker image (`gudhi/pip_for_gudhi` based on `Dockerfile_for_pip` file) to make it faster.

### Update docker image

C++ third parties installation is done thanks to yum on an image based on `quay.io/pypa/manylinux2014_x86_64`.

Docker image needs to be rebuilt and pushed each time `.github/build-requirements`, when a new third party is added, when a new CGAL version improves gudhi performances, ...
As `.github/test-requirements` is not installed, no need to rebuild image when this file is modified.

```bash
docker build -f Dockerfile_for_pip -t gudhi/pip_for_gudhi:latest .
docker login # requires some specific rights on https://hub.docker.com/u/gudhi/repository/docker/gudhi
docker push gudhi/pip_for_gudhi:latest
```
