# Tests strategy

This document tries to sum up the tests strategy that has been put in place for gudhi continuous integration.

The aim is to help maintainers to anticipate third parties modifications, updates.

## Builds

### Linux

As all the third parties are already installed (thanks to docker), the compilations has been seperated by categories to be parallelized:

* examples (C++)
* tests (C++)
* utils (C++)
* doxygen (C++ documentation that is available in the artefacts)
* python (including documentation and code coverage that are available in the artefacts)

(cf. `.circleci/config.yml`)

These build categories are done with and without CGAL, and, with and without Eigen to be sure the users won't be annoyed if a third party is missing.

With CGAL and with Eigen builds are performed inside the docker image `gudhi/ci_for_gudhi` based on `Dockerfile_for_circleci_image` file.
Without CGAL, and, with or without Eigen builds are performed inside the docker image `gudhi/ci_for_gudhi_wo_cgal` based on `Dockerfile_for_circleci_image_without_cgal` file.

#### Update docker images

C++ third parties installation are done thanks to apt on Ubuntu latest LTS.

Docker images need to be rebuild and push each time `.github/build-requirements`, `.github/test-requirements`, when a new third party is added, when a new CGAL version improves gudhi performances, ...

```bash
docker build -f Dockerfile_for_circleci_image -t gudhi/ci_for_gudhi:latest .
docker build -f Dockerfile_for_circleci_image_without_cgal -t gudhi/ci_for_gudhi_wo_cgal:latest .
docker login # requires some specific rights on https://hub.docker.com/u/gudhi/repository/docker/gudhi
docker push gudhi/ci_for_gudhi:latest
docker push gudhi/ci_for_gudhi_wo_cgal:latest
```

### Windows

The compilations has been seperated by categories to be parallelized, but I don't know why builds are not run in parallel:

* examples (C++)
* tests (C++)
* utils (C++)
* python

Doxygen (C++) is not tested.
(cf. `.appveyor.yml`)

C++ third parties installation are done thanks to [vcpkg](https://github.com/microsoft/vcpkg/).
In case of installation issue, check in [vcpkg issues](https://github.com/microsoft/vcpkg/issues).

### OSx

The compilations has been seperated by categories to be parallelized:

* examples (C++)
* tests (C++)
* utils (C++)
* python
* Doxygen (C++)

(cf. `azure-pipelines.yml`)

C++ third parties installation are done thanks to [brew](https://formulae.brew.sh/formula/).
In case of installation issue, check in formula issues.

## Pip packaging

Pip packaging is done in 2 parts:

* on push and pull requests, the wheels are built (pip package dry-run)
* on releases, the wheels are built and sent to pypi.org (package)

Only the Linux pip package is based on a docker image (`gudhi/pip_for_gudhi` based on `Dockerfile_for_pip` file) to make it faster.

### Update docker image

C++ third parties installation are done thanks to yum on an image based on `quay.io/pypa/manylinux2014_x86_64`.

Docker image need to be rebuild and push each time `.github/build-requirements`, when a new third party is added, when a new CGAL version improves gudhi performances, ...
As `.github/test-requirements` is not installed, no need to rebuild image when this file is modified.

```bash
docker build -f Dockerfile_for_pip -t gudhi/pip_for_gudhi:latest .
docker login # requires some specific rights on https://hub.docker.com/u/gudhi/repository/docker/gudhi
docker push gudhi/pip_for_gudhi:latest
```
