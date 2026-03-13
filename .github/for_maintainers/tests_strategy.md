# Tests strategy

This document tries to sum up the tests strategy that has been put in place for gudhi continuous integration.

The aim is to help maintainers to anticipate third parties modifications, updates.

## CMake options

[CMake GUDHI options](../../src/cmake/modules/GUDHI_options.cmake) allows to activate/deactivate what should be built and tested.
Note the special option `WITH_GUDHI_THIRD_PARTY` that, when set to `OFF`, accelerates doxygen documentation generation or `user_version` for instance.

## Builds

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

### Linux

As all the third parties are already installed (thanks to docker), the compilations have been separated in categories to be parallelized:

* examples (C++)
* tests (C++)
* utils (C++)
* python (including documentation and code coverage that are available in the artefacts; here the WITH_GUDHI_REMOTE_TEST option is enabled which adds datasets fetching test)

(cf. `.circleci/config.yml`)

These build categories are done with and without CGAL, and, with and without Eigen to be sure the users won't be annoyed if a third party is missing.

With CGAL and with Eigen builds are performed inside the docker image `ghcr.io/gudhi/gudhi-deploy/gudhi/ci_for_gudhi` based on `https://github.com/GUDHI/gudhi-deploy/Dockerfile_for_circleci_image` file.
Without CGAL, and, with or without Eigen builds are performed inside the docker image `ghcr.io/gudhi/gudhi-deploy/gudhi/ci_for_gudhi_wo_cgal` based on `https://github.com/GUDHI/gudhi-deploy/Dockerfile_for_circleci_image_without_cgal` file.
These docker images are based on a Ubuntu latest LTS and C++ third parties installation is done thanks to apt.

cf. [docker images](tests_strategy#docker-images) section for more details

Another tests are done in parallel:
* doxygen (C++ documentation that is available in the artefacts)
* bibliography (tests that the *.bibtex files are compiling with latexmk)

These are performed inside the docker image `ghcr.io/gudhi/gudhi-deploy/gudhi/doxygen_for_gudhi` based on `https://github.com/GUDHI/gudhi-deploy/Dockerfile_for_doxygen_circleci_image` file.

cf. [docker images](tests_strategy#docker-images) section for more details

## Pip packaging

Pip packaging is done in 2 parts:

* on push and pull requests, the wheels are built (pip package dry-run)
* on releases, the wheels are built and sent to pypi.org (package)

The Linux pip package is built in a docker image (`ghcr.io/gudhi/gudhi-deploy/pip_for_gudhi` based on `https://github.com/GUDHI/gudhi-deploy/Dockerfile_for_pip` file) to make it faster.
This docker image is based on `quay.io/pypa/manylinux_2_28_x86_64` and C++ third parties installation is done thanks to yum.

cf. [docker images](tests_strategy#docker-images) section for more details

#### Docker images

The docker images are automatically generated on each releases of `https://github.com/GUDHI/gudhi-deploy/`. These images are stored on the GitHub docker registry:
* `ghcr.io/gudhi/gudhi-deploy/ci_for_gudhi`
* `ghcr.io/gudhi/gudhi-deploy/ci_for_gudhi_wo_cgal`
* `ghcr.io/gudhi/gudhi-deploy/doxygen_for_gudhi`
* `ghcr.io/gudhi/gudhi-deploy/pip_for_gudhi`

Docker images need to be rebuilt and pushed each time `https://github.com/GUDHI/gudhi-deploy/build-requirements.txt`, `https://github.com/GUDHI/gudhi-deploy/test-requirements.txt`, when a new third party is added, when a new CGAL version improves gudhi performances, ...

To update these docker images, make a pull request on `https://github.com/GUDHI/gudhi-deploy/`, and once it is merged, make a new release (in the `YYYY.MM.NB` format where `YYYY` is the year, `MM` is the month and `NB` is an incremental number).
The docker images will be built on GitHub docker registry, triggered by the given tag, and will be available at `ghcr.io/gudhi/gudhi-deploy/[gudhi docker image]:YYYY.MM.NB`.

Once this is done, it is up to the developper to update on `https://github.com/GUDHI/gudhi-devel/` the use of the new tag `YYYY.MM.NB` in the following files:
* `.circleci/config.yml`
* `.github/workflows/pip-build-linux.yml`
* `.github/workflows/sanitize-unitary-tests-linux.yml`
