# Install a conda development environment to compile GUDHI

## Install a conda distribution

Install one of the [conda distribution](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).
Miniforge is the preferred version, [here is its installation guide](https://conda-forge.org/download/)

## Create a dedicated environment

If you went with Miniforge, `conda install -c conda-forge mamba` can be skipped as it is already installed,
and you can skip all `-c conda-forge` as it is the default channel with Miniforge.

```bash
conda install -c conda-forge mamba  # installation with mamba is faster
conda create --name gudhi
conda activate gudhi
mamba install -c conda-forge python cmake eigen cgal-cpp
```

If you do not have any C++ compiler installed - it is required - on your machine (g++, clang++, Visual Studio, ...), you can:
```bash
mamba install -c conda-forge cxx-compiler
```

If you want to modify the C++ documentation - it is optional - you will have to install doxygen to test your modifications:
```bash
mamba install -c conda-forge doxygen
```

Some of the requirements are in the gudhi-devel repository (please refer to
[how to use github to contribute to gudhi](how_to_use_github_to_contribute_to_gudhi.md)).
Once the gudhi-devel repository is cloned on your machine (`git clone...`) - let's call it `/workdir/gudhi-devel` i.e. -
and once the submodules are initialised (`git submodule update --init`):

To install mandatory packages to compile the GUDHI python module:
```bash
pip install -r ext/gudhi-deploy/build-requirements.txt 
```

If you want to modify the Python code or documentation and want to test it - it is optional - you will have to install:
```bash
pip install -r ext/gudhi-deploy/test-requirements.txt  # pytorch can be painful to install - not mandatory
```

## Compilation

It is not mandatory, and it can be quite long to compile everything ([ccache](https://ccache.dev/) can help here),
but in order to compile all c++ utilities, examples, benchmarks, unitary tests, and python module:
```bash
cd /workdir/gudhi-devel
rm -rf build; mkdir build  # /!\ any existing build folder will be removed
cd build
# To build all even examples and benchmarks
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=$CONDA_PREFIX -DWITH_GUDHI_EXAMPLE=ON -DWITH_GUDHI_BENCHMARK=ON ..
```

### Specific python compilation

In order to compile only the python module, it is the same process, but just change to `src/python` directory:
```bash
cd /workdir/gudhi-devel
rm -rf build; mkdir build  # /!\ any existing build folder will be removed
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=$CONDA_PREFIX ..
cd src/python
make -j 16  # 16 is the number of CPU that are used to compile the python module. Can be any other value.
# to clean the build
# make clean
```

In order to use freshly compiled gudhi python module:
```bash
PYTHONPATH=/workdir/gudhi-devel/build/src/python python # or ipython, jupyter, ...
```

### Specific C++ documentation generation

```bash
cd /workdir/gudhi-devel
rm -rf build; mkdir build  # /!\ any existing build folder will be removed
cd build
# python OFF to prevent python modules search makes cmake faster
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=$CONDA_PREFIX -DWITH_GUDHI_PYTHON=OFF -DUSER_VERSION_DIR=version ..
# C++ documentation generation must be performed inside the "user" version of the library
make user_version
cd version
mkdir build
cd build
# python OFF to prevent python modules search makes cmake faster
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=$CONDA_PREFIX -DWITH_GUDHI_PYTHON=OFF  ..
make doxygen
grep warning doxygen.log # Warnings can be lost with parallel doxygen
firefox html/index.html # [optional] To display the c++ documentation. Anything else than firefox can be used.
```

### Specific python documentation generation

```bash
cd /workdir/gudhi-devel
rm -rf build; mkdir build  # /!\ any existing build folder will be removed
cd build
# python OFF to prevent python modules search makes cmake faster - it is the next cmake call in user version that matters
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=$CONDA_PREFIX -DWITH_GUDHI_PYTHON=OFF -DUSER_VERSION_DIR=version ..
# Python documentation generation must be performed inside the "user" version of the library
make user_version
cd version
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=$CONDA_PREFIX ..
cd python
# To build python module in parallel
make sphinx
firefox sphinx/index.html # [optional] To display the python documentation. Anything else than firefox can be used.
```
