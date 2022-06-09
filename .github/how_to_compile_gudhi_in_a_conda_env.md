# Install a conda development environment to compile GUDHI

## Install miniconda

Download the [installer](https://docs.conda.io/en/latest/miniconda.html) required by your system and follow the [instructions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

## Create a dedicated environment

```bash
conda install -c conda-forge mamba  # installation with mamba is faster
conda create --name gudhi
conda activate gudhi
mamba install -c conda-forge python cmake doxygen eigen cgal-cpp
```

Some of the requirements are in the gudhi-devel repository (please refer to
[how to use github to contribute to gudhi](how_to_use_github_to_contribute_to_gudhi.md)).
Once the gudhi-devel repository is cloned on your machine (`git clone...`) - let's call it `/workdir/gudhi-devel` i.e. -
and once the submodules are initialised (`git submodule update --init`):

```bash
pip install -r ext/gudhi-deploy/build-requirements.txt 
pip install -r ext/gudhi-deploy/test-requirements.txt  # pytorch can be painful to install - not mandatory
```

## Compilation

In order to compile all c++ utilities, examples, benchmarks, unitary tests, and python module:
```bash
cd /workdir/gudhi-devel
rm -rf build; mkdir build  # /!\ any existing build folder will be removed
cd build
# To build all even examples and benchmarks
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=$CONDA_PREFIX -DWITH_GUDHI_EXAMPLE=ON -DWITH_GUDHI_BENCHMARK=ON ..
```

### Specific python compilation

In order to compile only python module
```bash
cd /workdir/gudhi-devel
rm -rf build; mkdir build  # /!\ any existing build folder will be removed
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=$CONDA_PREFIX ..
cd src/python
# To build python module in parallel
python setup.py build_ext -j 16 --inplace  # 16 is the number of CPU that are used to compile the python module. Can be any other value.
# to clean the build
# python setup.py clean --all
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
make user_version;
cd version
mkdir build
cd build
# python OFF to prevent python modules search makes cmake faster
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=$CONDA_PREFIX -DWITH_GUDHI_PYTHON=OFF  ..
make doxygen 2>&1 | tee dox.log
grep warning dox.log # Warnings can be lost with parallel doxygen
firefox html/index.html # [optional] To display the c++ documentation. Anything else than firefox can be used.
```

### Specific python documentation generation

```bash
cd /workdir/gudhi-devel
rm -rf build; mkdir build  # /!\ any existing build folder will be removed
cd build
# python OFF to prevent python modules search makes cmake faster - it is the next cmake call in user version that matters
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=$CONDA_PREFIX -DWITH_GUDHI_PYTHON=OFF -DUSER_VERSION_DIR=version ..
make user_version;
cd version
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=$CONDA_PREFIX ..
cd python
# To build python module in parallel
python setup.py build_ext -j 16 --inplace  # 16 is the number of CPU that are used to compile the python module. Can be any other value.
firefox sphinx/index.html # [optional] To display the python documentation. Anything else than firefox can be used.
```