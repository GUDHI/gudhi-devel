# Naming conventions

## C++

### In the code:
* The classes and functions of a package should be in a sub-namespace of the `Gudhi` namespace. The sub-namespace names are in lowercase and use underscore separators. E.g. `Gudhi::package_name::`
* Concepts are named with camel case starting with uppercase. E.g. `PersistentHomology` for the concept of Persitence homology.
* Classes start with an uppercase letter and use underscore separators. E.g. `Skeleton_blocker_contractor`.
* Member functions and free functions are in lowercase and use underscore separators. E.g. `int num_vertices()`.
* Constants and macros are in uppercase.
* Macros should begin with the prefix `GUDHI_`.

### File names:
* All headers are named *.h and all sources are named *.cpp.
* If a single class or function is provided in a file, its name (with the same letter case) should be used for the file name.
* If a file does not contain a single class, its name should not begin with a capital letter.
* Test files should be called `test_[what_is_tested].cpp`. E.g. `test_sparsify_point_set.cpp`
* Example files should be called `example_[what_it_is].cpp`. E.g. `example_sparsify_point_set.cpp`

### In CMakeLists.txt files:
* The name of the "project" should be in this form: `Package_[tests|examples|…]`. E.g. `project(Simplex_tree_examples)`.
* The name if each "target" (first parameter of add_executable) should be in this form: `Package_{name of the cpp file without extension}`. E.g `add_executable(Subsampling_test_sparsify_point_set test_sparsify_point_set.cpp)`.

### Code style
We are using [google c++ style guide](https://google.github.io/styleguide/cppguide.html) recommendations with 120 characters per line of code.
[clang-format](https://clang.llvm.org/docs/ClangFormat.html) can be used to format automatically your code:
```bash
cd src # there is a .clang-format file with these specifications
clang-format -style=file -i Simplex_tree/include/gudhi/Simplex_tree.h # -i means in place, your file will be modified
```

## Python

In progress...

### Code style
We are using [PEP8 Python style guide](https://www.python.org/dev/peps/pep-0008/) recommendations with 120 characters per line of code.
[black](https://black.readthedocs.io/en/stable/) can be used to format automatically your code:
```bash
black -l 120 src/python/example/bottleneck_basic_example.py
```
