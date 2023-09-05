:orphan:

.. To get rid of WARNING: document isn't included in any toctree

Installation
############

Packages
********
The easiest way to install the Python version of GUDHI is using pre-built packages.
We recommend `conda <https://gudhi.inria.fr/conda/>`_

.. code-block:: bash

    conda install -c conda-forge gudhi

Gudhi is also available on `PyPI <https://pypi.org/project/gudhi/>`_

.. code-block:: bash

    pip install gudhi

Third party packages are also available, for instance on Debian or Ubuntu

.. code-block:: bash

    apt install python3-gudhi

In all cases, you may still want to install some of the optional `run time dependencies`_.

Compiling
*********
These instructions are for people who want to compile gudhi from source, they are
unnecessary if you installed a binary package of Gudhi as above. They assume that
you have downloaded a `release <https://github.com/GUDHI/gudhi-devel/releases>`_,
with a name like `gudhi.3.X.Y.tar.gz`, then run `tar xf gudhi.3.X.Y.tar.gz`, which
created a directory `gudhi.3.X.Y`, hereinafter referred to as `/path-to-gudhi/`.
If you are instead using a git checkout, beware that the paths are a bit
different, and in particular the `python/` subdirectory is actually `src/python/`
there.

The library uses c++17 and requires `Boost <https://www.boost.org/>`_ :math:`\geq` 1.71.0,
`CMake <https://www.cmake.org/>`_ :math:`\geq` 3.8,
Python :math:`\geq` 3.5, `NumPy <http://numpy.org>`_ :math:`\geq` 1.15.0, `Cython <https://www.cython.org/>`_
:math:`\geq` 0.27 and `pybind11 <https://github.com/pybind/pybind11>`_ to compile the GUDHI Python module.
It is a multi-platform library and compiles on Linux, Mac OSX and Visual Studio 2017 or later.

If you have several Python/python installed, the version 2.X may be used by default, but you can force it by adding
:code:`-DPython_ADDITIONAL_VERSIONS=3` to the cmake command.

GUDHI Python module compilation
===============================

After making sure that the `Compilation dependencies`_ are properly installed,
one can build the GUDHI Python module, by running the following commands in a terminal:

.. code-block:: bash

    cd /path-to-gudhi/
    mkdir build
    cd build/
    cmake -DCMAKE_BUILD_TYPE=Release ..
    cd python
    make

.. note::

    :code:`make python` (or :code:`make` in python directory) is only a
    `CMake custom targets <https://cmake.org/cmake/help/latest/command/add_custom_target.html>`_
    to shortcut :code:`python setup.py build_ext --inplace` command.
    No specific other options  (:code:`-j8` for parallel, or even :code:`make clean`, ...) are
    available.
    But one can use :code:`python setup.py ...` specific options in the python directory:

.. code-block:: bash

    python setup.py clean --all               # Clean former compilation
    python setup.py build_ext -j 8 --inplace  # Build in parallel

GUDHI Python module installation
================================

Once the compilation succeeds, one can add the GUDHI Python module path to the
PYTHONPATH:

.. code-block:: bash

    # For windows, you have to set PYTHONPATH environment variable
    export PYTHONPATH='$PYTHONPATH:/path-to-gudhi/build/python'

Or install it definitely in your Python packages folder:

.. code-block:: bash

    cd /path-to-gudhi/build/python
    python setup.py install # add --user to the command if you do not have the permission
    # Or 'pip install .'

.. note::

    It does not take into account :code:`CMAKE_INSTALL_PREFIX`.
    But one can use
    `alternate location installation <https://docs.python.org/3/install/#alternate-installation>`_.

Test suites
===========

To test your build, `py.test <http://doc.pytest.org>`_ is required. Run the
following `Ctest <https://cmake.org/cmake/help/latest/manual/ctest.1.html>`_
(CMake test driver program) command in a terminal:

.. code-block:: bash

    cd /path-to-gudhi/build/python
    # For windows, you have to set PYTHONPATH environment variable
    export PYTHONPATH='$PYTHONPATH:/path-to-gudhi/build/python'
    ctest

.. note::

    One can use :code:`ctest` specific options in the python directory:

.. code-block:: bash

    # Launch tests in parallel on 8 cores and set failing tests in verbose mode
    ctest -j 8 --output-on-failure

Debugging issues
================

If tests fail, please check your PYTHONPATH and try to :code:`import gudhi`
and check the errors.
The problem can come from a third-party library bad link or installation.

If :code:`import gudhi` succeeds, please have a look to debug information:

.. code-block:: python

    import gudhi as gd
    print(gd.__debug_info__)
    print("+ Installed modules are: " + gd.__available_modules)
    print("+ Missing modules are: " + gd.__missing_modules)

You shall have something like:

.. code-block:: none

    Pybind11 version 2.8.1 
    Python version 3.7.12
    Cython version 0.29.25 
    Numpy version 1.21.4 
    Boost version 1.77.0
    + Installed modules are: off_utils;simplex_tree;rips_complex;cubical_complex;periodic_cubical_complex;
        persistence_graphical_tools;reader_utils;witness_complex;strong_witness_complex;
    + Missing modules are: bottleneck;nerve_gic;subsampling;tangential_complex;alpha_complex;euclidean_witness_complex;
        euclidean_strong_witness_complex;

Here, you can see that the modules that need CGAL are missing, because CGAL is not installed.
:code:`persistence_graphical_tools` is installed, but 
`its functions <https://gudhi.inria.fr/python/latest/persistence_graphical_tools_ref.html>`_ will produce an error as
matplotlib is not available.
Unitary tests cannot be run as pytest is missing.

A complete configuration would be :

.. code-block:: none

    Pybind11 version 2.8.1 
    Python version 3.9.7
    Cython version 0.29.24 
    Pytest version 6.2.5 
    Matplotlib version 3.5.0 
    Numpy version 1.21.4 
    Scipy version 1.7.3 
    Scikit-learn version 1.0.1 
    POT version 0.8.0 
    HNSWlib found
    PyKeOps version [pyKeOps]: 2.1
    EagerPy version 0.30.0 
    TensorFlow version 2.7.0 
    Sphinx version 4.3.0 
    Sphinx-paramlinks version 0.5.2 
    pydata_sphinx_theme version 0.13.1 
    NetworkX version 3.0 
    Eigen3 version 3.4.0
    Boost version 1.74.0
    CGAL version 5.3
    GMP_LIBRARIES = /usr/lib/x86_64-linux-gnu/libgmp.so
    GMPXX_LIBRARIES = /usr/lib/x86_64-linux-gnu/libgmpxx.so
    MPFR_LIBRARIES = /usr/lib/x86_64-linux-gnu/libmpfr.so
    TBB version 2021.7.0 found and used
    TBB_LIBRARY = /usr/lib/x86_64-linux-gnu/libtbb.so
    TBB_MALLOC_LIBRARY = /usr/lib/x86_64-linux-gnu/libtbbmalloc.so
    + Installed modules are: bottleneck;off_utils;simplex_tree;rips_complex;cubical_complex;periodic_cubical_complex;
        persistence_graphical_tools;reader_utils;witness_complex;strong_witness_complex;nerve_gic;subsampling;
        tangential_complex;alpha_complex;euclidean_witness_complex;euclidean_strong_witness_complex;
    + Missing modules are: 


Documentation
=============

To build the documentation, `sphinx-doc <http://www.sphinx-doc.org>`_,
`sphinxcontrib-bibtex <https://sphinxcontrib-bibtex.readthedocs.io>`_,
`sphinxcontrib-paramlinks <https://github.com/sqlalchemyorg/sphinx-paramlinks>`_ and
`pydata-sphinx-theme <https://github.com/pydata/pydata-sphinx-theme>`_ :math:`\geq` 0.8.0 are
required. As the documentation is auto-tested, `CGAL`_, `Eigen`_,
`Matplotlib`_, `NumPy`_, `POT`_, `Scikit-learn`_ and `SciPy`_ are
also mandatory to build the documentation.

Run the following commands in a terminal:

.. code-block:: bash

    cd /path-to-gudhi/build/python
    make sphinx

Optional third-party library
****************************

Compilation dependencies
========================

These third party dependencies are detected by `CMake <https://www.cmake.org/>`_.
They have to be installed before performing the `GUDHI Python module compilation`_.

CGAL
----

Some GUDHI modules (cf. :doc:`modules list </index>`), and few examples
require `CGAL <https://www.cgal.org/>`_, a C++ library that provides easy
access to efficient and reliable geometric algorithms.


The procedure to install this library
according to your operating system is detailed
`here <http://doc.cgal.org/latest/Manual/installation.html>`_.

The following examples require CGAL version :math:`\geq` 4.11.0:

.. only:: builder_html

    * :download:`alpha_complex_diagram_persistence_from_off_file_example.py <../example/alpha_complex_diagram_persistence_from_off_file_example.py>`
    * :download:`alpha_complex_from_points_example.py <../example/alpha_complex_from_points_example.py>`
    * :download:`bottleneck_basic_example.py <../example/bottleneck_basic_example.py>`
    * :download:`tangential_complex_plain_homology_from_off_file_example.py <../example/tangential_complex_plain_homology_from_off_file_example.py>`
    * :download:`euclidean_strong_witness_complex_diagram_persistence_from_off_file_example.py <../example/euclidean_strong_witness_complex_diagram_persistence_from_off_file_example.py>`
    * :download:`euclidean_witness_complex_diagram_persistence_from_off_file_example.py <../example/euclidean_witness_complex_diagram_persistence_from_off_file_example.py>`

Eigen
-----

Some GUDHI modules (cf. :doc:`modules list </index>`), and few examples
require `Eigen <http://eigen.tuxfamily.org/>`_, a C++ template
library for linear algebra: matrices, vectors, numerical solvers, and related
algorithms.

The following examples require `Eigen <http://eigen.tuxfamily.org/>`_ version :math:`\geq` 3.1.0:

.. only:: builder_html

    * :download:`alpha_complex_diagram_persistence_from_off_file_example.py <../example/alpha_complex_diagram_persistence_from_off_file_example.py>`
    * :download:`alpha_complex_from_points_example.py <../example/alpha_complex_from_points_example.py>`
    * :download:`tangential_complex_plain_homology_from_off_file_example.py <../example/tangential_complex_plain_homology_from_off_file_example.py>`
    * :download:`euclidean_strong_witness_complex_diagram_persistence_from_off_file_example.py <../example/euclidean_strong_witness_complex_diagram_persistence_from_off_file_example.py>`
    * :download:`euclidean_witness_complex_diagram_persistence_from_off_file_example.py <../example/euclidean_witness_complex_diagram_persistence_from_off_file_example.py>`

oneAPI Threading Building Blocks
--------------------------------

`Intel® oneAPI TBB <https://github.com/oneapi-src/oneTBB>`_ lets you easily write parallel C++ programs that take
full advantage of multicore performance, that are portable and composable, and that have future-proof scalability.

Having Intel® oneAPI TBB  (version 20.19.7 or higher) installed is recommended to parallelize and accelerate some GUDHI
computations.

Run time dependencies
=====================

These third party dependencies are detected by Python `import` mechanism at run time.
They can be installed when required.

EagerPy
-------

Some Python functions can handle automatic differentiation (possibly only when
a flag `enable_autodiff=True` is used). In order to reduce code duplication, we
use `EagerPy <https://eagerpy.jonasrauber.de/>`_ which wraps arrays from
PyTorch, TensorFlow and JAX in a common interface.

Joblib
------

`Joblib <https://joblib.readthedocs.io/>`_ is used both as a dependency of `Scikit-learn`_,
and directly for parallelism in some modules (:class:`~gudhi.point_cloud.knn.KNearestNeighbors`,
:func:`~gudhi.representations.metrics.pairwise_persistence_diagram_distances`).

Hnswlib
-------

:class:`~gudhi.point_cloud.knn.KNearestNeighbors` can use the Python package
`Hnswlib <https://github.com/nmslib/hnswlib>`_ as a backend if explicitly
requested, to speed-up queries.

LaTeX
~~~~~

If a sufficiently complete LaTeX toolchain is available (including dvipng and ghostscript), the LaTeX option of
matplotlib is enabled for prettier captions (cf.
`matplotlib text rendering with LaTeX <https://matplotlib.org/3.3.0/tutorials/text/usetex.html>`_).
It also requires `type1cm` LaTeX package (not detected by matplotlib).

If you are facing issues with LaTeX rendering, like this one:

.. code-block:: none

    Traceback (most recent call last):
      File "/usr/lib/python3/dist-packages/matplotlib/texmanager.py", line 302, in _run_checked_subprocess
        report = subprocess.check_output(command,
    ...
    ! LaTeX Error: File `type1cm.sty' not found.
    ...

This is because the LaTeX package is not installed on your system. On Ubuntu systems you can install texlive-full
(for all LaTeX packages), or more specific packages like texlive-latex-extra, cm-super.

You can still deactivate LaTeX rendering by saying:

.. code-block:: python

    import gudhi as gd
    gd.persistence_graphical_tools._gudhi_matplotlib_use_tex=False

Matplotlib
----------

The :doc:`persistence graphical tools </persistence_graphical_tools_user>`
module requires `Matplotlib <http://matplotlib.org>`_, a Python 2D plotting
library which produces publication quality figures in a variety of hardcopy
formats and interactive environments across platforms.

The following examples require the `Matplotlib <http://matplotlib.org>`_:

.. only:: builder_html

    * :download:`alpha_complex_diagram_persistence_from_off_file_example.py <../example/alpha_complex_diagram_persistence_from_off_file_example.py>`
    * :download:`gudhi_graphical_tools_example.py <../example/gudhi_graphical_tools_example.py>`
    * :download:`periodic_cubical_complex_barcode_persistence_from_perseus_file_example.py <../example/periodic_cubical_complex_barcode_persistence_from_perseus_file_example.py>`
    * :download:`rips_complex_diagram_persistence_from_off_file_example.py <../example/rips_complex_diagram_persistence_from_off_file_example.py>`
    * :download:`rips_persistence_diagram.py <../example/rips_persistence_diagram.py>`
    * :download:`rips_complex_diagram_persistence_from_distance_matrix_file_example.py <../example/rips_complex_diagram_persistence_from_distance_matrix_file_example.py>`
    * :download:`tangential_complex_plain_homology_from_off_file_example.py <../example/tangential_complex_plain_homology_from_off_file_example.py>`
    * :download:`euclidean_strong_witness_complex_diagram_persistence_from_off_file_example.py <../example/euclidean_strong_witness_complex_diagram_persistence_from_off_file_example.py>`
    * :download:`euclidean_witness_complex_diagram_persistence_from_off_file_example.py <../example/euclidean_witness_complex_diagram_persistence_from_off_file_example.py>`

NetworkX
--------

The :class:`~gudhi.cover_complex.GraphInducedComplex`, :class:`~gudhi.cover_complex.MapperComplex` and
:class:`~gudhi.cover_complex.NerveComplex` `get_networkx` methods require `NetworkX <https://networkx.org/>`_,
a Python library for studying and visualizing graphs and networks.

PyKeOps
-------

:class:`~gudhi.point_cloud.knn.KNearestNeighbors` can use the Python package
`PyKeOps <https://www.kernel-operations.io/keops/python/>`_ as a backend if
explicitly requested, to speed-up queries using a GPU.

Python Optimal Transport
------------------------

The :doc:`Wasserstein distance </wasserstein_distance_user>`
module requires `POT <https://pythonot.github.io/>`_, a library that provides
several solvers for optimization problems related to Optimal Transport.

PyTorch
-------

`PyTorch <https://pytorch.org/>`_ is currently only used as a dependency of
`PyKeOps`_, and in some tests.

Scikit-learn
------------

The :doc:`persistence representations </representations>` module requires
`scikit-learn <https://scikit-learn.org/>`_, a Python-based ecosystem of
open-source software for machine learning.

Some classes of the :doc:`cover complex </cover_complex_sklearn_isk_ref>` module
also follow and require the `scikit-learn <https://scikit-learn.org/>`_ interface.

:class:`~gudhi.point_cloud.knn.KNearestNeighbors` can use the Python package
`scikit-learn <https://scikit-learn.org/>`_ as a backend if explicitly
requested.

SciPy
-----

The :doc:`persistence graphical tools </persistence_graphical_tools_user>` and
:doc:`Wasserstein distance </wasserstein_distance_user>` modules require `SciPy
<http://scipy.org>`_, a Python-based ecosystem of open-source software for
mathematics, science, and engineering.

:class:`~gudhi.point_cloud.knn.KNearestNeighbors` can use the Python package
`SciPy <http://scipy.org>`_ :math:`\geq` 1.6.0 as a backend if explicitly requested.

TensorFlow
----------

:class:`~gudhi.tensorflow.perslay.Perslay` from the :doc:`persistence representations </representations>` module
requires `TensorFlow <https://www.tensorflow.org/>`_.
The :doc:`cubical complex </cubical_complex_tflow_itf_ref>`, :doc:`simplex tree </ls_simplex_tree_tflow_itf_ref>`
and :doc:`Rips complex </rips_complex_tflow_itf_ref>` modules require `TensorFlow`_
for incorporating them in neural nets. 

`TensorFlow`_ is also used in some automatic differentiation tests.

Bug reports and contributions
*****************************

Please help us improving the quality of the GUDHI library.
You may `report bugs <https://github.com/GUDHI/gudhi-devel/issues>`_ or
`contact us <https://gudhi.inria.fr/contact/>`_ for any suggestions.

GUDHI is open to external contributions. If you want to join our development team, please take some time to read our
`contributing guide <https://github.com/GUDHI/gudhi-devel/blob/master/.github/CONTRIBUTING.md>`_.
