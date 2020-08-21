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
with a name like `gudhi.3.2.0.tar.gz`, then run `tar xf gudhi.3.2.0.tar.gz`, which
created a directory `gudhi.3.2.0`, hereinafter referred to as `/path-to-gudhi/`.
If you are instead using a git checkout, beware that the paths are a bit
different, and in particular the `python/` subdirectory is actually `src/python/`
there.

The library uses c++14 and requires `Boost <https://www.boost.org/>`_ :math:`\geq` 1.56.0,
`CMake <https://www.cmake.org/>`_ :math:`\geq` 3.1  to generate makefiles,
`NumPy <http://numpy.org>`_, `Cython <https://www.cython.org/>`_ and
`pybind11 <https://github.com/pybind/pybind11>`_ to compile
the GUDHI Python module.
It is a multi-platform library and compiles on Linux, Mac OSX and Visual
Studio 2017 or later.

On `Windows <https://wiki.python.org/moin/WindowsCompilers>`_ , only Python
:math:`\geq` 3.5 are available because of the required Visual Studio version.

On other systems, if you have several Python/python installed, the version 2.X
will be used by default, but you can force it by adding
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
    # May require sudo or administrator privileges
    make install

.. note::

    :code:`make install` is only a
    `CMake custom targets <https://cmake.org/cmake/help/latest/command/add_custom_target.html>`_
    to shortcut :code:`python setup.py install` command.
    It does not take into account :code:`CMAKE_INSTALL_PREFIX`.
    But one can use :code:`python setup.py install ...` specific options in the python directory:

.. code-block:: bash

    python setup.py install --prefix /home/gudhi  # Install in /home/gudhi directory

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

    import gudhi
    print(gudhi.__debug_info__)

You shall have something like:

.. code-block:: none

    Python version 2.7.15
    Cython version 0.26.1
    Numpy version 1.14.1
    Eigen3 version 3.1.1
    Installed modules are: off_reader;simplex_tree;rips_complex;
        cubical_complex;periodic_cubical_complex;reader_utils;witness_complex;
        strong_witness_complex;alpha_complex;
    Missing modules are: bottleneck_distance;nerve_gic;subsampling;
        tangential_complex;persistence_graphical_tools;
        euclidean_witness_complex;euclidean_strong_witness_complex;
    CGAL version 4.7.1000
    GMP_LIBRARIES = /usr/lib/x86_64-linux-gnu/libgmp.so
    GMPXX_LIBRARIES = /usr/lib/x86_64-linux-gnu/libgmpxx.so
    TBB version 9107 found and used

Here, you can see that bottleneck_distance, nerve_gic, subsampling and
tangential_complex are missing because of the CGAL version.
persistence_graphical_tools is not available as matplotlib is not
available.
Unitary tests cannot be run as pytest is missing.

A complete configuration would be :

.. code-block:: none

    Python version 3.6.5
    Cython version 0.28.2
    Pytest version 3.3.2
    Matplotlib version 2.2.2
    Numpy version 1.14.5
    Eigen3 version 3.3.4
    Installed modules are: off_reader;simplex_tree;rips_complex;
        cubical_complex;periodic_cubical_complex;persistence_graphical_tools;
        reader_utils;witness_complex;strong_witness_complex;
        persistence_graphical_tools;bottleneck_distance;nerve_gic;subsampling;
        tangential_complex;alpha_complex;euclidean_witness_complex;
        euclidean_strong_witness_complex;
    CGAL header only version 4.11.0
    GMP_LIBRARIES = /usr/lib/x86_64-linux-gnu/libgmp.so
    GMPXX_LIBRARIES = /usr/lib/x86_64-linux-gnu/libgmpxx.so
    TBB version 9107 found and used

Documentation
=============

To build the documentation, `sphinx-doc <http://www.sphinx-doc.org>`_ and
`sphinxcontrib-bibtex <https://sphinxcontrib-bibtex.readthedocs.io>`_ are
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

Threading Building Blocks
-------------------------

`Intel® TBB <https://www.threadingbuildingblocks.org/>`_ lets you easily write
parallel C++ programs that take full advantage of multicore performance, that
are portable and composable, and that have future-proof scalability.

Having Intel® TBB installed is recommended to parallelize and accelerate some
GUDHI computations.

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

    import gudhi
    gudhi.persistence_graphical_tools._gudhi_matplotlib_use_tex=False

PyKeOps
-------

:class:`~gudhi.point_cloud.knn.KNearestNeighbors` can use the Python package
`PyKeOps <https://www.kernel-operations.io/keops/python/>`_ as a backend if
explicitly requested, to speed-up queries using a GPU.

Python Optimal Transport
------------------------

The :doc:`Wasserstein distance </wasserstein_distance_user>`
module requires `POT <https://pot.readthedocs.io/>`_, a library that provides
several solvers for optimization problems related to Optimal Transport.

PyTorch
-------

`PyTorch <https://pytorch.org/>`_ is currently only used as a dependency of
`PyKeOps`_, and in some tests.

Scikit-learn
------------

The :doc:`persistence representations </representations>` module require
`scikit-learn <https://scikit-learn.org/>`_, a Python-based ecosystem of
open-source software for machine learning.

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
`SciPy <http://scipy.org>`_ as a backend if explicitly requested.

Bug reports and contributions
*****************************

Please help us improving the quality of the GUDHI library. You may report bugs or suggestions to:

    Contact: gudhi-users@lists.gforge.inria.fr

GUDHI is open to external contributions. If you want to join our development team, please contact us.
