:orphan:

.. To get rid of WARNING: document isn't included in any toctree

Installation
############

Conda
*****
The easiest way to install the Python version of GUDHI is using
`conda <https://gudhi.inria.fr/conda/>`_.

Compiling
*********
The library uses c++14 and requires `Boost <https://www.boost.org/>`_ ≥ 1.56.0,
`CMake <https://www.cmake.org/>`_ ≥ 3.1  to generate makefiles,
`NumPy <http://numpy.org>`_ and `Cython <https://www.cython.org/>`_ to compile
the GUDHI Python module.
It is a multi-platform library and compiles on Linux, Mac OSX and Visual
Studio 2015.

On `Windows <https://wiki.python.org/moin/WindowsCompilers>`_ , only Python
≥ 3.5 are available because of the required Visual Studio version.

On other systems, if you have several Python/python installed, the version 2.X
will be used by default, but you can force it by adding
:code:`-DPython_ADDITIONAL_VERSIONS=3` to the cmake command.

GUDHI Python module compilation
===============================

To build the GUDHI Python module, run the following commands in a terminal:

.. code-block:: bash

    cd /path-to-gudhi/
    mkdir build
    cd build/
    cmake ..
    cd python
    make

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


Test suites
===========

To test your build, `py.test <http://doc.pytest.org>`_ is optional. Run the
following command in a terminal:

.. code-block:: bash

    cd /path-to-gudhi/build/python
    # For windows, you have to set PYTHONPATH environment variable
    export PYTHONPATH='$PYTHONPATH:/path-to-gudhi/build/python'
    make test

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
`Matplotlib`_, `NumPy`_ and `SciPy`_ are also mandatory to build the
documentation.

Run the following commands in a terminal:

.. code-block:: bash

    cd /path-to-gudhi/build/python
    make sphinx

Optional third-party library
****************************

CGAL
====

Some GUDHI modules (cf. :doc:`modules list </index>`), and few examples
require CGAL, a C++ library that provides easy access to efficient and
reliable geometric algorithms.


The procedure to install this library
according to your operating system is detailed
`here <http://doc.cgal.org/latest/Manual/installation.html>`_.

The following examples requires CGAL version ≥ 4.11.0:

.. only:: builder_html

    * :download:`alpha_complex_diagram_persistence_from_off_file_example.py <../example/alpha_complex_diagram_persistence_from_off_file_example.py>`
    * :download:`alpha_complex_from_points_example.py <../example/alpha_complex_from_points_example.py>`
    * :download:`bottleneck_basic_example.py <../example/bottleneck_basic_example.py>`
    * :download:`tangential_complex_plain_homology_from_off_file_example.py <../example/tangential_complex_plain_homology_from_off_file_example.py>`
    * :download:`euclidean_strong_witness_complex_diagram_persistence_from_off_file_example.py <../example/euclidean_strong_witness_complex_diagram_persistence_from_off_file_example.py>`
    * :download:`euclidean_witness_complex_diagram_persistence_from_off_file_example.py <../example/euclidean_witness_complex_diagram_persistence_from_off_file_example.py>`

Eigen
=====

Some GUDHI modules (cf. :doc:`modules list </index>`), and few examples
require `Eigen <http://eigen.tuxfamily.org/>`_, a C++ template
library for linear algebra: matrices, vectors, numerical solvers, and related
algorithms.

The following examples require `Eigen <http://eigen.tuxfamily.org/>`_ version ≥ 3.1.0:

.. only:: builder_html

    * :download:`alpha_complex_diagram_persistence_from_off_file_example.py <../example/alpha_complex_diagram_persistence_from_off_file_example.py>`
    * :download:`alpha_complex_from_points_example.py <../example/alpha_complex_from_points_example.py>`
    * :download:`tangential_complex_plain_homology_from_off_file_example.py <../example/tangential_complex_plain_homology_from_off_file_example.py>`
    * :download:`euclidean_strong_witness_complex_diagram_persistence_from_off_file_example.py <../example/euclidean_strong_witness_complex_diagram_persistence_from_off_file_example.py>`
    * :download:`euclidean_witness_complex_diagram_persistence_from_off_file_example.py <../example/euclidean_witness_complex_diagram_persistence_from_off_file_example.py>`

Matplotlib
==========

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

Python Optimal Transport
========================

The :doc:`Wasserstein distance </wasserstein_distance_user>`
module requires `POT <https://pot.readthedocs.io/>`_, a library that provides
several solvers for optimization problems related to Optimal Transport.

SciPy
=====

The :doc:`persistence graphical tools </persistence_graphical_tools_user>` and
:doc:`Wasserstein distance </wasserstein_distance_user>` modules require `SciPy
<http://scipy.org>`_, a Python-based ecosystem of open-source software for
mathematics, science, and engineering.

Threading Building Blocks
=========================

`Intel® TBB <https://www.threadingbuildingblocks.org/>`_ lets you easily write
parallel C++ programs that take full advantage of multicore performance, that
are portable and composable, and that have future-proof scalability.

Having Intel® TBB installed is recommended to parallelize and accelerate some
GUDHI computations.

Bug reports and contributions
*****************************

Please help us improving the quality of the GUDHI library. You may report bugs or suggestions to:

    Contact: gudhi-users@lists.gforge.inria.fr

GUDHI is open to external contributions. If you want to join our development team, please contact us.
