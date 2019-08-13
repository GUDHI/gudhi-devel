:orphan:

.. To get rid of WARNING: document isn't included in any toctree

Installation
############

Compiling
*********
The library uses c++11 and requires `Boost <https://www.boost.org/>`_ ≥ 1.56.0,
`CMake <https://www.cmake.org/>`_   ≥ 3.1  to generate makefiles, and
`Cython <https://www.cython.org/>`_ to compile the GUDHI Python module.
It is a multi-platform library and compiles on Linux, Mac OSX and Visual
Studio 2015.

On `Windows <https://wiki.python.org/moin/WindowsCompilers>`_ , only Python
3.5 and 3.6 are available because of the required Visual Studio version.

On other systems, if you have several Python/Cython installed, the version 2.X
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
    cd cython
    make

GUDHI Python module installation
================================

Once the compilation succeeds, one can add the GUDHI Python module path to the
PYTHONPATH:

.. code-block:: bash

    # For windows, you have to set PYTHONPATH environment variable
    export PYTHONPATH='$PYTHONPATH:/path-to-gudhi/build/cython'

Or install it definitely in your Python packages folder:

.. code-block:: bash

    cd /path-to-gudhi/build/cython
    # May require sudo or administrator privileges
    make install


Test suites
===========

To test your build, `py.test <http://doc.pytest.org>`_ is optional. Run the
following command in a terminal:

.. code-block:: bash

    cd /path-to-gudhi/build/cython
    # For windows, you have to set PYTHONPATH environment variable
    export PYTHONPATH='$PYTHONPATH:/path-to-gudhi/build/cython'
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
    Eigen3 version 3.1.1
    Installed modules are: off_reader;simplex_tree;rips_complex;cubical_complex;periodic_cubical_complex;
        persistence_graphical_tools;reader_utils;witness_complex;strong_witness_complex;alpha_complex;
        euclidean_witness_complex;euclidean_strong_witness_complex;
    Missing modules are: bottleneck_distance;nerve_gic;subsampling;tangential_complex;persistence_graphical_tools;
    CGAL version 4.7.1000
    GMP_LIBRARIES = /usr/lib/x86_64-linux-gnu/libgmp.so
    GMPXX_LIBRARIES = /usr/lib/x86_64-linux-gnu/libgmpxx.so
    TBB version 9107 found and used

Here, you can see that bottleneck_distance, nerve_gic, subsampling and
tangential_complex are missing because of the CGAL version.
persistence_graphical_tools is not available as numpy and matplotlib are not
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
    Installed modules are: off_reader;simplex_tree;rips_complex;cubical_complex;periodic_cubical_complex;
        persistence_graphical_tools;reader_utils;witness_complex;strong_witness_complex;persistence_graphical_tools;
        bottleneck_distance;nerve_gic;subsampling;tangential_complex;alpha_complex;euclidean_witness_complex;
        euclidean_strong_witness_complex;
    CGAL header only version 4.11.0
    GMP_LIBRARIES = /usr/lib/x86_64-linux-gnu/libgmp.so
    GMPXX_LIBRARIES = /usr/lib/x86_64-linux-gnu/libgmpxx.so
    TBB version 9107 found and used

Documentation
=============

To build the documentation, `sphinx-doc <http://www.sphinx-doc.org>`_ and
`sphinxcontrib-bibtex <https://sphinxcontrib-bibtex.readthedocs.io>`_ are
required. As the documentation is auto-tested, `CGAL`_, `Eigen3`_,
`Matplotlib`_, `NumPy`_ and `SciPy`_ are also mandatory to build the
documentation.

Run the following commands in a terminal:

.. code-block:: bash

    cd /path-to-gudhi/build/cython
    make sphinx

Optional third-party library
****************************

CGAL
====

The :doc:`Alpha complex </alpha_complex_user>`,
:doc:`Tangential complex </tangential_complex_user>` and
:doc:`Witness complex </witness_complex_user>` data structures, and
:doc:`Bottleneck distance </bottleneck_distance_user>` requires CGAL, which is a
C++ library which provides easy access to efficient and reliable geometric
algorithms.

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

Eigen3
======

The :doc:`Alpha complex </alpha_complex_user>`,
:doc:`Tangential complex </tangential_complex_user>` and
:doc:`Witness complex </witness_complex_user>` data structures and few
examples requires `Eigen3 <http://eigen.tuxfamily.org/>`_, a C++ template
library for linear algebra: matrices, vectors, numerical solvers, and related
algorithms.

The following examples require the `Eigen3 <http://eigen.tuxfamily.org/>`_:

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

NumPy
=====

The :doc:`persistence graphical tools </persistence_graphical_tools_user>`
module requires `NumPy <http://numpy.org>`_, a fundamental package for
scientific computing with Python.

The following examples require the `NumPy <http://numpy.org>`_:

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

SciPy
=====

The :doc:`persistence graphical tools </persistence_graphical_tools_user>`
module requires `SciPy <http://scipy.org>`_, a Python-based ecosystem of
open-source software for mathematics, science, and engineering.

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
