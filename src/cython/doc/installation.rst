Installation
############

Compiling
*********
The library uses c++11 and requires `Boost <https://www.boost.org/>`_ ≥ 1.48.0
and `CMake <https://www.cmake.org/>`_   ≥ 3.1.
It is a multi-platform library and compiles on Linux, Mac OSX and Visual
Studio 2015.
It also requires cmake to generate makefiles, and cython to compile the
library.

On `Windows <https://wiki.python.org/moin/WindowsCompilers>`_ , only Python
3.5 and 3.6 are available because of the required Visual Studio version.

On other systems, if you have several Python/cython installed, the version 2.X
will be used by default, but you can force it by adding
:code:`-DPython_ADDITIONAL_VERSIONS=3` to the cmake command.

GUDHI Cythonization
===================

To build the GUDHI cython module, run the following commands in a terminal:

.. code-block:: bash

    cd /path-to-gudhi/
    mkdir build
    cd build/
    cmake ..
    make cython

Test suites
===========

To test your build, `py.test <http://doc.pytest.org>`_ is optional. Run the
following command in a terminal:

.. code-block:: bash

    cd /path-to-gudhi/build/cython
    # For windows, you have to set PYTHONPATH environment variable
    export PYTHONPATH='$PYTHONPATH:/path-to-gudhi/build/cython'
    ctest -R py_test

If tests fail, please try to :code:`import gudhi` and check the errors.
The problem can come from a third-party library bad link or installation.

Documentation
=============

To build the documentation, `sphinx-doc <http://http://www.sphinx-doc.org>`_ is
required. Please refer to *conf.py* file to see which
`sphinx-doc <http://http://www.sphinx-doc.org>`_ modules are required to
generate the documentation. Run the following commands in a terminal:

.. code-block:: bash

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

Having CGAL,  the Computational Geometry Algorithms Library, version 4.7.0 or
higher installed is recommended. The procedure to install this library
according to your operating system is detailed
`here <http://doc.cgal.org/latest/Manual/installation.html>`_.

The following examples requires CGAL version ≥ 4.7.0:

.. only:: builder_html

    * :download:`alpha_complex_diagram_persistence_from_off_file_example.py <../example/alpha_complex_diagram_persistence_from_off_file_example.py>`
    * :download:`alpha_complex_from_points_example.py <../example/alpha_complex_from_points_example.py>`

The following examples requires CGAL version ≥ 4.8.0:

.. only:: builder_html

    * :download:`bottleneck_basic_example.py <../example/bottleneck_basic_example.py>`
    * :download:`tangential_complex_plain_homology_from_off_file_example.py <../example/tangential_complex_plain_homology_from_off_file_example.py>`

The following examples requires CGAL version ≥ 4.8.1:

.. only:: builder_html

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

Numpy
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
