Installation
############

Compiling
*********

The library uses c++11 and requires `Boost <http://www.boost.org/>`_ with version
1.48.0 or more recent. It is a multi-platform library and compiles on Linux,
Mac OSX and Visual Studio 2015.
It also requires cmake to generate makefiles, and cython to compile the
library.

Demos and examples
==================

To build the demos and examples, run the following commands in a terminal:

.. code-block:: bash

    cd /path-to-gudhi/
    mkdir build
    cd build/
    cmake ..
    make cython

A list of examples is available here.

Test suites
===========

To test your build, `py.test <http://doc.pytest.org>`_ is required. Run the
following command in a terminal:

.. code-block:: bash

    cd /path-to-gudhi/build/src/cython
    # For windows, you have to set PYTHONPATH environment variable
    export PYTHONPATH='$PYTHONPATH:/path-to-gudhi/build/src/cython'
    py.test

Documentation
=============

To build the documentation, `sphinx-doc <http://http://www.sphinx-doc.org>`_ is
required. Please refer to *conf.py* file to see which
`sphinx-doc <http://http://www.sphinx-doc.org>`_ modules are required to
generate the documentation. Run the following commands in a terminal:

.. code-block:: bash

    make sphynx

Optional third-party library
****************************

CGAL
====

The :doc:`Alpha complex </alpha_complex_user>` data structure and
:doc:`Tangential complex </tangential_complex_user>` data structures,
:doc:`Bottleneck distance </bottleneck_distance_user>` requires CGAL, which is a
C++ library which provides easy access to efficient and reliable geometric
algorithms.

Having CGAL version 4.7.0 or higher installed is recommended. The procedure to
install this library according to your operating system is detailed
`here <http://doc.cgal.org/latest/Manual/installation.html>`_.

The following examples require the Computational Geometry Algorithms Library:

.. only:: builder_html

    * :download:`alpha_complex_diagram_persistence_from_off_file_example.py <../example/alpha_complex_diagram_persistence_from_off_file_example.py>`
    * :download:`alpha_complex_from_points_example.py <../example/alpha_complex_from_points_example.py>`

The following example requires CGAL version ≥ 4.8.0:

.. only:: builder_html

    * :download:`bottleneck_basic_example.py <../example/bottleneck_basic_example.py>`
    * :download:`tangential_complex_plain_homology_from_off_file_example.py <../example/tangential_complex_plain_homology_from_off_file_example.py>`

Eigen3
======

The :doc:`Alpha complex </alpha_complex_user>` and
:doc:`Tangential complex </tangential_complex_user>` data structures and few
examples requires `Eigen3 <http://eigen.tuxfamily.org/>`_, a C++ template
library for linear algebra: matrices, vectors, numerical solvers, and related
algorithms.

The following examples require the `Eigen3 <http://eigen.tuxfamily.org/>`_:

.. only:: builder_html

    * :download:`alpha_complex_diagram_persistence_from_off_file_example.py <../example/alpha_complex_diagram_persistence_from_off_file_example.py>`
    * :download:`alpha_complex_from_points_example.py <../example/alpha_complex_from_points_example.py>`
    * :download:`tangential_complex_plain_homology_from_off_file_example.py <../example/tangential_complex_plain_homology_from_off_file_example.py>`

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
