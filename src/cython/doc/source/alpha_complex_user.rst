=========================
Alpha complex user manual
=========================
Definition
----------

.. include:: alpha_complex_sum.rst

Alpha_complex is constructing a Simplex_tree using Delaunay Triangulation from
CGAL (the Computational Geometry Algorithms Library).

The complex is a template class requiring an Epick_d dD Geometry Kernel [16] from CGAL as template parameter.

Remarks
    When Alpha_complex is constructed with an infinite value of alpha, the complex is a Delaunay complex.