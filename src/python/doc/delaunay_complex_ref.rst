:orphan:

.. To get rid of WARNING: document isn't included in any toctree

=================================
Delaunay complex reference manual
=================================

.. autofunction:: gudhi.delaunay_complex

======================================
Delaunay Čech complex reference manual
======================================

.. autofunction:: gudhi.delaunay_cech_complex

==============================
Alpha complex reference manual
==============================

.. autofunction:: gudhi.alpha_complex

.. autofunction:: gudhi.weighted_alpha_complex

.. autoclass:: gudhi.AlphaComplex

   .. deprecated:: 3.11.0
      Use :py:func:`gudhi.alpha_complex`, :py:func:`gudhi.weighted_alpha_complex` or :py:func:`gudhi.delaunay_complex`
      instead.

   .. automethod:: gudhi.AlphaComplex.create_simplex_tree
   .. automethod:: gudhi.AlphaComplex.get_point

==================
Relative precision
==================

.. autofunction:: gudhi.DelaunayComplex.set_float_relative_precision

.. autofunction:: gudhi.DelaunayComplex.get_float_relative_precision
