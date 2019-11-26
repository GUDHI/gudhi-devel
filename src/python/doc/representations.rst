:orphan:

.. To get rid of WARNING: document isn't included in any toctree

======================
Representations manual
======================

.. include:: representations_sum.inc

This module, originally named sklearn_tda, aims at bridging the gap between persistence diagrams and machine learning tools, in particular scikit-learn. It provides tools, using the scikit-learn standard interface, to compute distances and kernels on diagrams, and to convert diagrams into vectors.

A diagram is represented as a numpy array of shape (n,2), as can be obtained from :func:`gudhi.SimplexTree.persistence_intervals_in_dimension` for instance. Points at infinity are represented as a numpy array of shape (n,1), storing only the birth time.

A small example is provided

.. only:: builder_html

    * :download:`diagram_vectorizations_distances_kernels.py <../example/diagram_vectorizations_distances_kernels.py>`


Preprocessing
-------------
.. automodule:: gudhi.representations.preprocessing
   :members:
   :special-members:
   :show-inheritance:

Vector methods
--------------
.. automodule:: gudhi.representations.vector_methods
   :members:
   :special-members:
   :show-inheritance:

Kernel methods
--------------
.. automodule:: gudhi.representations.kernel_methods
   :members:
   :special-members:
   :show-inheritance:

Metrics
-------
.. automodule:: gudhi.representations.metrics
   :members:
   :special-members:
   :show-inheritance:
