:orphan:

.. To get rid of WARNING: document isn't included in any toctree

======================
Representations manual
======================

.. include:: representations_sum.inc

This module, originally available at https://github.com/MathieuCarriere/sklearn-tda and named sklearn_tda, aims at bridging the gap between persistence diagrams and machine learning, by providing implementations of most of the vector representations for persistence diagrams in the literature, in a scikit-learn format. More specifically, it provides tools, using the scikit-learn standard interface, to compute distances and kernels on persistence diagrams, and to convert these diagrams into vectors in Euclidean space.

A diagram is represented as a numpy array of shape (n,2), as can be obtained from :func:`~gudhi.SimplexTree.persistence_intervals_in_dimension` for instance. Points at infinity are represented as a numpy array of shape (n,1), storing only the birth time.

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

Basic example
-------------

This example computes the first two Landscapes associated to a persistence diagram with four points. The landscapes are evaluated on ten samples, leading to two vectors with ten coordinates each, that are eventually concatenated in order to produce a single vector representation.

.. testcode::

    import numpy as np
    from gudhi.representations import Landscape
    # A single diagram with 4 points
    D = np.array([[0.,4.],[1.,2.],[3.,8.],[6.,8.]])
    diags = [D]
    l=Landscape(num_landscapes=2,resolution=10).fit_transform(diags)
    print(l) 

The output is:

.. testoutput::

    [[1.02851895 2.05703791 2.57129739 1.54277843 0.89995409 1.92847304
      2.95699199 3.08555686 0.         0.64282435 0.         0.
      0.51425948 0.         0.         0.        ]]
