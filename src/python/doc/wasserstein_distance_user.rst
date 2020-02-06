:orphan:

.. To get rid of WARNING: document isn't included in any toctree

Wasserstein distance user manual
================================
Definition
----------

.. include:: wasserstein_distance_sum.inc

Functions
---------
This implementation uses the Python Optimal Transport library and is based on
ideas from "Large Scale Computation of Means and Cluster for Persistence
Diagrams via Optimal Transport" :cite:`10.5555/3327546.3327645`.

.. autofunction:: gudhi.wasserstein.wasserstein_distance

This other implementation comes from `Hera
<https://bitbucket.org/grey_narn/hera/src/master/>`_ (BSD-3-Clause) which is
based on "Geometry Helps to Compare Persistence Diagrams"
:cite:`Kerber:2017:GHC:3047249.3064175` by Michael Kerber, Dmitriy
Morozov, and Arnur Nigmetov.

.. autofunction:: gudhi.hera.wasserstein_distance

Basic example
-------------

This example computes the 1-Wasserstein distance from 2 persistence diagrams with euclidean ground metric.
Note that persistence diagrams must be submitted as (n x 2) numpy arrays and must not contain inf values.

.. testcode::

    import gudhi.wasserstein
    import numpy as np

    diag1 = np.array([[2.7, 3.7],[9.6, 14.],[34.2, 34.974]])
    diag2 = np.array([[2.8, 4.45],[9.5, 14.1]])

    message = "Wasserstein distance value = " + '%.2f' % gudhi.wasserstein.wasserstein_distance(diag1, diag2, order=1., internal_p=2.)
    print(message)

The output is:

.. testoutput::

    Wasserstein distance value = 1.45
