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

This example computes the 1-Wasserstein distance from 2 persistence diagrams with Euclidean ground metric.
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

We can also have access to the optimal matching by letting `matching=True`. 
It is encoded as a list of indices (i,j), meaning that the i-th point in X
is mapped to the j-th point in Y. 
An index of -1 represents the diagonal.

.. testcode::

    import gudhi.wasserstein
    import numpy as np

    dgm1 = np.array([[2.7, 3.7],[9.6, 14.],[34.2, 34.974]])
    dgm2 = np.array([[2.8, 4.45], [5, 6], [9.5, 14.1]])
    cost, matchings = gudhi.wasserstein.wasserstein_distance(diag1, diag2, matching=True, order=1., internal_p=2.)

    message_cost = "Wasserstein distance value = %.2f" %cost
    print(message_cost)
    dgm1_to_diagonal = matchings[np.where(matchings[:,0] == -1)][:,1]
    dgm2_to_diagonal = matchings[np.where(matchings[:,1] == -1)][:,0]
    off_diagonal_match = np.delete(matchings, np.where(matchings == -1)[0], axis=0)

    for i,j in off_diagonal_match:
        print("point %s in dgm1 is matched to point %s in dgm2" %(i,j))
    for i in dgm1_to_diagonal:
        print("point %s in dgm1 is matched to the diagonal" %i)
    for j in dgm2_to_diagonal:
        print("point %s in dgm2 is matched to the diagonal" %j)

The output is:

.. testoutput::
    
    Wasserstein distance value = 2.15
    point 0 in dgm1 is matched to point 0 in dgm2
    point 1 in dgm1 is matched to point 2 in dgm2
    point 2 in dgm1 is matched to the diagonal
    point 1 in dgm2 is matched to the diagonal
