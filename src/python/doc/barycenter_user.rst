:orphan:

.. To get rid of WARNING: document isn't included in any toctree

Barycenter user manual
================================
Definition
----------

.. include:: barycenter_sum.inc

This implementation is based on ideas from "Frechet means for distribution of 
persistence diagrams", Turner et al. 2014.

Function
--------
.. autofunction:: gudhi.barycenter.lagrangian_barycenter


Basic example
-------------

This example computes the Frechet mean (aka Wasserstein barycenter) between 
four persistence diagrams.
It is initialized on the 4th diagram.
As the algorithm is not convex, its output depends on the initialization and 
is only a local minimum of the objective function.
Initialization can be either given as an integer (in which case the i-th 
diagram of the list is used as initial estimate) or as a diagram. 
If None, it will randomly select one of the diagram of the list 
as initial estimate.
Note that persistence diagrams must be submitted as 
(n x 2) numpy arrays and must not contain inf values.

.. testcode::

    import gudhi.barycenter
    import numpy as np

    dg1 = np.array([[0.2, 0.5]])
    dg2 = np.array([[0.2, 0.7]])
    dg3 = np.array([[0.3, 0.6], [0.7, 0.8], [0.2, 0.3]])
    dg4 = np.array([])
    pdiagset = [dg1, dg2, dg3, dg4]
    bary = gudhi.barycenter.lagrangian_barycenter(pdiagset=pdiagset,init=3)

    message = "Wasserstein barycenter estimated:"    
    print(message)
    print(bary)

The output is:

.. testoutput::

    Wasserstein barycenter estimated:
    [[0.27916667 0.55416667]
     [0.7375     0.7625    ]
     [0.2375     0.2625    ]]


