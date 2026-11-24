:orphan:

.. To get rid of WARNING: document isn't included in any toctree

Bottleneck distance user manual
===============================
Definition
----------

.. include:: bottleneck_distance_sum.inc

This implementation by François Godi is based on ideas from "Geometry Helps in Bottleneck Matching and Related Problems"
:cite:`DBLP:journals/algorithmica/EfratIK01` and requires `CGAL <installation.html#cgal>`_ (`GPL v3 </licensing/>`_).

.. autofunction:: gudhi.diagram_distances.bottleneck_distance

This other implementation comes from `Hera
<https://bitbucket.org/grey_narn/hera/src/master/>`_ (BSD-3-Clause) which is
based on "Geometry Helps to Compare Persistence Diagrams"
:cite:`Kerber:2017:GHC:3047249.3064175` by Michael Kerber, Dmitriy
Morozov, and Arnur Nigmetov.

.. warning::
   Beware that its approximation allows for a multiplicative error, while the function above uses an additive error.

.. autofunction:: gudhi.diagram_distances.hera.bottleneck_distance


Distance computation
--------------------

The following example explains how the distance is computed:

.. testcode::

    import gudhi
    import numpy as np

    print(f"Bottleneck distance = {gudhi.diagram_distances.bottleneck_distance([[0., 0.]], [[0., 13.]]):.1f}")

.. testoutput::

    Bottleneck distance = 6.5

.. figure::
    ../../doc/Bottleneck_distance/bottleneck_distance_example.png
    :figclass: align-center
    
    The point (0, 13) is at distance 6.5 from the diagonal and more
    specifically from the point (6.5, 6.5).


Basic example
-------------

This other example computes the bottleneck distance from 2 persistence diagrams:

.. testcode::

    import gudhi
    import numpy as np

    diag1 = np.array([[2.7, 3.7],[9.6, 14.],[34.2, 34.974], [3.,float('Inf')]])
    diag2 = np.array([[2.8, 4.45],[9.5, 14.1],[3.2,float('Inf')]])

    approximation = gudhi.diagram_distances.bottleneck_distance(diag1, diag2, 0.1)
    print(f"Bottleneck distance approximation = {approximation:.2f}")

    distance = gudhi.diagram_distances.bottleneck_distance(diag1, diag2)
    print(f"Bottleneck distance value = {distance:.2f}")

The output is:

.. testoutput::

    Bottleneck distance approximation = 0.72
    Bottleneck distance value = 0.75
