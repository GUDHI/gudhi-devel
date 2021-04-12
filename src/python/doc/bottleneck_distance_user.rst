:orphan:

.. To get rid of WARNING: document isn't included in any toctree

Bottleneck distance user manual
===============================
Definition
----------

.. include:: bottleneck_distance_sum.inc

This implementation by François Godi is based on ideas from "Geometry Helps in Bottleneck Matching and Related Problems"
:cite:`DBLP:journals/algorithmica/EfratIK01` and requires `CGAL <installation.html#cgal>`_ (`GPL v3 </licensing/>`_).

.. autofunction:: gudhi.bottleneck_distance

This other implementation comes from `Hera
<https://bitbucket.org/grey_narn/hera/src/master/>`_ (BSD-3-Clause) which is
based on "Geometry Helps to Compare Persistence Diagrams"
:cite:`Kerber:2017:GHC:3047249.3064175` by Michael Kerber, Dmitriy
Morozov, and Arnur Nigmetov.

.. warning::
   Beware that its approximation allows for a multiplicative error, while the function above uses an additive error.

.. autofunction:: gudhi.hera.bottleneck_distance


Distance computation
--------------------

The following example explains how the distance is computed:

.. testcode::

    import gudhi

    message = "Bottleneck distance = " + '%.1f' % gudhi.bottleneck_distance([[0., 0.]], [[0., 13.]])
    print(message)

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

    diag1 = [[2.7, 3.7],[9.6, 14.],[34.2, 34.974], [3.,float('Inf')]]
    diag2 = [[2.8, 4.45],[9.5, 14.1],[3.2,float('Inf')]]

    message = "Bottleneck distance approximation = " + '%.2f' % gudhi.bottleneck_distance(diag1, diag2, 0.1)
    print(message)

    message = "Bottleneck distance value = " + '%.2f' % gudhi.bottleneck_distance(diag1, diag2)
    print(message)

The output is:

.. testoutput::

    Bottleneck distance approximation = 0.72
    Bottleneck distance value = 0.75

