Bottleneck distance user manual
===============================
Definition
----------

.. include:: bottleneck_distance_sum.rst

Function
--------
.. autofunction:: gudhi.bottleneck_distance


Basic example
-------------

This example computes the bottleneck distance from 2 persistence diagrams:

.. testcode::

    import gudhi

    diag1 = [[2.7, 3.7],[9.6, 14.],[34.2, 34.974], [3.,float('Inf')]]
    diag2 = [[2.8, 4.45],[9.5, 14.1],[3.2,float('Inf')]]

    message = "Bottleneck distance approximation=" + '%.2f' % gudhi.bottleneck_distance(diag1, diag2, 0.1)
    print(message)

    message = "Bottleneck distance exact value=" + '%.2f' % gudhi.bottleneck_distance(diag1, diag2, 0)
    print(message)

The output is:

.. testoutput::

    Bottleneck distance approximation=0.81
    Bottleneck distance exact value=0.75
