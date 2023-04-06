:orphan:

.. To get rid of WARNING: document isn't included in any toctree

TensorFlow layer for Vietoris-Rips persistence
##############################################

.. include:: differentiation_sum.inc

Example of gradient computed from Vietoris-Rips persistence
-----------------------------------------------------------

.. testsetup::

    import numpy
    numpy.set_printoptions(precision=4)

.. testcode::

    from gudhi.tensorflow import RipsLayer
    import tensorflow as tf

    X = tf.Variable([[1.,1.],[2.,2.]], dtype=tf.float32, trainable=True)
    rl = RipsLayer(maximum_edge_length=2., homology_dimensions=[0])

    with tf.GradientTape() as tape:
        dgm = rl.call(X)[0][0]
        loss = tf.math.reduce_sum(tf.square(.5*(dgm[:,1]-dgm[:,0])))

    grads = tape.gradient(loss, [X])
    print(grads[0].numpy())

.. testcleanup::

    numpy.set_printoptions(precision=8)

.. testoutput::

    [[-0.5 -0.5]
     [ 0.5  0.5]]

Documentation for RipsLayer
---------------------------

.. autoclass:: gudhi.tensorflow.RipsLayer
   :members:
   :show-inheritance:
