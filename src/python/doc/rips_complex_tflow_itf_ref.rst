:orphan:

.. To get rid of WARNING: document isn't included in any toctree

TensorFlow layer for Vietoris-Rips persistence
##############################################

.. include:: differentiation_sum.inc

Example of gradient computed from Vietoris-Rips persistence
-----------------------------------------------------------

.. code-block:: python
    from gudhi.tensorflow import *
    import numpy as np
    import tensorflow as tf

    Xinit = np.array([[1.,1.],[2.,2.]], dtype=np.float32)
    X = tf.Variable(initial_value=Xinit, trainable=True)
    rl = RipsLayer(maximum_edge_length=2., dimension=0)

    with tf.GradientTape() as tape:
        dgm = rl.call(X)
        loss = tf.math.reduce_sum(tf.square(.5*(dgm[:,1]-dgm[:,0])))
    grads = tape.gradient(loss, [X])

Documentation for RipsLayer
---------------------------

.. autoclass:: gudhi.tensorflow.RipsLayer
   :members:
   :special-members: __init__
   :show-inheritance:
