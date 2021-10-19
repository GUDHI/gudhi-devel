:orphan:

.. To get rid of WARNING: document isn't included in any toctree

TensorFlow layer for cubical persistence
########################################

.. include:: differentiation_sum.inc

Example of gradient computed from cubical persistence
-----------------------------------------------------

.. code-block:: python

    from gudhi.tensorflow import *
    import numpy as np
    import tensorflow as tf

    Xinit = np.array([[0.,2.,2.],[2.,2.,2.],[2.,2.,1.]], dtype=np.float32)
    X = tf.Variable(initial_value=Xinit, trainable=True)
    cl = CubicalLayer(dimension=0)

    with tf.GradientTape() as tape:
        dgm = cl.call(X)
        loss = tf.math.reduce_sum(tf.square(.5*(dgm[:,1]-dgm[:,0])))
    grads = tape.gradient(loss, [X])
    print(grads[0].numpy())

Documentation for CubicalLayer
------------------------------

.. autoclass:: gudhi.tensorflow.CubicalLayer
   :members:
   :special-members: __init__
   :show-inheritance:
