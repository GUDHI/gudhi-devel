:orphan:

.. To get rid of WARNING: document isn't included in any toctree

TensorFlow layer for computing Rips persistence
###############################################

.. list-table::
   :widths: 40 30 30
   :header-rows: 0

   * - :Since: GUDHI 3.5.0
     - :License: MIT
     - :Requires: `TensorFlow <installation.html#tensorflow>`_

TensorFlow layer for computing Rips persistence example
-------------------------------------------------------

.. code-block:: python

    from gudhi.differentiation import *
    import numpy as np
    import tensorflow as tf
    import gudhi as gd

    Xinit = np.array([[1.,1.],[2.,2.]], dtype=np.float32)
    X = tf.Variable(initial_value=Xinit, trainable=True)
    rl = RipsLayer(maximum_edge_length=2., dimension=0)

    with tf.GradientTape() as tape:
        dgm = rl.call(X)
        loss = tf.math.reduce_sum(tf.square(.5*(dgm[:,1]-dgm[:,0])))
    grads = tape.gradient(loss, [X])
    assert np.abs(grads[0].numpy()-np.array([[-.5,-.5],[.5,.5]])).sum() <= 1e-6


TensorFlow layer for computing Rips persistence reference
---------------------------------------------------------

.. autoclass:: gudhi.differentiation.RipsLayer
   :members:
   :special-members: __init__
   :show-inheritance:
