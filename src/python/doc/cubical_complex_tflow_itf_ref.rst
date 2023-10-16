:orphan:

.. To get rid of WARNING: document isn't included in any toctree

TensorFlow layer for cubical persistence
########################################

.. include:: differentiation_sum.inc

Example of gradient computed from cubical persistence
-----------------------------------------------------

.. testcode::

    from gudhi.tensorflow import CubicalLayer
    import tensorflow as tf

    X = tf.Variable([[0.,2.,2.],[2.,2.,2.],[2.,2.,1.]], dtype=tf.float32, trainable=True)
    cl = CubicalLayer(homology_dimensions=[0])

    with tf.GradientTape() as tape:
        dgm = cl.call(X)[0][0]
        loss = tf.math.reduce_sum(tf.square(.5*(dgm[:,1]-dgm[:,0])))

    grads = tape.gradient(loss, [X])
    print(grads[0].numpy())

.. testoutput::

    [[ 0.   0.   0. ]
     [ 0.   0.5  0. ]
     [ 0.   0.  -0.5]]

Documentation for CubicalLayer
------------------------------

.. autoclass:: gudhi.tensorflow.CubicalLayer
   :members:
   :show-inheritance:
