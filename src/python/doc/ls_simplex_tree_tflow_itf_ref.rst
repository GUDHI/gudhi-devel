:orphan:

.. To get rid of WARNING: document isn't included in any toctree

TensorFlow layer for lower-star persistence on simplex trees
############################################################

.. include:: differentiation_sum.inc

Example of gradient computed from lower-star filtration of a simplex tree
-------------------------------------------------------------------------

.. testcode::

    from gudhi.tensorflow import LowerStarSimplexTreeLayer
    import tensorflow as tf
    import gudhi as gd

    st = gd.SimplexTree()
    st.insert([0, 1]) 
    st.insert([1, 2]) 
    st.insert([2, 3]) 
    st.insert([3, 4]) 
    st.insert([4, 5]) 
    st.insert([5, 6]) 
    st.insert([6, 7]) 
    st.insert([7, 8]) 
    st.insert([8, 9]) 
    st.insert([9, 10]) 

    F = tf.Variable([6.,4.,3.,4.,5.,4.,3.,2.,3.,4.,5.], dtype=tf.float32, trainable=True)
    sl = LowerStarSimplexTreeLayer(simplextree=st, homology_dimensions=[0])

    with tf.GradientTape() as tape:
        dgm = sl.call(F)[0][0]
        loss = tf.math.reduce_sum(tf.square(.5*(dgm[:,1]-dgm[:,0])))

    grads = tape.gradient(loss, [F])
    print(grads[0].indices.numpy())
    print(grads[0].values.numpy())

.. testoutput::

    [2 4]
    [-1.  1.]

Documentation for LowerStarSimplexTreeLayer
-------------------------------------------

.. autoclass:: gudhi.tensorflow.LowerStarSimplexTreeLayer
   :members:
   :show-inheritance:
