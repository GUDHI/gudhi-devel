from gudhi.differentiation import *
import numpy as np
import tensorflow as tf
import gudhi as gd

def test_rips_diff():

    Xinit = np.array([[1.,1.],[2.,2.]], dtype=np.float32)
    X = tf.Variable(initial_value=Xinit, trainable=True)
    rl = RipsLayer(maximum_edge_length=2., dimension=0)

    with tf.GradientTape() as tape:
        dgm = rl.call(X)
        loss = tf.math.reduce_sum(tf.square(.5*(dgm[:,1]-dgm[:,0])))
    grads = tape.gradient(loss, [X])
    assert np.abs(grads[0].numpy()-np.array([[-.5,-.5],[.5,.5]])).sum() <= 1e-6


def test_cubical_diff():

    Xinit = np.array([[0.,2.,2.],[2.,2.,2.],[2.,2.,1.]], dtype=np.float32)
    X = tf.Variable(initial_value=Xinit, trainable=True)
    cl = CubicalLayer(dimension=0)

    with tf.GradientTape() as tape:
        dgm = cl.call(X)
        loss = tf.math.reduce_sum(tf.square(.5*(dgm[:,1]-dgm[:,0])))
    grads = tape.gradient(loss, [X])
    assert np.abs(grads[0].numpy()-np.array([[0.,0.,0.],[0.,.5,0.],[0.,0.,-.5]])).sum() <= 1e-6

def test_st_diff():

    st = gd.SimplexTree()
    st.insert([0])
    st.insert([1]) 
    st.insert([2]) 
    st.insert([3]) 
    st.insert([4]) 
    st.insert([5]) 
    st.insert([6]) 
    st.insert([7]) 
    st.insert([8]) 
    st.insert([9]) 
    st.insert([10]) 
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

    Finit = np.array([6.,4.,3.,4.,5.,4.,3.,2.,3.,4.,5.], dtype=np.float32)
    F = tf.Variable(initial_value=Finit, trainable=True)
    sl = LowerStarSimplexTreeLayer(simplextree=st, dimension=0)

    with tf.GradientTape() as tape:
        dgm = sl.call(F)
        loss = tf.math.reduce_sum(tf.square(.5*(dgm[:,1]-dgm[:,0])))
    grads = tape.gradient(loss, [F])

    assert np.array_equal(np.array(grads[0].indices), np.array([2,4])) 
    assert np.array_equal(np.array(grads[0].values), np.array([-1,1]))

