from gudhi.differentiation import *
import numpy as np
import tensorflow as tf

def test_rips_diff():

    Xinit = np.array([[1.,1.],[2.,2.]], dtype=np.float32)
    X = tf.Variable(initial_value=Xinit, trainable=True)
    model = RipsModel(X=X, mel=2., dim=0, card=10)

    with tf.GradientTape() as tape:
        dgm = model.call()
        loss = tf.math.reduce_sum(tf.square(.5*(dgm[:,1]-dgm[:,0])))
        grads = tape.gradient(loss, [X])
        assert np.abs(grads[0].numpy()-np.array([[-.5,-.5],[.5,.5]])).sum() <= 1e-6


def test_cubical_diff():

    Xinit = np.array([[0.,2.,2.],[2.,2.,2.],[2.,2.,1.]], dtype=np.float32)
    X = tf.Variable(initial_value=Xinit, trainable=True)
    model = CubicalModel(X, dim=0, card=10)

    with tf.GradientTape() as tape:
        dgm = model.call()
        loss = tf.math.reduce_sum(tf.square(.5*(dgm[:,1]-dgm[:,0])))
        grads = tape.gradient(loss, [X])
        assert np.abs(grads[0].numpy()-np.array([[0.,0.,0.],[0.,.5,0.],[0.,0.,-.5]])).sum() <= 1e-6

def test_st_diff():

    Finit = np.array([6.,4.,3.,4.,5.,4.,3.,2.,3.,4.,5.], dtype=np.float32)
    F = tf.Variable(initial_value=Finit, trainable=True)
    model = SimplexTreeModel(F, stbase="simplextree.txt", dim=0, card=10)

    with tf.GradientTape() as tape:
        dgm = model.call()
        loss = tf.math.reduce_sum(tf.square(.5*(dgm[:,1]-dgm[:,0])))
        grads = tape.gradient(loss, [F])
        assert np.array_equal(np.array(grads[0].indices), np.array([2,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])) 
        assert np.array_equal(np.array(grads[0].values), np.array([-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]))
