:orphan:

.. To get rid of WARNING: document isn't included in any toctree

==============
PersLay manual
==============

.. include:: perslay_sum.inc

This module, originally available at https://github.com/MathieuCarriere/perslay, aims at bridging the gap between persistence diagrams and deep learning, by providing the implementation in tensorflow of a neural network layer for persistence diagrams, called PersLay, defined in https://arxiv.org/abs/1904.09378. More specifically, it provides code for reproducing most of the persistence diagram representations in the literature, and optimize their parameters with a neural network.

Since tensorflow requires its inputs to be fixed-size matrices, persistence diagrams are padded  to a fixed number of points, with a third boolean coordinate for each points assessing whether the point is real or just a dummy point added after padding. A diagram is thus represented as a tensorflow tensor of shape (n,3), as can be obtained from :func:`~gudhi.SimplexTree.persistence_intervals_in_dimension` and :func:`~gudhi.representations.preprocessing.Padding` for instance. Points at infinity are represented as a tensorflow tensors of shape (n,2), storing only the birth times.

A small example is provided

.. only:: builder_html

    * :download:`perslay_visu.py <../example/perslay_visu.py>`

Function
--------
.. autofunction:: gudhi.perslay.perslay_channel

Basic example
-------------

This example computes the first two Landscapes associated to a persistence diagram with four points. The landscapes are evaluated on ten samples, leading to two vectors with ten coordinates each, that are eventually concatenated in order to produce a single vector representation.

.. testcode::

    import numpy as np
    import tensorflow as tf
    from gudhi.representations.preprocessing import Padding 
    from gudhi.perslay import perslay_channel
    # A single diagram with 4 points
    diags = [np.array([[0.,4.],[1.,2.],[3.,8.],[6.,8.]])]
    diags = Padding(use=True).fit_transform(diags)
    D = np.stack(diags, 0)
    tf.reset_default_graph()
    diagram = tf.placeholder(tf.float32, shape=D.shape)
    feed = {diagram: D}
    list_v = []
    perslay_channel(output=list_v, name="perslay", diag=diagram, "persistence_weight"=None, "perm_op"="topk", "keep"=3, "layer"="ls", "num_samples"=3000, \
                   "sample_init"=np.array([[ np.arange(-1.,2.,.001) ]], dtype=np.float32), "sample_const"= True, "fc_layers"=[])
    vector = tf.concat(list_v, 1)
    init = tf.global_variables_initializer()
    with tf.Session() as sess:
        sess.run(init)
        L = vector.eval(feed_dict=feed)[0,:]
        print(V)

The output is:

.. testoutput::

    [[1.02851895 2.05703791 2.57129739 1.54277843 0.89995409 1.92847304
      2.95699199 3.08555686 2.05703791 1.02851895 0.         0.64282435
      0.         0.         0.51425948 0.         0.         0.
      0.77138922 1.02851895]]
