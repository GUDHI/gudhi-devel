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

This example computes the first landscape associated to a persistence diagram with four points. The landscape is evaluated evenly on [-1,11].

.. testcode::

    import numpy as np
    import tensorflow as tf
    from gudhi.representations.preprocessing import Padding
    from gudhi.representations import Landscape, PersistenceImage
    from gudhi.perslay import perslay_channel
    import matplotlib.pyplot as plt
    # A single diagram with 4 points
    diags = [np.array([[0.,4.],[1.,2.],[3.,8.],[6.,8.]])]
    diags = Padding(use=True).fit_transform(diags)
    D = np.stack(diags, 0)
    tf.reset_default_graph()
    diagram = tf.placeholder(tf.float32, shape=D.shape)
    feed = {diagram: D}
    list_v = []
    samples = np.array(np.arange(-1.,11.,.5), dtype=np.float32)
    perslay_channel(output=list_v, name="perslay", diag=diagram, persistence_weight=None, 
                    perm_op="topk", keep=1, layer="ls", num_samples=len(samples),
                    sample_init=samples, sample_const=True, fc_layers=[])
    vector = tf.concat(list_v, 1)
    init = tf.global_variables_initializer()
    with tf.Session() as sess:
        sess.run(init)
        L = vector.eval(feed_dict=feed)[0,:]
        print(L)

The output is:

.. testoutput::

    [0.         0.         0.         0.70710677 1.4142135  2.1213202
     2.828427   2.1213202  1.4142135  0.70710677 1.4142135  2.1213202
     2.828427   3.535534   2.828427   2.1213202  1.4142135  0.70710677
     0.         0.         0.         0.         0.         0.        ]

