:orphan:

.. To get rid of WARNING: document isn't included in any toctree

======================
Representations manual
======================

.. include:: representations_sum.inc

This module aims at bridging the gap between persistence diagrams and machine learning, by providing implementations of most of the vector representations for persistence diagrams in the literature, in a scikit-learn format. More specifically, it provides tools, using the scikit-learn standard interface, to compute distances and kernels on persistence diagrams, and to convert these diagrams into vectors in Euclidean space. Moreover, this module also contains `PersLay <http://proceedings.mlr.press/v108/carriere20a.html>`_, which is a general neural network layer for performing deep learning with persistence diagrams, implemented in TensorFlow. 

A diagram is represented as a numpy array of shape (n,2), as can be obtained from :func:`~gudhi.SimplexTree.persistence_intervals_in_dimension` for instance. Points at infinity are represented as a numpy array of shape (n,1), storing only the birth time. The classes in this module can handle several persistence diagrams at once. In that case, the diagrams are provided as a list of numpy arrays. Note that it is not necessary for the diagrams to have the same number of points, i.e., for the corresponding arrays to have the same number of rows: all classes can handle arrays with different shapes.

This `notebook <https://github.com/GUDHI/TDA-tutorial/blob/master/Tuto-GUDHI-representations.ipynb>`_ explains how to
efficiently combine machine learning and topological data analysis with the
:doc:`representations module<representations>` in a scikit-learn fashion. This `notebook <https://github.com/MathieuCarriere/tda-tutorials/blob/perslay/Tuto-GUDHI-perslay-expe.ipynb>`_ 
and `this one <https://github.com/MathieuCarriere/tda-tutorials/blob/perslay/Tuto-GUDHI-perslay-visu.ipynb>`_ explain how to use PersLay.


Examples
--------

Landscapes
^^^^^^^^^^

This example computes the first two Landscapes associated to a persistence diagram with four points. The landscapes are evaluated on ten samples, leading to two vectors with ten coordinates each, that are eventually concatenated in order to produce a single vector representation.

.. testcode::

    import numpy as np
    from gudhi.representations import Landscape
    # A single diagram with 4 points
    D = np.array([[0.,4.],[1.,2.],[3.,8.],[6.,8.]])
    diags = [D]
    l=Landscape(num_landscapes=2,resolution=10).fit_transform(diags)
    print(l) 

.. testoutput::

    [[1.02851895 2.05703791 2.57129739 1.54277843 0.89995409 1.92847304
      2.95699199 3.08555686 2.05703791 1.02851895 0.         0.64282435
      0.         0.         0.51425948 0.         0.         0.
      0.77138922 1.02851895]]

Various kernels
^^^^^^^^^^^^^^^

This small example is also provided
:download:`diagram_vectorizations_distances_kernels.py <../example/diagram_vectorizations_distances_kernels.py>`

PersLay
^^^^^^^

.. testcode::

    import numpy             as np
    import tensorflow        as tf
    from sklearn.preprocessing import MinMaxScaler
    import gudhi.representations as gdr
    import gudhi.tensorflow as gdtf

    diagrams = [np.array([[0.,4.],[1.,2.],[3.,8.],[6.,8.]])]
    diagrams = gdr.DiagramScaler(use=True, scalers=[([0,1], MinMaxScaler())]).fit_transform(diagrams)
    diagrams = tf.RaggedTensor.from_tensor(tf.constant(diagrams, dtype=tf.float32))

    rho = tf.identity 
    phi = gdtf.GaussianPerslayPhi((100, 100), ((-.5, 1.5), (-.5, 1.5)), .1)
    weight = gdtf.PowerPerslayWeight(1.,0.)
    perm_op = tf.math.reduce_sum
    
    perslay = gdtf.Perslay(phi=phi, weight=weight, perm_op=perm_op, rho=rho)
    vectors = perslay(diagrams)
    print(vectors)

.. testoutput::

    tf.Tensor(
    [[[[1.7266072e-16]
       [4.1706043e-09]
       [1.1336876e-08]
       [8.5738821e-12]
       [2.1243891e-14]]

      [[4.1715076e-09]
       [1.0074080e-01]
       [2.7384272e-01]
       [3.0724244e-02]
       [7.6157507e-05]]

      [[8.0382870e-06]
       [1.5802664e+00]
       [8.2997030e-01]
       [1.2395413e+01]
       [3.0724116e-02]]

      [[8.0269419e-06]
       [1.3065740e+00]
       [9.0923014e+00]
       [6.1664842e-02]
       [1.3949171e-06]]

      [[9.0331329e-13]
       [1.4954816e-07]
       [1.5145997e-04]
       [1.0205092e-06]
       [7.8093526e-16]]]], shape=(1, 5, 5, 1), dtype=float32)

Preprocessing
-------------
.. automodule:: gudhi.representations.preprocessing
   :members:
   :special-members:
   :show-inheritance:

Vector methods
--------------
.. automodule:: gudhi.representations.vector_methods
   :members:
   :special-members:
   :show-inheritance:

Kernel methods
--------------
.. automodule:: gudhi.representations.kernel_methods
   :members:
   :special-members:
   :show-inheritance:

Metrics
-------
.. automodule:: gudhi.representations.metrics
   :members:
   :special-members:
   :show-inheritance:

PersLay
-------
.. automodule:: gudhi.tensorflow.perslay
   :members:
   :special-members:
   :show-inheritance:

