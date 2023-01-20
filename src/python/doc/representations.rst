:orphan:

.. To get rid of WARNING: document isn't included in any toctree

======================
Representations manual
======================

.. include:: representations_sum.inc

This module aims at bridging the gap between persistence diagrams and machine learning, by providing implementations of most of the vector representations for persistence diagrams in the literature, in a scikit-learn format. More specifically, it provides tools, using the scikit-learn standard interface, to compute distances and kernels on persistence diagrams, and to convert these diagrams into vectors in Euclidean space. Moreover, this module also contains `PersLay <http://proceedings.mlr.press/v108/carriere20a.html>`_, which is a general neural network layer for performing deep learning with persistence diagrams, implemented in TensorFlow.

A diagram is represented as a numpy array of shape (n,2), as can be obtained from :func:`~gudhi.SimplexTree.persistence_intervals_in_dimension` for instance. Points at infinity are represented as a numpy array of shape (n,1), storing only the birth time. The classes in this module can handle several persistence diagrams at once. In that case, the diagrams are provided as a list of numpy arrays. Note that it is not necessary for the diagrams to have the same number of points, i.e., for the corresponding arrays to have the same number of rows: all classes can handle arrays with different shapes.

This `notebook <https://github.com/GUDHI/TDA-tutorial/blob/master/Tuto-GUDHI-representations.ipynb>`__ explains how to
efficiently combine machine learning and topological data analysis with the
:doc:`representations module<representations>` in a scikit-learn fashion. This `notebook <https://github.com/MathieuCarriere/tda-tutorials/blob/perslay/Tuto-GUDHI-perslay-expe.ipynb>`__
and `this one <https://github.com/MathieuCarriere/tda-tutorials/blob/perslay/Tuto-GUDHI-perslay-visu.ipynb>`__ explain how to use PersLay.


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

.. testsetup:: perslay

    import numpy
    numpy.set_printoptions(precision=5)

.. testcode:: perslay

    import numpy             as np
    import tensorflow        as tf
    from sklearn.preprocessing import MinMaxScaler
    import gudhi.representations as gdr
    import gudhi.tensorflow.perslay as prsl

    diagrams = [np.array([[0.,4.],[1.,2.],[3.,8.],[6.,8.]])]
    diagrams = gdr.DiagramScaler(use=True, scalers=[([0,1], MinMaxScaler())]).fit_transform(diagrams)
    diagrams = tf.RaggedTensor.from_tensor(tf.constant(diagrams, dtype=tf.float32))

    rho = tf.identity
    phi = prsl.GaussianPerslayPhi((5, 5), ((-.5, 1.5), (-.5, 1.5)), .1)
    weight = prsl.PowerPerslayWeight(1.,0.)
    perm_op = tf.math.reduce_sum

    perslay = prsl.Perslay(phi=phi, weight=weight, perm_op=perm_op, rho=rho)
    vectors = perslay(diagrams)
    print(vectors)

.. testcleanup:: perslay

    numpy.set_printoptions(precision=8)

.. testoutput:: perslay

    tf.Tensor(
    [[[[1.72661e-16]
       [4.17060e-09]
       [1.13369e-08]
       [8.57388e-12]
       [2.12439e-14]]

      [[4.17151e-09]
       [1.00741e-01]
       [2.73843e-01]
       [3.07242e-02]
       [7.61575e-05]]

      [[8.03829e-06]
       [1.58027e+00]
       [8.29970e-01]
       [1.23954e+01]
       [3.07241e-02]]

      [[8.02694e-06]
       [1.30657e+00]
       [9.09230e+00]
       [6.16648e-02]
       [1.39492e-06]]

      [[9.03313e-13]
       [1.49548e-07]
       [1.51460e-04]
       [1.02051e-06]
       [7.80935e-16]]]], shape=(1, 5, 5, 1), dtype=float32)

Utils
-----
.. autoclass:: gudhi.representations.preprocessing.Clamping
   :members:
   :special-members:
   :show-inheritance:

Preprocessing
-------------
.. automodule:: gudhi.representations.preprocessing
   :members:
   :special-members:
   :show-inheritance:
   :exclude-members: Clamping

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
.. autoclass:: gudhi.tensorflow.perslay.Perslay
   :members:
   :special-members:
   :show-inheritance:

Weight functions
^^^^^^^^^^^^^^^^
.. autoclass:: gudhi.tensorflow.perslay.GaussianMixturePerslayWeight
   :members:
   :special-members:
   :show-inheritance:

.. autoclass:: gudhi.tensorflow.perslay.GridPerslayWeight
   :members:
   :special-members:
   :show-inheritance:

.. autoclass:: gudhi.tensorflow.perslay.PowerPerslayWeight
   :members:
   :special-members:
   :show-inheritance:

Phi functions
^^^^^^^^^^^^^
.. autoclass:: gudhi.tensorflow.perslay.FlatPerslayPhi
   :members:
   :special-members:
   :show-inheritance:

.. autoclass:: gudhi.tensorflow.perslay.GaussianPerslayPhi
   :members:
   :special-members:
   :show-inheritance:

.. autoclass:: gudhi.tensorflow.perslay.TentPerslayPhi
   :members:
   :special-members:
   :show-inheritance:
