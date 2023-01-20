:orphan:

.. To get rid of WARNING: document isn't included in any toctree

====================================
TensorFlow layer for representations
====================================

.. include:: representations_sum.inc

`PersLay <http://proceedings.mlr.press/v108/carriere20a.html>`_ is a layer for neural network architectures that allows
to automatically learn the best representation to use for persistence diagrams in supervised machine learning during
training time. Its parameters allow to reproduce most of the known finite-dimensional representations (such as, e.g.,
landscapes and images), and can be combined to create even new ones, that are best suited for a given supervised
machine learning task with persistence diagrams. PersLay is implemented in TensorFlow.

This `notebook <https://github.com/GUDHI/TDA-tutorial/blob/master/Tuto-GUDHI-perslay-visu.ipynb>`__ explains how to use
PersLay.

Example
-------

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


Perslay reference manual
------------------------

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
