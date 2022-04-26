:orphan:

.. To get rid of WARNING: document isn't included in any toctree

RipsNet user manual
=========================
Definition
----------

.. include:: ripsnet.inc

:class:`~gudhi.ripsnet` constructs a Tensorflow model for fast and robust estimation of persistent homology of
point clouds.
RipsNet is based on a Deep Sets architecture :cite:`DeepSets17`, for details see the paper RipsNet :cite:`RipsNet_arXiv`.

Example
-------------------

This example instantiates a RipsNet model which can then be trained as any tensorflow model.

.. testcode::

    from gudhi.tensorflow import *
    import tensorflow as tf
    from tensorflow.keras import regularizers, layers

    ragged_layers_size = [20, 10]
    dense_layers_size = [10, 20]
    output_units = 25
    activation_fct = 'gelu'
    output_activation = 'sigmoid'
    dropout = 0
    kernel_regularization = 0

    ragged_layers = []
    dense_layers = []

    for n_units in ragged_layers_size:
        ragged_layers.append(DenseRagged(units=n_units, use_bias=True, activation=activation_fct))

    for n_units in dense_layers_size:
        dense_layers.append(layers.Dense(n_units, activation=activation_fct,
                                         kernel_regularizer=regularizers.l2(kernel_regularization)))
        dense_layers.append(layers.Dropout(dropout))

    dense_layers.append(layers.Dense(output_units, activation=output_activation))

    phi_1 = DenseRaggedBlock(ragged_layers)
    perm_op = 'mean'  # can also be 'sum'.
    phi_2 = TFBlock(dense_layers)
    input_dim = 2

    RN = RipsNet(phi_1, phi_2, input_dim, perm_op=perm_op)

    data_test = [[[-7.04493841, 9.60285858],
                  [-13.14389003, -13.21854157],
                  [-3.21137961, -1.28593644]],
                 [[10.40324933, -0.80540584],
                  [16.54752459, 0.70355361],
                  [6.410207, -10.63175183],
                  [2.96613799, -11.97463568]],
                 [[4.85041719, -2.93820024],
                  [2.15379915, -5.39669696],
                  [5.83968556, -5.67350982],
                  [5.25955172, -6.36860269]]]

    tf_data_test = tf.ragged.constant(data_test, ragged_rank=1)

    RN.predict(tf_data_test)

Once RN is properly trained (which we skip in this documentation) it can be used to make predictions.
In this example RipsNet estimates persistence vectorizations (of output size 25) of a list of point clouds (of 3 points) in 2D.
A possible output is:

.. code-block::

    [[0.58554363 0.6054868  0.44672886 0.5216672  0.5814481  0.48068565
      0.49626726 0.5285395  0.4805212  0.37918684 0.49745193 0.49247316
      0.4706078  0.5491477  0.47016636 0.55804974 0.46501246 0.4065692
      0.5386659  0.5660226  0.52014357 0.5329493  0.52178216 0.5156043
      0.48742113]
     [0.9446074  0.99024785 0.1316272  0.3013248  0.98174655 0.52285945
      0.33727515 0.997285   0.3711884  0.00388432 0.63181967 0.5377489
      0.22074646 0.7681194  0.04337704 0.80116796 0.02139336 0.04605395
      0.8911999  0.9570045  0.5789719  0.8221929  0.7742506  0.4596561
      0.08529088]
     [0.8230771  0.9320036  0.25120026 0.48027694 0.8988322  0.5789062
      0.38307947 0.9252455  0.39485127 0.06090912 0.5786307  0.51115406
      0.28706372 0.70552015 0.16929033 0.7028084  0.12379596 0.1867683
      0.6969584  0.84437454 0.6172329  0.66728634 0.630455   0.47643042
      0.27172992]]

Detailed documentation
----------------------
.. automodule:: gudhi.tensorflow.ripsnet
   :members:
   :special-members:
   :show-inheritance:
