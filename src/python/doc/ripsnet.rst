:orphan:

.. To get rid of WARNING: document isn't included in any toctree

RipsNet user manual
=========================
Definition
----------

.. include:: ripsnet.inc

:class:`~gudhi.ripsnet` constructs a Tensorflow model for fast and robust estimation of persistent homology of
point clouds.
RipsNet is based on a `Deep Sets  <https://papers.nips.cc/paper/2017/file/f22e4747da1aa27e363d86d40ff442fe-Paper.pdf>`_
architecture, for details see `RipsNet  <https://arxiv.org/abs/2202.01725>`_.

Example
-------------------

This example instantiates a RipsNet model which can then be trained as any tensorflow model.

.. testcode::
    from gudhi import ripsnet
    from tensorflow.keras import regularizers, layers

    ragged_layers_size    = [30,20,10]
    dense_layers_size     = [50,100,200]
    output_units          = 2500
    activation_fct        = 'gelu'
    output_activation     = 'sigmoid'
    dropout               = 0
    kernel_regularization = 0

    ragged_layers = []
    dense_layers = []

    for n_units in ragged_layers_size:
            ragged_layers.append(DenseRagged(units=n_units, use_bias=True, activation=activation_fct))

    for n_units in dense_layers_size:
        dense_layers.append(tf.keras.layers.Dense(n_units, activation=activation_fct,
                                kernel_regularizer=regularizers.l2(kernel_regularization)))
        dense_layers.append(tf.keras.layers.Dropout(dropout))

    dense_layers.append(tf.keras.layers.Dense(output_units, activation=output_activation))

    phi_1     = DenseRaggedBlock(ragged_layers)
    perm_op   = 'mean' # can also be 'sum'.
    phi_2     = TFBlock(dense_layers)
    input_dim = 2

    RN = RipsNet(phi_1, phi_2, input_dim, perm_op=perm_op)
