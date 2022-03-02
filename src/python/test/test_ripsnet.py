""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Felix Hensel
    Copyright (C) 2022 Inria
    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

import numpy             as np
import tensorflow        as tf
from gudhi.tensorflow.ripsnet import *


def test_ripsnet():
    ragged_layers_size = [4]
    dense_layers_size = [4]
    output_units = 4
    activation_fct = 'gelu'
    output_activation = 'sigmoid'
    dropout = 0
    kernel_regularization = 0
    initializer = None# tf.keras.initializers.Constant(value=1)

    ragged_layers = []
    dense_layers = []
    for n_units in ragged_layers_size:
        ragged_layers.append(DenseRagged(units=n_units, use_bias=True, activation=activation_fct,
                                         kernel_initializer=initializer, bias_initializer=initializer))

    for n_units in dense_layers_size:
        dense_layers.append(tf.keras.layers.Dense(n_units, activation=activation_fct,
                                                  kernel_initializer=initializer, bias_initializer=initializer))

    dense_layers.append(tf.keras.layers.Dense(output_units, activation=output_activation,
                                              kernel_initializer=initializer, bias_initializer=initializer))

    weights_vect = [np.array([[-0.3868327 ,  0.5431584 ,  0.7523476 ,  0.80209386],
                              [ 0.22491306,  0.4626178 ,  0.34193814, -0.04737851]]),
                    np.array([ 0.5047069 , -0.11543324, -0.03882882, -0.16129738]),
                    np.array([[-0.7956421  ,  0.2326832  , -0.5405302  ,  0.096256964],
                              [ 0.06973686 ,  0.0251764  , -0.05733281 ,  0.3528394  ],
                              [-0.77462643 ,  0.03330394 , -0.8688136  , -0.22296508 ],
                              [-0.5054477  ,  0.7201048  ,  0.1857564  ,  0.65894866 ]]),
                    np.array([-0.30565566 , -0.77507186 , -0.049963538,  0.5765676  ]),
                    np.array([[-0.25560755  ,  0.71504813  ,  0.0047909063, -0.1595783   ],
                              [-0.71575665  ,  0.6139034   , -0.47060093  ,  0.087501734 ],
                              [ 0.1588738   , -0.593038    ,  0.48378325  , -0.777213    ],
                              [ 0.6206032   , -0.20880768  ,  0.14528894  ,  0.18696047  ]]),
                    np.array([-0.17761804, -0.6905532 ,  0.64367545, -0.2173939 ])]

    phi_1 = DenseRaggedBlock(ragged_layers)
    perm_op = 'mean'
    phi_2 = TFBlock(dense_layers)
    input_dim = 2

    model = RipsNet(phi_1, phi_2, input_dim, perm_op=perm_op)

    test_input_raw = [np.array([[1.,2.],[3.,4.]])]

    test_input = tf.ragged.constant([
        [list(c) for c in list(test_input_raw[i])] for i in range(len(test_input_raw))], ragged_rank=1)

    model.predict(test_input)

    model.set_weights(weights_vect)


    clean_data_test = [np.array([[ -7.04493841,   9.60285858],
                                 [-13.14389003, -13.21854157],
                                 [ -3.21137961,  -1.28593644]]),
                       np.array([[ 10.40324933,  -0.80540584],
                                 [ 16.54752459,   0.70355361],
                                 [  6.410207  , -10.63175183],
                                 [  2.96613799, -11.97463568]]),
                       np.array([[ 4.85041719, -2.93820024],
                                 [ 2.15379915, -5.39669696],
                                 [ 5.83968556, -5.67350982],
                                 [ 5.25955172, -6.36860269]])]

    noisy_data_test = [np.array([[ -8.93311026,   1.52317533],
                                 [-16.80344139,  -3.76871298],
                                 [-11.58448573,  -2.76311122],
                                 [-15.06107796,   5.05253587]]),
                       np.array([[-3.834947  , -5.1897498 ],
                                 [-3.51701182, -4.23539191],
                                 [-2.68678747, -1.63902703],
                                 [-4.65070816, -3.96363227]]),
                       np.array([[ 4.7841113 , 19.2922069 ],
                                 [10.5164214 ,  5.50246605],
                                 [-9.38163622,  7.03682948]])]

    tf_clean_data_test = tf.ragged.constant([
        [list(c) for c in list(clean_data_test[i])] for i in range(len(clean_data_test))], ragged_rank=1)
    tf_noisy_data_test = tf.ragged.constant([
        [list(c) for c in list(noisy_data_test[i])] for i in range(len(noisy_data_test))], ragged_rank=1)

    clean_prediction = np.array([[0.5736222, 0.3047213, 0.6746019, 0.49468565],
                                 [0.4522748, 0.7521156, 0.3061385, 0.77519494],
                                 [0.5349713, 0.49940312, 0.51753736, 0.6435147]])
    noisy_prediction = np.array([[0.53986603, 0.33934325, 0.64809155, 0.49939266],
                                 [0.5637899, 0.28744557, 0.67097586, 0.50240993],
                                 [0.5689339,  0.55574864, 0.47079712, 0.68721807]])

    #print(np.linalg.norm(clean_prediction - model.predict(tf_clean_data_test)))

    assert(clean_prediction.shape == model.predict(tf_clean_data_test).shape)
    assert(noisy_prediction.shape == model.predict(tf_noisy_data_test).shape)
    assert(np.linalg.norm(clean_prediction - model.predict(tf_clean_data_test)) <= 1e-7)
    assert(np.linalg.norm(noisy_prediction - model.predict(tf_noisy_data_test)) <= 1e-7)
    return
