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
    input_dim = 2
    ragged_layers_size = [10]
    dense_layers_size = [10]
    output_units = 2
    activation_fct = 'gelu'
    output_activation = 'sigmoid'
    dropout = 0
    kernel_regularization = 0

    ragged_layers = []
    dense_layers = []
    for n_units in ragged_layers_size:
        ragged_layers.append(DenseRagged(units=n_units, use_bias=True, activation=activation_fct))

    for n_units in dense_layers_size:
        dense_layers.append(tf.keras.layers.Dense(n_units, activation=activation_fct))

    dense_layers.append(tf.keras.layers.Dense(output_units, activation=output_activation))

    weights_vect = [np.array([[-0.82675004, 0.082853, -0.7522571, -0.86455333],
                              [0.01191106, 0.4168273, 0.13942415, -0.95728004]]),
                    np.array([-0.48964924, 0.07440309, -0.84914917, -0.13673241]),
                    np.array([[0.6735323, -0.86710656, -0.6601294, 0.7874811],
                              [-0.5002684, -0.36409453, 0.5563092, -0.16215993],
                              [0.31057122, -0.79817116, -0.40550056, 0.8173016],
                              [0.64593965, 0.5079483, -0.5176215, 0.5098056]]),
                    np.array([0.68127185, -0.03681624, -0.46403527, 0.42004192]),
                    np.array([[-0.23722765, 0.18152483, 0.42457554, 0.55383813],
                              [0.11030384, -0.72890997, 0.81839615, -0.26585487],
                              [-0.16284865, -0.22054555, -0.03428887, -0.10710541],
                              [-0.6388146, -0.6412578, -0.05625373, 0.6325129]]),
                    np.array([0.31346196, 0.32782924, 0.18332976, -0.4718462])]

    phi_1 = DenseRaggedBlock(ragged_layers)
    perm_op = 'mean'
    phi_2 = TFBlock(dense_layers)
    input_dim = 2

    model = RipsNet(phi_1, phi_2, input_dim, perm_op=perm_op)

    model.set_weights(weights_vect)

    clean_data_test = [np.array([[  8.00690107,   5.84065138],
                                 [-10.01871508,   5.34973559],
                                 [ -9.65317915,   5.53776047],
                                 [  7.78598811,   7.76091515]]),
                       np.array([[  9.95347807,  13.72474594],
                                 [ -3.40168368, -15.23953774],
                                 [ -8.73604688,   6.76905796]]),
                       np.array([[ 3.25849098, -4.01998912],
                                 [ 2.23853692, -2.96113938],
                                 [ 3.38496086, -4.47234084],
                                 [ 2.13325439, -2.91673025]])]

    noisy_data_test = [np.array([[-15.49749451,  -4.81412853],
                                 [-10.02162025,   0.42705654],
                                 [-14.41883655,  10.24555917],
                                 [ -8.6865668 ,   6.12261606]]),
                       np.array([[-3.14624955,  2.59008951],
                                 [-3.88721637,  0.65461877],
                                 [-3.16411703,  0.44306022],
                                 [-2.75523434,  1.57533588]]),
                       np.array([[ 10.14006214,  10.92637587],
                                 [-13.75312184,  -2.09110089],
                                 [ -7.83379779,   0.76465675]])]

    tf_clean_data_test = tf.ragged.constant([
        [list(c) for c in list(clean_data_test[i])] for i in range(len(clean_data_test))], ragged_rank=1)
    tf_noisy_data_test = tf.ragged.constant([
        [list(c) for c in list(noisy_data_test[i])] for i in range(len(noisy_data_test))], ragged_rank=1)

    clean_prediction = np.array([[0.3942722,  0.95260733, 0.9447353,  0.00872418],
                        [0.70753145, 0.74972904, 0.45573318, 0.25137717],
                        [0.6459818,  0.5658171,  0.19202137, 0.6256107 ]])
    noisy_prediction = np.array([[0.5073361,  0.93115664, 0.89942205, 0.01895535],
                        [0.4248793,  0.8612615,  0.7404877,  0.09027195],
                        [0.36387908, 0.9295354,  0.90001035, 0.01976481]])

    assert(clean_prediction.shape == model.predict(tf_clean_data_test).shape)
    assert(noisy_prediction.shape == model.predict(tf_noisy_data_test).shape)
    assert((clean_prediction == model.predict(tf_clean_data_test)).all())
    assert((noisy_prediction == model.predict(tf_noisy_data_test)).all())
    return