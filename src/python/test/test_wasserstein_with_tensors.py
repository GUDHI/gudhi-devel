""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Mathieu Carriere

    Copyright (C) 2020 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

from gudhi.wasserstein import wasserstein_distance as pot
import numpy as np

def test_wasserstein_distance_grad_tensorflow():
    import tensorflow as tf

    with tf.GradientTape() as tape:
        diag4 = tf.convert_to_tensor(tf.Variable(initial_value=np.array([[0., 10.]]),           trainable=True))
        diag5 = tf.convert_to_tensor(tf.Variable(initial_value=np.array([[1., 11.], [3., 4.]]), trainable=True))
        dist45 = pot(diag4, diag5, internal_p=1, order=1, enable_autodiff=True)
        assert dist45 == 3.

    grads = tape.gradient(dist45, [diag4, diag5])
    assert np.array_equal(grads[0].values, [[-1., -1.]])
    assert np.array_equal(grads[1].values, [[1., 1.], [-1., 1.]])