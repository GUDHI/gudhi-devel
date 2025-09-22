""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Mathieu Carriere

    Copyright (C) 2020 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""


import numpy as np
import torch
import tensorflow as tf

from gudhi.wasserstein import wasserstein_distance as pot


def test_wasserstein_distance_grad():
    diag1 = torch.tensor([[2.7, 3.7], [9.6, 14.0], [34.2, 34.974]], requires_grad=True)
    diag2 = torch.tensor([[2.8, 4.45], [9.5, 14.1]], requires_grad=True)
    diag3 = torch.tensor([[2.8, 4.45], [9.5, 14.1]], requires_grad=True)
    assert diag1.grad is None and diag2.grad is None and diag3.grad is None
    dist12 = pot(diag1, diag2, internal_p=2, order=2, enable_autodiff=True, keep_essential_parts=False)
    dist30 = pot(diag3, torch.tensor([]), internal_p=2, order=2, enable_autodiff=True, keep_essential_parts=False)
    dist12.backward()
    dist30.backward()
    assert (
        not torch.isnan(diag1.grad).any()
        and not torch.isnan(diag2.grad).any()
        and not torch.isnan(diag3.grad).any()
    )
    diag4 = torch.tensor([[0.0, 10.0]], requires_grad=True)
    diag5 = torch.tensor([[1.0, 11.0], [3.0, 4.0]], requires_grad=True)
    dist45 = pot(diag4, diag5, internal_p=1, order=1, enable_autodiff=True, keep_essential_parts=False)
    assert dist45 == 3.0
    dist45.backward()
    assert np.array_equal(diag4.grad, [[-1.0, -1.0]])
    assert np.array_equal(diag5.grad, [[1.0, 1.0], [-1.0, 1.0]])
    diag6 = torch.tensor([[5.0, 10.0]], requires_grad=True)
    pot(diag6, diag6, internal_p=2, order=2, enable_autodiff=True, keep_essential_parts=False).backward()
    # https://github.com/jonasrauber/eagerpy/issues/6
    # assert np.array_equal(diag6.grad, [[0., 0.]])


def test_wasserstein_distance_grad_tensorflow():
    with tf.GradientTape() as tape:
        diag4 = tf.convert_to_tensor(
            tf.Variable(initial_value=np.array([[0.0, 10.0]]), trainable=True)
        )
        diag5 = tf.convert_to_tensor(
            tf.Variable(initial_value=np.array([[1.0, 11.0], [3.0, 4.0]]), trainable=True)
        )
        dist45 = pot(diag4, diag5, internal_p=1, order=1, enable_autodiff=True, keep_essential_parts=False)
        assert dist45 == 3.0

    grads = tape.gradient(dist45, [diag4, diag5])
    assert np.array_equal(grads[0].values, [[-1.0, -1.0]])
    assert np.array_equal(grads[1].values, [[1.0, 1.0], [-1.0, 1.0]])
