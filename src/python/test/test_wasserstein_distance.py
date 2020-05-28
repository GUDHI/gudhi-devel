""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Theo Lacombe, Marc Glisse

    Copyright (C) 2019 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

from gudhi.wasserstein.wasserstein import _proj_on_diag
from gudhi.wasserstein import wasserstein_distance as pot
from gudhi.hera import wasserstein_distance as hera
import numpy as np
import pytest

__author__ = "Theo Lacombe"
__copyright__ = "Copyright (C) 2019 Inria"
__license__ = "MIT"


def test_proj_on_diag():
    dgm = np.array([[1., 1.], [1., 2.], [3., 5.]])
    assert np.array_equal(_proj_on_diag(dgm), [[1., 1.], [1.5, 1.5], [4., 4.]])
    empty = np.empty((0, 2))
    assert np.array_equal(_proj_on_diag(empty), empty)


def _basic_wasserstein(wasserstein_distance, delta, test_infinity=True, test_matching=True):
    diag1 = np.array([[2.7, 3.7], [9.6, 14.0], [34.2, 34.974]])
    diag2 = np.array([[2.8, 4.45], [9.5, 14.1]])
    diag3 = np.array([[0, 2], [4, 6]])
    diag4 = np.array([[0, 3], [4, 8]])
    emptydiag = np.array([])

    # We just need to handle positive numbers here
    def approx(x):
        return pytest.approx(x, rel=delta)

    assert wasserstein_distance(emptydiag, emptydiag, internal_p=2.,     order=1.) == 0.
    assert wasserstein_distance(emptydiag, emptydiag, internal_p=np.inf, order=1.) == 0.
    assert wasserstein_distance(emptydiag, emptydiag, internal_p=np.inf, order=2.) == 0.
    assert wasserstein_distance(emptydiag, emptydiag, internal_p=2.,     order=2.) == 0.

    assert wasserstein_distance(diag3, emptydiag, internal_p=np.inf,     order=1.) == approx(2.)
    assert wasserstein_distance(diag3, emptydiag, internal_p=1.,         order=1.) == approx(4.)

    assert wasserstein_distance(diag4, emptydiag, internal_p=1.,     order=2.) == approx(5.)  # thank you Pythagorician triplets
    assert wasserstein_distance(diag4, emptydiag, internal_p=np.inf, order=2.) == approx(2.5)
    assert wasserstein_distance(diag4, emptydiag, internal_p=2.,     order=2.) == approx(3.5355339059327378)

    assert wasserstein_distance(diag1, diag2, internal_p=2.,   order=1.)   == approx(1.4453593023967701)
    assert wasserstein_distance(diag1, diag2, internal_p=2.35, order=1.74) == approx(0.9772734057168739)

    assert wasserstein_distance(diag1, emptydiag, internal_p=2.35, order=1.7863) == approx(3.141592214572228)

    assert wasserstein_distance(diag3, diag4, internal_p=1.,     order=1.) == approx(3.)
    assert wasserstein_distance(diag3, diag4, internal_p=np.inf, order=1.) == approx(3.)  # no diag matching here
    assert wasserstein_distance(diag3, diag4, internal_p=np.inf, order=2.) == approx(np.sqrt(5))
    assert wasserstein_distance(diag3, diag4, internal_p=1.,     order=2.) == approx(np.sqrt(5))
    assert wasserstein_distance(diag3, diag4, internal_p=4.5,    order=2.) == approx(np.sqrt(5))

    if test_infinity:
        diag5 = np.array([[0, 3], [4, np.inf]])
        diag6 = np.array([[7, 8], [4, 6], [3, np.inf]])

        assert wasserstein_distance(diag4, diag5) == np.inf
        assert wasserstein_distance(diag5, diag6, order=1, internal_p=np.inf) == approx(4.)


    if test_matching:
        match = wasserstein_distance(emptydiag, emptydiag, matching=True, internal_p=1., order=2)[1]
        assert np.array_equal(match, [])
        match = wasserstein_distance(emptydiag, emptydiag, matching=True, internal_p=np.inf, order=2.24)[1]
        assert np.array_equal(match, [])
        match = wasserstein_distance(emptydiag, diag2, matching=True, internal_p=np.inf, order=2.)[1]
        assert np.array_equal(match , [[-1, 0], [-1, 1]])
        match = wasserstein_distance(diag2, emptydiag, matching=True, internal_p=np.inf, order=2.24)[1]
        assert np.array_equal(match , [[0, -1], [1, -1]])
        match = wasserstein_distance(diag1, diag2, matching=True, internal_p=2., order=2.)[1]
        assert np.array_equal(match, [[0, 0], [1, 1], [2, -1]])

    if test_entropy:
        assert wasserstein_distance(diag1, diag2, internal_p=2., order=2., reg=1) == approx(1.2251725975696444)
        with pytest.raises(FloatingPointError):
            wasserstein_distance(diag1, diag2, internal_p=2., order=2., reg=1E-20)
        res1 = wasserstein_distance(diag1, diag3, internal_p=2., order=2.)
        res2 = wasserstein_distance(diag1, diag3, internal_p=2., order=2., reg=0)
        assert res1 == res2  # these two should be very much the same.
        with pytest.raises(ValueError):
            wasserstein_distance(diag1, diag2, reg=-1)

        cost_res, P_res = 3.4637979748620684, np.array([[0.46505399, 0.0042368 , 0.53070921],
                                                        [0.03878662, 0.47330083, 0.48791255],
                                                        [0.49615939, 0.52246237, 0.98137825]])
        cost, P = wasserstein_distance(diag3, diag4, internal_p=2., order=2., reg = 10, matching=True)
        assert cost == approx(cost_res)
        assert np.linalg.norm(P - P_res) < 1E-5



def hera_wrap(**extra):
    def fun(*kargs,**kwargs):
        return hera(*kargs,**kwargs,**extra)
    return fun

def pot_wrap(**extra):
    def fun(*kargs,**kwargs):
        return pot(*kargs,**kwargs,**extra)
    return fun

def test_wasserstein_distance_pot():
    _basic_wasserstein(pot, 1e-15, test_infinity=False, test_matching=True)
    _basic_wasserstein(pot_wrap(enable_autodiff=True), 1e-15, test_infinity=False, test_matching=False)

def test_wasserstein_distance_hera():
    _basic_wasserstein(hera_wrap(delta=1e-12), 1e-12, test_matching=False)
    _basic_wasserstein(hera_wrap(delta=.1), .1, test_matching=False)

def test_wasserstein_distance_grad():
    import torch

    diag1 = torch.tensor([[2.7, 3.7], [9.6, 14.0], [34.2, 34.974]], requires_grad=True)
    diag2 = torch.tensor([[2.8, 4.45], [9.5, 14.1]], requires_grad=True)
    diag3 = torch.tensor([[2.8, 4.45], [9.5, 14.1]], requires_grad=True)
    assert diag1.grad is None and diag2.grad is None and diag3.grad is None
    dist12 = pot(diag1, diag2, internal_p=2, order=2, enable_autodiff=True)
    dist30 = pot(diag3, torch.tensor([]), internal_p=2, order=2, enable_autodiff=True)
    dist12.backward()
    dist30.backward()
    assert not torch.isnan(diag1.grad).any() and not torch.isnan(diag2.grad).any() and not torch.isnan(diag3.grad).any()
    diag4 = torch.tensor([[0., 10.]], requires_grad=True)
    diag5 = torch.tensor([[1., 11.], [3., 4.]], requires_grad=True)
    dist45 = pot(diag4, diag5, internal_p=1, order=1, enable_autodiff=True)
    assert dist45 == 3.
    dist45.backward()
    assert np.array_equal(diag4.grad, [[-1., -1.]])
    assert np.array_equal(diag5.grad, [[1., 1.], [-1., 1.]])
    diag6 = torch.tensor([[5., 10.]], requires_grad=True)
    pot(diag6, diag6, internal_p=2, order=2, enable_autodiff=True).backward()
    # https://github.com/jonasrauber/eagerpy/issues/6
    # assert np.array_equal(diag6.grad, [[0., 0.]])
