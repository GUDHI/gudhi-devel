""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Theo Lacombe, Marc Glisse

    Copyright (C) 2019 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

import numpy as np
from gudhi.quantization import kmeans_quantization, balanced_quantization
import pytest



# TO BE CONTINUED.

def test_balanced_quantization():
    dgm1 = np.array([[0,1], [1, 2], [2, 3], [2,3]])
    dgm1 = np.array([[0,1], [1, 2], [2, 3], [2,3]])
    dgm1 = np.array([[0,1], [1, 2], [2, 3], [2,3]])
    dgm1 = np.array([[0,1], [1, 2], [2, 3], [2,3]])
    dgm1 = np.array([[0,1], [1, 2], [2, 3], [2,3]])
    dgm1 = np.array([[0,1], [1, 2], [2, 3], [2,3]])

    assert balanced_quantization(dgm1, k=3, gamma=0.) == pytest.approx(0.)
    assert balanced_quantization(dgm1, k=3, gamma=0.1) == pytest.approx(0.)
    assert balanced_quantization(dgm1, k=3, gamma=1.) == pytest.approx(0.)
    assert balanced_quantization(dgm1, k=10, gamma=0.) == pytest.approx(0.)
    assert balanced_quantization(dgm1, k=1, gamma=0.) == pytest.approx(0.)
    assert balanced_quantization(dgm1, k=0, gamma=0.1) == pytest.approx(0.)
    assert balanced_quantization(dgm1, k=2, gamma=0.1) == pytest.approx(0.)


def test_kmeans_quantization():
    
    dgm1 = np.array([[0,1], [1, 2], [2, 3], [2,3]])

    assert kmeans_quantization(dgm1, k=3, gamma=0.) == pytest.approx(0.)

