"""This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
Author(s):       Henrique Ennes & Raphaël Tinarrage

Copyright (C) 2016 Inria & Institute of Science and Technology Austria
"""

import numpy as np
import pytest
from gudhi.liedetect import Orthonormalize

"""
Fixes data set.
"""


@pytest.fixture
def fix_fit_pts():
    """
    Builds data set to fit orthonormalizer.
    """
    thetas = np.linspace(0, 2 * np.pi, 100)
    return np.array([(np.cos(theta), 7 * np.sin(theta)) for theta in thetas])


@pytest.fixture
def fix_transform_pts(fix_fit_pts):
    """
    Builds data set to apply orthonormalizer.
    Simply the twice expansion of the original set.
    """
    return 2 * fix_fit_pts


@pytest.fixture
def fix_orhtonormalize(fix_fit_pts):
    return Orthonormalize(fix_fit_pts)


"""
Test functions
"""


def test_normalizer(fix_orhtonormalize):
    pts = fix_orhtonormalize.fit_transform()
    assert abs(np.mean(np.linalg.norm(pts, axis=1)) - 1) < 1e-4


def test_transform(fix_orhtonormalize, fix_transform_pts):
    ortho = fix_orhtonormalize
    ortho.fit()
    transformed_pts = ortho.transform(fix_transform_pts)

    assert abs(np.mean(np.linalg.norm(transformed_pts, axis=1) - 2)) < 1e-4

    inverse = ortho.inverse_transform(transformed_pts)

    assert abs(np.mean(np.linalg.norm(inverse - fix_transform_pts, axis=1))) < 1e-4
