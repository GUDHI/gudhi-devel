""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2016 Inria

    Modification(s):
      - 2024/07 Vincent Rouvreau: Separate Delaunay and Alpha tests in a different python file
      - 2025/04 Hannah Schreiber: Add tests to verify possibility of tensor input
      - YYYY/MM Author: Description of the modification
"""

__author__ = "Vincent Rouvreau"
__maintainer__ = "Vincent Rouvreau, Hannah Schreiber"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "GPL v3"


import numpy as np
import pytest

from gudhi import AlphaComplex


def _alpha_get_point(precision):
    point_list = [[0, 0], [1, 0], [0, 1], [1, 1]]
    alpha_complex = AlphaComplex(points=point_list, precision=precision)

    assert point_list[0] == alpha_complex.get_point(0)
    assert point_list[1] == alpha_complex.get_point(1)
    assert point_list[2] == alpha_complex.get_point(2)
    assert point_list[3] == alpha_complex.get_point(3)

    with pytest.raises(IndexError):
        alpha_complex.get_point(len(point_list))


def test_alpha_get_point():
    for precision in ["fast", "safe", "exact"]:
        _alpha_get_point(precision)


def _3d_tetrahedrons(precision):
    points = 10 * np.random.rand(10, 3)
    alpha = AlphaComplex(points=points, precision=precision)
    st_alpha = alpha.create_simplex_tree(default_filtration_value=False)
    # New AlphaComplex for get_point to work
    delaunay = AlphaComplex(points=points, precision=precision)
    st_delaunay = delaunay.create_simplex_tree(default_filtration_value=True)

    delaunay_tetra = []
    for sk in st_delaunay.get_skeleton(4):
        if len(sk[0]) == 4:
            tetra = [
                delaunay.get_point(sk[0][0]),
                delaunay.get_point(sk[0][1]),
                delaunay.get_point(sk[0][2]),
                delaunay.get_point(sk[0][3]),
            ]
            delaunay_tetra.append(sorted(tetra, key=lambda tup: tup[0]))

    alpha_tetra = []
    for sk in st_alpha.get_skeleton(4):
        if len(sk[0]) == 4:
            tetra = [
                alpha.get_point(sk[0][0]),
                alpha.get_point(sk[0][1]),
                alpha.get_point(sk[0][2]),
                alpha.get_point(sk[0][3]),
            ]
            alpha_tetra.append(sorted(tetra, key=lambda tup: tup[0]))

    # Check the tetrahedrons from one list are in the second one
    assert len(alpha_tetra) == len(delaunay_tetra)
    for tetra_from_del in delaunay_tetra:
        assert tetra_from_del in alpha_tetra


def test_3d_tetrahedrons():
    for precision in ["fast", "safe", "exact"]:
        _3d_tetrahedrons(precision)


def test_inconsistency_points_and_weights():
    points = [
        [1.0, 1.0, 0.0],
        [7.0, 0.0, 0.0],
        [4.0, 6.0, 0.0],
        [9.0, 6.0, 0.0],
        [0.0, 14.0, 0.0],
        [2.0, 19.0, 0.0],
        [9.0, 17.0, 0.0],
    ]
    with pytest.raises(ValueError):
        # 7 points, 8 weights, on purpose
        alpha = AlphaComplex(points=points, weights=[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])

    with pytest.raises(ValueError):
        # 7 points, 6 weights, on purpose
        alpha = AlphaComplex(points=points, weights=[1.0, 2.0, 3.0, 4.0, 5.0, 6.0])


def _weighted_doc_example(precision):
    stree = AlphaComplex(
        points=[
            [1.0, -1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
            [1.0, 1.0, 1.0],
            [2.0, 2.0, 2.0],
        ],
        weights=[4.0, 4.0, 4.0, 4.0, 1.0],
        precision=precision,
    ).create_simplex_tree()

    assert stree.filtration([0, 1, 2, 3]) == pytest.approx(-1.0)
    assert stree.filtration([0, 1, 3, 4]) == pytest.approx(95.0)
    assert stree.filtration([0, 2, 3, 4]) == pytest.approx(95.0)
    assert stree.filtration([1, 2, 3, 4]) == pytest.approx(95.0)


def test_weighted_doc_example():
    for precision in ["fast", "safe", "exact"]:
        _weighted_doc_example(precision)


def test_float_relative_precision():
    assert AlphaComplex.get_float_relative_precision() == 1e-5
    # Must be > 0.
    with pytest.raises(ValueError):
        AlphaComplex.set_float_relative_precision(0.0)
    # Must be < 1.
    with pytest.raises(ValueError):
        AlphaComplex.set_float_relative_precision(1.0)

    points = [[1, 1], [7, 0], [4, 6], [9, 6], [0, 14], [2, 19], [9, 17]]
    st = AlphaComplex(points=points).create_simplex_tree()
    filtrations = list(st.get_filtration())

    # Get a better precision
    AlphaComplex.set_float_relative_precision(1e-15)
    assert AlphaComplex.get_float_relative_precision() == 1e-15

    st = AlphaComplex(points=points).create_simplex_tree()
    filtrations_better_resolution = list(st.get_filtration())

    assert len(filtrations) == len(filtrations_better_resolution)
    for idx in range(len(filtrations)):
        # check simplex is the same
        assert filtrations[idx][0] == filtrations_better_resolution[idx][0]
        # check filtration is about the same with a relative precision of the worst case
        assert filtrations[idx][1] == pytest.approx(
            filtrations_better_resolution[idx][1], rel=1e-5
        )


def test_numpy_arrays():
    points = np.array(
        [
            [1.0, -1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
            [1.0, 1.0, 1.0],
            [2.0, 2.0, 2.0],
        ]
    )
    weights = np.array([4.0, 4.0, 4.0, 4.0, 1.0])
    alpha = AlphaComplex(points=points, weights=weights)


def test_tensors():
    try:
        import torch

        points = (torch.rand((5, 2)) * 2 - 1).requires_grad_()
        weights = (torch.rand((5)) * 2 - 1).requires_grad_()
        alpha = AlphaComplex(points=points, weights=weights)
    except ImportError:
        pass

    try:
        import tensorflow as tf

        points = tf.random.uniform(shape=[5, 2])
        weights = tf.random.uniform(shape=[5])
        alpha = AlphaComplex(points=points, weights=weights)
    except ImportError:
        pass
