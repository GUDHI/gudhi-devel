""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2021 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

import gudhi as gd
import numpy as np
import matplotlib as mplt
import pytest
import warnings


def test_array_handler():
    diags = np.array([[1, 2], [3, 4], [5, 6]], float)
    arr_diags, input_type = gd.persistence_graphical_tools._format_handler(diags)
    assert input_type == 1
    for idx in range(len(diags)):
        assert arr_diags[idx][0] == 0
        np.testing.assert_array_equal(arr_diags[idx][1], diags[idx])

    diags = [(1.0, 2.0), (3.0, 4.0), (5.0, 6.0)]
    arr_diags, input_type = gd.persistence_graphical_tools._format_handler(diags)
    assert input_type == 1
    for idx in range(len(diags)):
        assert arr_diags[idx][0] == 0
        assert arr_diags[idx][1] == diags[idx]

    diags = [(0, (1.0, 2.0)), (0, (3.0, 4.0)), (0, (5.0, 6.0))]
    arr_diags, input_type = gd.persistence_graphical_tools._format_handler(diags)
    assert input_type == 0
    assert arr_diags == diags

    diags = [[(1.0, 2.0), (3.0, 4.0), (5.0, 6.0)], [(1.0, 2.0), (3.0, 4.0), (5.0, 6.0)]]
    arr_diags, input_type = gd.persistence_graphical_tools._format_handler(diags)
    assert input_type == 2
    assert len(arr_diags) == 6
    for idx in range(3):
        assert arr_diags[idx][0] == 0
        assert arr_diags[idx + 3][0] == 1
        np.testing.assert_array_equal(arr_diags[idx][1], arr_diags[idx + 3][1])


def test_min_birth_max_death():
    diags = [
        (0, (0.0, float("inf"))),
        (0, (0.0983494, float("inf"))),
        (0, (0.0, 0.122545)),
        (0, (0.0, 0.12047)),
        (0, (0.0, 0.118398)),
        (0, (0.118398, 1.0)),
        (0, (0.0, 0.117908)),
        (0, (0.0, 0.112307)),
        (0, (0.0, 0.107535)),
        (0, (0.0, 0.106382)),
    ]
    assert gd.persistence_graphical_tools._min_birth_max_death(diags) == (0.0, 1.0)
    assert gd.persistence_graphical_tools._min_birth_max_death(diags, band=4.0) == (0.0, 5.0)


def test_limit_min_birth_max_death():
    diags = [
        (0, (2.0, float("inf"))),
        (0, (2.0, float("inf"))),
    ]
    assert gd.persistence_graphical_tools._min_birth_max_death(diags) == (2.0, 3.0)
    assert gd.persistence_graphical_tools._min_birth_max_death(diags, band=4.0) == (2.0, 6.0)


def test_limit_to_max_intervals():
    diags = [
        (0, (0.0, float("inf"))),
        (0, (0.0983494, float("inf"))),
        (0, (0.0, 0.122545)),
        (0, (0.0, 0.12047)),
        (0, (0.0, 0.118398)),
        (0, (0.118398, 1.0)),
        (0, (0.0, 0.117908)),
        (0, (0.0, 0.112307)),
        (0, (0.0, 0.107535)),
        (0, (0.0, 0.106382)),
    ]
    # check no warnings if max_intervals equals to the diagrams number
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        truncated_diags = gd.persistence_graphical_tools._limit_to_max_intervals(
            diags, 10, key=lambda life_time: life_time[1][1] - life_time[1][0]
        )
        # check diagrams are not sorted
        assert truncated_diags == diags

    # check warning if max_intervals lower than the diagrams number
    with pytest.warns(UserWarning) as record:
        truncated_diags = gd.persistence_graphical_tools._limit_to_max_intervals(
            diags, 5, key=lambda life_time: life_time[1][1] - life_time[1][0]
        )
        # check diagrams are truncated and sorted by life time
        assert truncated_diags == [
            (0, (0.0, float("inf"))),
            (0, (0.0983494, float("inf"))),
            (0, (0.118398, 1.0)),
            (0, (0.0, 0.122545)),
            (0, (0.0, 0.12047)),
        ]
    assert len(record) == 1


def _limit_plot_persistence(function):
    pplot = function(persistence=[])
    assert isinstance(pplot, mplt.axes.SubplotBase)
    pplot = function(persistence=[], legend=True)
    assert isinstance(pplot, mplt.axes.SubplotBase)
    pplot = function(persistence=[(0, float("inf"))])
    assert isinstance(pplot, mplt.axes.SubplotBase)
    pplot = function(persistence=[(0, float("inf"))], legend=True)
    assert isinstance(pplot, mplt.axes.SubplotBase)


def test_limit_plot_persistence():
    for function in [gd.plot_persistence_barcode, gd.plot_persistence_diagram, gd.plot_persistence_density]:
        _limit_plot_persistence(function)


def _non_existing_persistence_file(function):
    with pytest.raises(FileNotFoundError):
        function(persistence_file="a_filename_that_should_not_exist.weird_extension")


def test_non_existing_persistence_file():
    for function in [gd.plot_persistence_barcode, gd.plot_persistence_diagram, gd.plot_persistence_density]:
        _non_existing_persistence_file(function)


def _sklearn_one_homology_dim_plot_persistence(function):
    # from gudhi.sklearn.rips_persistence import RipsPersistence
    # X = [[1., 1.], [7., 0.], [4., 6.], [9., 6.], [0., 14.], [2., 19.], [9., 17.]]
    # diag = RipsPersistence(homology_dimensions=1).fit_transform([X])
    # plot_persistence_diagram(diag[0]) # should work
    diags = [np.array([[11.0, 12.0], [6.0, 7.0]]), np.array([[10.0, 11.0], [5.0, 6.0]])]
    for diag in diags:
        ax = function(diag)


def test_sklearn_one_homology_dim_plot_persistence():
    for function in [gd.plot_persistence_barcode, gd.plot_persistence_diagram]:
        _sklearn_one_homology_dim_plot_persistence(function)


def _sklearn_several_homology_dim_plot_persistence(function):
    # from gudhi.sklearn.rips_persistence import RipsPersistence
    # X = [[1., 1.], [7., 0.], [4., 6.], [9., 6.], [0., 14.], [2., 19.], [9., 17.]]
    # diag = RipsPersistence(homology_dimensions=[1,0]).fit_transform([X])
    # plot_persistence_diagram(diag[0]) # should work
    diags = [
        [np.array([[11.0, 12.0], [6.0, 7.0]]), np.array([[0.0, 5.0], [0.0, 6.0], [0.0, float("inf")]])],
        [np.array([[11.0, 12.0], [6.0, 7.0]]), np.array([[0.0, 5.0], [0.0, 6.0], [0.0, float("inf")]])],
    ]
    for diag in diags:
        ax = function(diag)
        assert ax.get_legend().get_title().get_text() == "Range"


def test_sklearn_several_homology_dim_plot_persistence():
    for function in [gd.plot_persistence_barcode, gd.plot_persistence_diagram]:
        _sklearn_several_homology_dim_plot_persistence(function)
