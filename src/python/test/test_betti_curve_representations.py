import numpy as np
import scipy.interpolate
import pytest

from gudhi.representations.vector_methods import BettiCurve

def test_betti_curve_is_irregular_betti_curve_followed_by_interpolation():
    m = 10
    n = 1000
    pinf = 0.05
    pzero = 0.05
    res = 100
    
    pds = []
    for i in range(0, m):
        pd = np.zeros((n, 2))
        pd[:, 0] = np.random.uniform(0, 10, n)
        pd[:, 1] = np.random.uniform(pd[:, 0], 10, n)
        pd[np.random.uniform(0, 1, n) < pzero, 0] = 0
        pd[np.random.uniform(0, 1, n) < pinf, 1] = np.inf
        pds.append(pd)
    
    bc = BettiCurve(resolution=None, predefined_grid=None)
    bc.fit(pds)
    bettis = bc.transform(pds)

    bc2 = BettiCurve(resolution=None, predefined_grid=None)
    bettis2 = bc2.fit_transform(pds)
    assert((bc2.grid_ == bc.grid_).all())
    assert((bettis2 == bettis).all())

    for i in range(0, m):
        grid = np.linspace(pds[i][np.isfinite(pds[i])].min(), pds[i][np.isfinite(pds[i])].max() + 1, res)
        bc_gridded = BettiCurve(predefined_grid=grid)
        bc_gridded.fit([])
        bettis_gridded = bc_gridded(pds[i])

        interp = scipy.interpolate.interp1d(bc.grid_, bettis[i, :], kind="previous", fill_value="extrapolate")
        bettis_interp = np.array(interp(grid), dtype=int)
        assert((bettis_interp == bettis_gridded).all())


def test_empty_with_predefined_grid():
    random_grid = np.sort(np.random.uniform(0, 1, 100))
    bc = BettiCurve(predefined_grid=random_grid)
    bettis = bc.fit_transform([])
    assert((bc.grid_ == random_grid).all())
    assert((bettis == 0).all())

    
def test_empty():
    bc = BettiCurve(resolution=None, predefined_grid=None)
    bettis = bc.fit_transform([])
    assert(bc.grid_ == [-np.inf])
    assert((bettis == 0).all())

def test_wrong_value_of_predefined_grid():
    with pytest.raises(ValueError):
        BettiCurve(predefined_grid=[1, 2, 3])
