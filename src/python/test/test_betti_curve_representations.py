import numpy as np
import scipy.interpolate

from gudhi.representations.vector_methods import BettiCurve2

def test_betti_curve_is_irregular_betti_curve_followed_by_interpolation():
    m = 10
    n = 1000
    pinf = 0.05
    pzero = 0.05
    res = 100
    success = True
    
    pds = []
    for i in range(0, m):
        pd = np.zeros((n, 2))
        pd[:, 0] = np.random.uniform(0, 10, n)
        pd[:, 1] = np.random.uniform(pd[:, 0], 10, n)
        pd[np.random.uniform(0, 1, n) < pzero, 0] = 0
        pd[np.random.uniform(0, 1, n) < pinf, 1] = np.inf
        pds.append(pd)
    
    bc = BettiCurve2(None)
    bc.fit(pds)
    bettis = bc.transform(pds)

    bc2 = BettiCurve2(None)
    bettis2 = bc2.fit_transform(pds)
    success = success and (bc2.grid_ == bc.grid_).all()
    success = success and (bettis2 == bettis).all()

    for i in range(0, m):
        grid = np.linspace(pds[i].min(), pds[i].max() + 1, res)
        bc_gridded = BettiCurve2(grid)
        bettis_gridded = bc_gridded(pds[i])

        interp = scipy.interpolate.interp1d(bc.grid_, bettis[i, :], kind="previous", fill_value="extrapolate")
        bettis_interp = np.array(interp(grid), dtype=int)
        success = success and (bettis_interp == bettis_gridded).all()

    assert(success)
