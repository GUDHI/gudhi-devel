from gudhi.wasserstein.quantization import _get_cells, _from_batch, _init_c, quantization
import numpy as np

__author__ = "Theo Lacombe"
__license__ = "MIT"


def test_get_cells():
    # Basic test
    X = np.array([[0, 0.01], [0, 1.], [0, 1.2], [0, 5.]])
    c = np.array([[0, 1.], [0, 5.]])
    res = _get_cells(X, c, withdiag=True, internal_p=2.)
    answer = [np.array([[0., 1.], [0, 1.2]]),
                                np.array([[0.,5.]]),
                                np.array([[0,0.01]])]
    assert all([np.array_equal(x,y) for x,y in zip(res,answer)])

    # Test withdiag=False
    X = np.array([[0, 0.01], [0, 1.], [0, 1.2], [0, 5.]])
    c = np.array([[0, 1.], [0, 5.]])
    res = _get_cells(X, c, withdiag=False, internal_p=2.)
    answer = [np.array([[0, 0.01], [0., 1.], [0, 1.2]]),
                                np.array([[0.,5.]])
              ]
    assert all([np.array_equal(x,y) for x,y in zip(res,answer)])

    # Test one cell has no points
    X = np.array([[0, 0.01], [0, 1.], [0, 1.2], [0, 5.]])
    c = np.array([[0, 1.], [0, 5.], [7, 9]])
    res = _get_cells(X, c, withdiag=True, internal_p=2.)
    print(res)
    answer = [np.array([[0., 1.], [0, 1.2]]),
              np.array([[0., 5.]]),
              np.empty((0,2)),
              np.array([[0, 0.01]])]
    assert all([np.array_equal(x, y) for x, y in zip(res, answer)])

    # Test input is empty
    X = np.array([])
    c = np.array([[0, 1.], [0, 5.]])
    res = _get_cells(X, c, withdiag=True, internal_p=2.)
    answer = [[]] * (c.shape[0] + 1)
    assert all([np.array_equal(x, y) for x, y in zip(res, answer)])

    # Test input wrong dimensions
    X = np.array([[0,1,2],[4,5,6]])
    c = np.array([[0, 1.], [0, 5.]])
    try:
        _get_cells(X, c, withdiag=True, internal_p=2.)
    except ValueError as e:
        assert str(e) == "Input batch must be of shape (n x 2), not (2,3)"


def test_from_batches():
    d1 = np.array([[0,1], [0,2], [0,3]])
    d2 = np.array([[0,1], [0,3], [0,6]])
    d3 = np.array([[0,2], [0,5], [0,3]])
    d4 = np.array([[0,1], [0,1], [0,1]])
    d5 = np.array([[1,2], [1,4], [1,8]])
    pdiagset = [d1, d2, d3, d4, d5]
    n = len(pdiagset)
    # test batches for odd number of diagrams
    batch_size = 2
    batches = np.array_split(np.arange(0, n, dtype=int), n // batch_size + (n % batch_size != 0))
    assert all([np.array_equal(x,y) for x,y in zip(batches, [np.array([0,1]), np.array([2,3]), np.array([4])])])
    # Test getting batch
    res = _from_batch(pdiagset, batches[0])
    answer = np.array([[0,1], [0,2], [0,3], [0,1], [0,3], [0,6]])
    assert np.array_equal(res, answer)
    res = _from_batch(pdiagset, batches[-1])
    answer = np.array([[1,2], [1,4], [1,8]])
    assert np.array_equal(res, answer)

    d6 = np.array([])
    pdiagset = [d1, d2, d3, d4, d5, d6]
    n = len(pdiagset)
    batch_size = 2
    batches = np.array_split(np.arange(0, n, dtype=int), n // batch_size + (n % batch_size != 0))
    assert all([np.array_equal(x, y) for x, y in zip(batches, [np.array([0, 1]),
                                                               np.array([2, 3]),
                                                               np.array([4, 5])])])
    # Test getting batch
    res = _from_batch(pdiagset, batches[1])
    answer = np.array([[0,2], [0,5], [0,3], [0,1], [0,1], [0,1]])
    assert np.array_equal(res, answer)
    # Test when one of the diagram is empty with wrong shape (get ignored in concat)
    res = _from_batch(pdiagset, batches[-1])
    answer = np.array([[1,2], [1,4], [1,8]])
    assert np.array_equal(res, answer)

    d7 = np.empty((0,2))
    d8 = np.array([[0,1],[0,42]])
    pdiagset = pdiagset + [d7, d8]
    n = len(pdiagset)
    batch_size = 3
    batches = np.array_split(np.arange(0, n, dtype=int), n // batch_size + (n % batch_size != 0))
    print(batches)
    assert all([np.array_equal(x, y) for x, y in zip(batches, [np.array([0, 1, 2]),
                                                               np.array([3, 4, 5]),
                                                               np.array([6, 7])])])

    res = _from_batch(pdiagset, batches[-1])
    answer = np.array([[0,1], [0,42]])
    assert np.array_equal(res, answer)

def test_init_c():
    d1 = np.array([[0, 0.001], [0, 1], [0, 2], [0, 3], [1, 1.01]])
    d2 = np.array([[0, 1], [0, 3], [0, 6]])
    d3 = np.array([[0, 2], [0, 5], [0, 3]])
    d4 = np.array([[0, 1], [0, 1], [0, 1]])
    d5 = np.array([[1, 2], [1, 4], [1, 8]])
    pdiagset = [d1, d2, d3, d4, d5]
    k = 2
    res = _init_c(pdiagset, k=k, internal_p=2)
    answer = np.array([[0, 2], [0, 3]])
    assert np.array_equal(res, answer)


def test_quantization():
    # test empty pdiagset
    pdiagset = []
    try:
        quantization(pdiagset)
    except ValueError as e:
        assert str(e) == "Input pdiagset is empty."

    # Test quantization on a vanilla example
    # Two clusters of points, around (0,1) and (0,3) ; and points close to the diagonal.
    d1 = np.array([[0, 1.001], [0, 3], [2, 2.001], [3, 3.001], [1, 1.01]])
    d2 = np.array([[0, 1], [0, 3.001], [0, 0.001]])
    d3 = np.array([[0, 0.999], [0, 3.002], [0, 2.998]])
    d4 = np.array([[0, 0.003], [0, 1.001], [0,3.004], [4, 4.01]])
    pdiagset = [d1, d2, d3, d4]
    c_final = quantization(pdiagset)
    answer = np.array([[0,1], [0,3]])
    my_tol = 2*1e-3  # result is something like 3.0015, so a bit greater than 1e-3.
    assert np.linalg.norm(c_final - answer) < my_tol

    # Same example but with empty diagrams in addition
    # Two clusters of points, around (0,1) and (0,3) ; and points close to the diagonal.
    d1 = np.array([[0, 1.001], [0, 3], [2, 2.001], [3, 3.001], [1, 1.01]])
    d2 = np.array([[0, 1], [0, 3.001], [0, 0.001]])
    d3 = np.array([[0, 0.999], [0, 3.002], [0, 2.998]])
    d4 = np.array([[0, 0.003], [0, 1.001], [0,3.004], [4, 4.01]])
    # Test with two different likely encoding of the empty diagram.
    d5 = np.empty((0,2))
    d6 = np.array([])
    pdiagset = [d1, d2, d3, d4, d5, d6]
    c_final = quantization(pdiagset)
    answer = np.array([[0,1], [0,3]])
    my_tol = 2*1e-3  # result is something like 3.0015, so a bit greater than 1e-3.
    assert np.linalg.norm(c_final - answer) < my_tol