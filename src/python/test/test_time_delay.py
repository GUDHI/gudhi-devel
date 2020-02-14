from gudhi.point_cloud.timedelay import TimeDelayEmbedding
import numpy as np

def test_normal():
    # Sample array
    ts = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    # Normal case.
    prep = TimeDelayEmbedding()
    attractor = prep(ts)
    assert (attractor[0] == np.array([1, 2, 3])
            print(attractor[0].all()))
    assert (attractor[1] == np.array([2, 3, 4]))
    assert (attractor[2] == np.array([3, 4, 5]))
    assert (attractor[3] == np.array([4, 5, 6]))
    assert (attractor[4] == np.array([5, 6, 7]))
    assert (attractor[5] == np.array([6, 7, 8]))
    assert (attractor[6] == np.array([7, 8, 9]))
    assert (attractor[7] == np.array([8, 9, 10]))
    # Delay = 3
    prep = TimeDelayEmbedding(delay=3)
    attractor = prep(ts)
    assert (attractor[0] == np.array([1, 4, 7]))
    assert (attractor[1] == np.array([2, 5, 8]))
    assert (attractor[2] == np.array([3, 6, 9]))
    assert (attractor[3] == np.array([4, 7, 10]))
    # Skip = 3
    prep = TimeDelayEmbedding(skip=3)
    attractor = prep(ts)
    assert (attractor[0] == np.array([1, 2, 3]))
    assert (attractor[1] == np.array([4, 5, 6]))
    assert (attractor[2] == np.array([7, 8, 9]))
    # Delay = 2 / Skip = 2
    prep = TimeDelayEmbedding(delay=2, skip=2)
    attractor = prep(ts)
    assert (attractor[0] == np.array([1, 3, 5]))
    assert (attractor[1] == np.array([3, 5, 7]))
    assert (attractor[2] == np.array([5, 7, 9]))
