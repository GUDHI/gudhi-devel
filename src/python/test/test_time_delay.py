from gudhi.point_cloud.timedelay import TimeDelayEmbedding
import numpy as np

def test_normal():
    # Sample array
    ts = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    # Normal case.
    prep = TimeDelayEmbedding()
    attractor = prep(ts)
    assert (attractor[0] == np.array([1, 2, 3])).all()
    assert (attractor[1] == np.array([2, 3, 4])).all()
    assert (attractor[2] == np.array([3, 4, 5])).all()
    assert (attractor[3] == np.array([4, 5, 6])).all()
    assert (attractor[4] == np.array([5, 6, 7])).all()
    assert (attractor[5] == np.array([6, 7, 8])).all()
    assert (attractor[6] == np.array([7, 8, 9])).all()
    assert (attractor[7] == np.array([8, 9, 10])).all()
    # Delay = 3
    prep = TimeDelayEmbedding(delay=3)
    attractor = prep(ts)
    assert (attractor[0] == np.array([1, 4, 7])).all()
    assert (attractor[1] == np.array([2, 5, 8])).all()
    assert (attractor[2] == np.array([3, 6, 9])).all()
    assert (attractor[3] == np.array([4, 7, 10])).all()
    # Skip = 3
    prep = TimeDelayEmbedding(skip=3)
    attractor = prep(ts)
    assert (attractor[0] == np.array([1, 2, 3])).all()
    assert (attractor[1] == np.array([4, 5, 6])).all()
    assert (attractor[2] == np.array([7, 8, 9])).all()
    # Delay = 2 / Skip = 2
    prep = TimeDelayEmbedding(delay=2, skip=2)
    attractor = prep(ts)
    assert (attractor[0] == np.array([1, 3, 5])).all()
    assert (attractor[1] == np.array([3, 5, 7])).all()
    assert (attractor[2] == np.array([5, 7, 9])).all()
