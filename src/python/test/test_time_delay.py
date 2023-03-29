from gudhi.point_cloud.timedelay import TimeDelayEmbedding
import numpy as np


def test_normal():
    # Sample array
    ts = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    # Normal case.
    prep = TimeDelayEmbedding()
    pointclouds = prep(ts)
    assert (pointclouds[0] == np.array([1, 2, 3])).all()
    assert (pointclouds[1] == np.array([2, 3, 4])).all()
    assert (pointclouds[2] == np.array([3, 4, 5])).all()
    assert (pointclouds[3] == np.array([4, 5, 6])).all()
    assert (pointclouds[4] == np.array([5, 6, 7])).all()
    assert (pointclouds[5] == np.array([6, 7, 8])).all()
    assert (pointclouds[6] == np.array([7, 8, 9])).all()
    assert (pointclouds[7] == np.array([8, 9, 10])).all()
    # Delay = 3
    prep = TimeDelayEmbedding(delay=3)
    pointclouds = prep(ts)
    assert (pointclouds[0] == np.array([1, 4, 7])).all()
    assert (pointclouds[1] == np.array([2, 5, 8])).all()
    assert (pointclouds[2] == np.array([3, 6, 9])).all()
    assert (pointclouds[3] == np.array([4, 7, 10])).all()
    # Skip = 3
    prep = TimeDelayEmbedding(skip=3)
    pointclouds = prep(ts)
    assert (pointclouds[0] == np.array([1, 2, 3])).all()
    assert (pointclouds[1] == np.array([4, 5, 6])).all()
    assert (pointclouds[2] == np.array([7, 8, 9])).all()
    # Delay = 2 / Skip = 2
    prep = TimeDelayEmbedding(delay=2, skip=2)
    pointclouds = prep(ts)
    assert (pointclouds[0] == np.array([1, 3, 5])).all()
    assert (pointclouds[1] == np.array([3, 5, 7])).all()
    assert (pointclouds[2] == np.array([5, 7, 9])).all()

    # Vector series
    ts = np.arange(0, 10).reshape(-1, 2)
    prep = TimeDelayEmbedding(dim=4)
    prep.fit([ts])
    assert (prep.transform([ts])[0] == [[0, 1, 2, 3], [2, 3, 4, 5], [4, 5, 6, 7], [6, 7, 8, 9]]).all()
