from gudhi.point_cloud.timedelay import TimeDelayEmbedding

def test_normal():
    # Sample array
    ts = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    # Normal case.
    prep = TimeDelayEmbedding()
    attractor = prep(ts)
    assert (attractor[0] == [1, 2, 3])
    assert (attractor[1] == [2, 3, 4])
    assert (attractor[2] == [3, 4, 5])
    assert (attractor[3] == [4, 5, 6])
    assert (attractor[4] == [5, 6, 7])
    assert (attractor[5] == [6, 7, 8])
    assert (attractor[6] == [7, 8, 9])
    assert (attractor[7] == [8, 9, 10])
    # Delay = 3
    prep = TimeDelayEmbedding(delay=3)
    attractor = prep(ts)
    assert (attractor[0] == [1, 4, 7])
    assert (attractor[1] == [2, 5, 8])
    assert (attractor[2] == [3, 6, 9])
    assert (attractor[3] == [4, 7, 10])
    # Skip = 3
    prep = TimeDelayEmbedding(skip=3)
    attractor = prep(ts)
    assert (attractor[0] == [1, 2, 3])
    assert (attractor[1] == [4, 5, 6])
    assert (attractor[2] == [7, 8, 9])
    # Delay = 2 / Skip = 2
    prep = TimeDelayEmbedding(delay=2, skip=2)
    attractor = prep(ts)
    assert (attractor[0] == [1, 3, 5])
    assert (attractor[1] == [3, 5, 7])
    assert (attractor[2] == [5, 7, 9])
