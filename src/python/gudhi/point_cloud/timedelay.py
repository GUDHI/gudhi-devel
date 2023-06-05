# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Martin Royer, Yuichi Ike, Masatoshi Takenouchi
#
# Copyright (C) 2020 Inria, Copyright (C) 2020 Fujitsu Laboratories Ltd.
# Modification(s):
#   - YYYY/MM Author: Description of the modification

import numpy as np


class TimeDelayEmbedding:
    """Point cloud transformation class. Embeds time-series data in the R^d according to
    `Takens' Embedding Theorem <https://en.wikipedia.org/wiki/Takens%27s_theorem>`_ and obtains the
    coordinates of each point.

    Example
    -------

    Given delay=3 and skip=2, a point cloud which is obtained by embedding
    a scalar time-series into R^3 is as follows::

        time-series = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        point cloud = [[1, 4, 7],
                       [3, 6, 9]]

    Given delay=1 and skip=1, a point cloud which is obtained by embedding
    a 2D vector time-series data into R^4 is as follows::

        time-series = [[0, 1], [2, 3], [4, 5], [6, 7], [8, 9]]
        point cloud = [[0, 1, 2, 3],
                       [2, 3, 4, 5],
                       [4, 5, 6, 7],
                       [6, 7, 8, 9]]
    """

    def __init__(self, dim=3, delay=1, skip=1):
        """
        Constructor for the TimeDelayEmbedding class.

        Parameters:
            dim (int): `d` of R^d to be embedded. Optional (default=3).
            delay (int): Time-Delay embedding. Optional (default=1).
            skip (int): How often to skip embedded points. Optional (default=1).
        """
        self._dim = dim
        self._delay = delay
        self._skip = skip

    def __call__(self, ts):
        """Transform method for single time-series data.

        Parameters
        ----------
        ts : Iterable[float] or Iterable[Iterable[float]]
            A single time-series data, with scalar or vector values.

        Returns
        -------
        point cloud : n x dim numpy arrays
            Makes point cloud from a single time-series data.
        """
        return self._transform(np.array(ts))

    def fit(self, ts, y=None):
        return self
    
    def _transform(self, ts):
        """Guts of transform method."""
        if ts.ndim == 1:
            repeat = self._dim
        else:
            assert self._dim % ts.shape[1] == 0
            repeat = self._dim // ts.shape[1]
        end = len(ts) - self._delay * (repeat - 1)
        short = np.arange(0, end, self._skip)
        vertical = np.arange(0, repeat * self._delay, self._delay)
        return ts[np.add.outer(short, vertical)].reshape(len(short), -1)

    def transform(self, ts):
        """Transform method for multiple time-series data.

        Parameters
        ----------
        ts : Iterable[Iterable[float]] or Iterable[Iterable[Iterable[float]]]
            Multiple time-series data, with scalar or vector values.

        Returns
        -------
        point clouds : list of n x dim numpy arrays
            Makes point cloud from each time-series data.
        """
        return [self._transform(np.array(s)) for s in ts]
