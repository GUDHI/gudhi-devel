# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Martin Royer, Yuichi Ike, Masatoshi Takenouchi
#
# Copyright (C) 2020 Inria, Copyright (C) 2020 Fujitsu Laboratories Ltd.
# Modification(s):
#   - YYYY/MM Author: Description of the modification

import numpy as np

class TimeDelayEmbedding:
    """Point cloud transformation class.
    Embeds time-series data in the R^d according to Takens' Embedding Theorem
    and obtains the coordinates of each point.
    Parameters
    ----------
    dim : int, optional (default=3)
        `d` of R^d to be embedded.
    delay : int, optional (default=1)
        Time-Delay embedding.
    skip : int, optional (default=1)
        How often to skip embedded points.
        Given delay=3 and skip=2, an point cloud which is obtained by embedding 
        a single time-series data into R^3 is as follows.
        
    .. code-block:: none
    
    time-series = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    point clouds = [[1, 4, 7], 
                    [3, 6, 9]]
                   
    """
    def __init__(self, dim=3, delay=1, skip=1):
        self._dim = dim
        self._delay = delay
        self._skip = skip

    def __call__(self, ts):
        """Transform method for single time-series data.
        Parameters
        ----------
        ts : list[float]
            A single time-series data.
        Returns
        -------
        point clouds : list of n x 2 numpy arrays
            Makes point cloud every a single time-series data.
        Raises
        -------
        TypeError
            If the parameter's type does not match the desired type.
        """
        ndts = np.array(ts)
        if ndts.ndim == 1:
            return self._transform(ndts)
        else:
            raise TypeError("Expects 1-dimensional array.")

    def fit(self, ts, y=None):
        return self
    
    def _transform(self, ts):
        """Guts of transform method."""
        return ts[
            np.add.outer(
                np.arange(0, len(ts)-self._delay*(self._dim-1), self._skip),
                np.arange(0, self._dim*self._delay, self._delay))
        ]

    def transform(self, ts):
        """Transform method for multiple time-series data.
        Parameters
        ----------
        ts : list[list[float]]
            Multiple time-series data.
        Attributes
        ----------
        ndts : 
             The ndts means that all time series need to have exactly 
             the same size. 
        Returns
        -------
        point clouds : list of n x 3 numpy arrays
            Makes point cloud every a single time-series data.
        Raises
        -------
        TypeError
            If the parameter's type does not match the desired type.
        """
        ndts = np.array(ts)
        if ndts.ndim == 2:
            return np.apply_along_axis(self._transform, 1, ndts)
        else:
            raise TypeError("Expects 2-dimensional array.")
