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

    """
    def __init__(self, dim=3, delay=1, skip=1):
        self._dim = dim
        self._delay = delay
        self._skip = skip

    def __call__(self, *args, **kwargs):
        return self.transform(*args, **kwargs)

    def _transform(self, ts):
        """Guts of transform method."""
        return ts[
            np.add.outer(
                np.arange(0, len(ts)-self._delay*(self._dim-1), self._skip),
                np.arange(0, self._dim*self._delay, self._delay))
        ]

    def transform(self, ts):
        """Transform method.

        Parameters
        ----------
        ts : list[float] or list[list[float]]
            A single or multiple time-series data.

        Returns
        -------
        point clouds : list[list[float, float, float]] or list[list[list[float, float, float]]]
            Makes point cloud every a single time-series data.
        """
        ndts = np.array(ts)
        if ndts.ndim == 1:
            # for single.
            return self._transform(ndts).tolist()
        else:
            # for multiple.
            return np.apply_along_axis(self._transform, 1, ndts).tolist()
