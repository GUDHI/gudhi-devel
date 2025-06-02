# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2024 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

from typing import Union, Iterable, Literal, Optional, Any
import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from joblib import Parallel, delayed

from .. import DelaunayCechComplex, AlphaComplex

# Mermaid sequence diagram - https://mermaid-js.github.io/mermaid-live-editor/
# sequenceDiagram
#   participant USER
#   participant R as CechPersistence
#   USER->>R: fit_transform(X)
#   Note right of R: homology_dimensions=[i,j]
#   R->>thread1: _tranform(X[0])
#   R->>thread2: _tranform(X[1])
#   Note right of R: ...
#   thread1->>R: [array( Hi(X[0]) ), array( Hj(X[0]) )]
#   thread2->>R: [array( Hi(X[1]) ), array( Hj(X[1]) )]
#   Note right of R: ...
#   R->>USER: [[array( Hi(X[0]) ), array( Hj(X[0]) )],<br/> [array( Hi(X[1]) ), array( Hj(X[1]) )],<br/>...]


class CechPersistence(BaseEstimator, TransformerMixin):
    """
    This is a class for computing the same persistent homology as the Čech complex, while being significantly smaller,
    using internally a :class:`~gudhi.DelaunayCechComplex`.
    """

    def __init__(
        self,
        homology_dimensions: Union[int, Iterable[int]],
        precision: Literal["fast", "safe", "exact"] = "safe",
        output_squared_values: bool = False,
        threshold: float = float("inf"),
        homology_coeff_field: int = 11,
        n_jobs: Optional[int] = None,
    ):
        """Constructor for the CechPersistence class.

        Parameters:
            homology_dimensions: The returned persistence diagrams dimension(s).
                Short circuit the use of :class:`~gudhi.representations.preprocessing.DimensionSelector` when only
                one dimension matters (in other words, when `homology_dimensions` is an int).
            precision: Complex precision can be 'fast', 'safe' or 'exact'. Default is 'safe'.
            output_squared_values: Square filtration values when `True`. Default is `False` (contrary to its
                default value in :class:`~gudhi.DelaunayCechComplex`).
            threshold: The maximum filtration value the simplices shall not exceed. Default is set to infinity, and
                there is very little point using anything else since it does not save time.
                Notice that the filtration values (equal to radii, by default, or squared radii) will be different in
                function of `output_squared_values` value (when `False`, by default, or `True`).
            homology_coeff_field: The homology coefficient field. Must be a prime number. Default value is 11.
            n_jobs: Number of jobs to run in parallel. `None` (default value) means `n_jobs = 1` unless in a
                joblib.parallel_backend context. `-1` means using all processors. cf.
                https://joblib.readthedocs.io/en/latest/generated/joblib.Parallel.html for more details.
        """
        self.homology_dimensions = homology_dimensions
        self.precision = precision
        self.output_squared_values = output_squared_values
        self.threshold = threshold
        self.homology_coeff_field = homology_coeff_field
        self.n_jobs = n_jobs

        if self.precision not in ["safe", "fast", "exact"]:
            raise ValueError("Unknown precision")

        # Done twice (in __init__ and fit), but exception is better the sooner
        dim_list = np.asarray(self.homology_dimensions, dtype=int)
        if dim_list.ndim not in [0, 1]:
            raise ValueError(f"Invalid dimension. Got {self.homology_dimensions=}, expected type=int|Iterable[int].")

    def fit(self, X, Y=None):
        """
        Fit the `CechPersistence` class in function of `homology_dimensions` type.
        """
        # Must be in the `fit` part, as `transform` should be const and as `__init__` is not called on a parallel grid
        # search for instance
        self._dim_list = np.asarray(self.homology_dimensions, dtype=int)
        self._unwrap = False
        if self._dim_list.ndim == 0:
            self._unwrap = True
            self._dim_list = self._dim_list.reshape(1)
        return self

    def __transform(self, inputs):
        max_dimension = np.max(self._dim_list) + 1

        # alpha complex points take floats
        pts = np.asarray(inputs, dtype=np.float64)
        delaunay_cech = DelaunayCechComplex(points=pts, precision=self.precision)
        stree = delaunay_cech.create_simplex_tree(
            max_alpha_square=self.threshold, output_squared_values=self.output_squared_values
        )

        persistence_dim_max = False
        # Specific case where, despite expansion(max_dimension), stree has a lower dimension
        if max_dimension > stree.dimension():
            persistence_dim_max = True

        stree.compute_persistence(
            homology_coeff_field=self.homology_coeff_field, persistence_dim_max=persistence_dim_max
        )

        return [stree.persistence_intervals_in_dimension(dim) for dim in self._dim_list]

    def transform(self, X, Y=None):
        """Compute all the Čech complexes and their associated persistence diagrams.

        :param X: list of point clouds as Euclidean coordinates.
        :type X: list of list of float OR list of numpy.ndarray

        :return: Persistence diagrams in the format:

              - If `homology_dimensions` was set to `n`: `[array( Hn(X[0]) ), array( Hn(X[1]) ), ...]`
              - If `homology_dimensions` was set to `[i, j]`:
                `[[array( Hi(X[0]) ), array( Hj(X[0]) )], [array( Hi(X[1]) ), array( Hj(X[1]) )], ...]`
        :rtype: list of numpy ndarray of shape (,2) or list of list of numpy ndarray of shape (,2)
        """
        # threads is preferred as Delaunay-Cech construction and persistence computation releases the GIL
        res = Parallel(n_jobs=self.n_jobs, prefer="threads")(delayed(self.__transform)(inputs) for inputs in X)

        # cf. `fit`
        if self._unwrap:
            res = [d[0] for d in res]
        return res

    def get_feature_names_out(self):
        """Provide column names for implementing sklearn's set_output API."""
        return [f"H{i}" for i in self._dim_list]


class WeightedCechPersistence(BaseEstimator, TransformerMixin):
    """
    This is a class for computing the same persistent homology as the Weighted Čech complex, while being significantly
    smaller, using internally a Weighted version of :class:`~gudhi.AlphaComplex`.
    """

    def __init__(
        self,
        homology_dimensions: Union[int, Iterable[int]],
        precision: Literal["fast", "safe", "exact"] = "safe",
        output_squared_values: bool = True,
        threshold: float = float("inf"),
        homology_coeff_field: int = 11,
        n_jobs: Optional[int] = None,
    ):
        """Constructor for the WeightedCechPersistence class.

        Parameters:
            homology_dimensions: The returned persistence diagrams dimension(s).
                Short circuit the use of :class:`~gudhi.representations.preprocessing.DimensionSelector` when only
                one dimension matters (in other words, when `homology_dimensions` is an int).
            precision: Complex precision can be 'fast', 'safe' or 'exact'. Default is 'safe'.
            output_squared_values: Square filtration values when `True`. Default is `True`  (contrary to the unweighted
                version :class:`~gudhi.sklearn.cech_persistence.CechPersistence`).
            threshold: The maximum filtration value the simplices shall not exceed. Default is set to infinity, and
                there is very little point using anything else since it does not save time.
                Notice that the filtration values (equal to squared radii, by default, or radii) will be different in
                function of `output_squared_values` value (when `True`, by default, or `False`).
            homology_coeff_field: The homology coefficient field. Must be a prime number. Default value is 11.
            n_jobs: Number of jobs to run in parallel. `None` (default value) means `n_jobs = 1` unless in a
                joblib.parallel_backend context. `-1` means using all processors. cf.
                https://joblib.readthedocs.io/en/latest/generated/joblib.Parallel.html for more details.
        """
        self.homology_dimensions = homology_dimensions
        self.precision = precision
        self.output_squared_values = output_squared_values
        self.threshold = threshold
        self.homology_coeff_field = homology_coeff_field
        self.n_jobs = n_jobs

        if self.precision not in ["safe", "fast", "exact"]:
            raise ValueError("Unknown precision")

        # Done twice (in __init__ and fit), but exception is better the sooner
        dim_list = np.asarray(self.homology_dimensions, dtype=int)
        if dim_list.ndim not in [0, 1]:
            raise ValueError(f"Invalid dimension. Got {self.homology_dimensions=}, expected type=int|Iterable[int].")


    def fit(self, X, Y=None):
        """
        Fit the `WeightedCechPersistence` class in function of `homology_dimensions` type.
        """
        # Must be in the `fit` part, as `transform` should be const and as `__init__` is not called on a parallel grid
        # search for instance
        self._dim_list = np.asarray(self.homology_dimensions, dtype=int)
        self._unwrap = False
        if self._dim_list.ndim == 0:
            self._unwrap = True
            self._dim_list = self._dim_list.reshape(1)
        return self

    def __transform(self, inputs):
        max_dimension = np.max(self._dim_list) + 1

        # alpha complex weights and points take floats
        wgts = np.asarray(inputs[:, -1], dtype=np.float64)
        pts = np.asarray(inputs[:, :-1], dtype=np.float64)
        weighted_alpha = AlphaComplex(points=pts, weights=wgts, precision=self.precision)

        stree = weighted_alpha.create_simplex_tree(
            max_alpha_square=self.threshold, output_squared_values=self.output_squared_values
        )

        persistence_dim_max = False
        # Specific case where, despite expansion(max_dimension), stree has a lower dimension
        if max_dimension > stree.dimension():
            persistence_dim_max = True

        stree.compute_persistence(
            homology_coeff_field=self.homology_coeff_field, persistence_dim_max=persistence_dim_max
        )

        return [stree.persistence_intervals_in_dimension(dim) for dim in self._dim_list]

    def transform(self, X, Y=None):
        """Compute all the Weighted Čech complexes and their associated persistence diagrams.

        :param X: list of point clouds as Euclidean coordinates, plus the weight (as the last column value).
        :type X: list of list of float OR list of numpy.ndarray

        :return: Persistence diagrams in the format:

              - If `homology_dimensions` was set to `n`: `[array( Hn(X[0]) ), array( Hn(X[1]) ), ...]`
              - If `homology_dimensions` was set to `[i, j]`:
                `[[array( Hi(X[0]) ), array( Hj(X[0]) )], [array( Hi(X[1]) ), array( Hj(X[1]) )], ...]`
        :rtype: list of numpy ndarray of shape (,2) or list of list of numpy ndarray of shape (,2)
        """
        # threads is preferred as Alpha complex construction and persistence computation releases the GIL
        res = Parallel(n_jobs=self.n_jobs, prefer="threads")(delayed(self.__transform)(inputs) for inputs in X)

        # cf. `fit`
        if self._unwrap:
            res = [d[0] for d in res]
        return res

    def get_feature_names_out(self):
        """Provide column names for implementing sklearn's set_output API."""
        return [f"H{i}" for i in self._dim_list]
