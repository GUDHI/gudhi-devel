# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2021 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

from sklearn.base import BaseEstimator, TransformerMixin
from .. import RipsComplex
from joblib import Parallel, delayed


class RipsPersistence(BaseEstimator, TransformerMixin):
    """
    This is a class for computing the persistence diagrams from a Rips complex.
    """

    def __init__(
        self,
        max_rips_dimension=2,
        max_edge_length=float("inf"),
        collapse_edges=True,
        max_persistence_dimension=0,
        only_this_dim=-1,
        homology_coeff_field=11,
        min_persistence=0.0,
        n_jobs=None,
    ):
        """
        Constructor for the RipsPersistence class.

        Parameters:
            points (Iterable of iterable of float): A point cloud.
            max_rips_dimension (int): Rips expansion until this dimension.
            max_edge_length (float): Maximal edge length to be taken into account.
            max_persistence_dimension (int): The returned persistence diagrams maximal dimension. Default value is `0`.
                Ignored if `only_this_dim` is set.
            only_this_dim (int): The returned persistence diagrams dimension. If `only_this_dim` is set,
                `max_persistence_dimension` will be ignored. 
                Short circuit the use of :class:`~gudhi.sklearn.post_processing.DimensionSelector` when only one
                dimension matters.
            homology_coeff_field (int): The homology coefficient field. Must be a prime number. Default value is 11.
            min_persistence (float): The minimum persistence value to take into account (strictly greater than
                `min_persistence`). Default value is `0.0`. Sets `min_persistence` to `-1.0` to see all values.
            n_jobs (int): cf. https://joblib.readthedocs.io/en/latest/generated/joblib.Parallel.html
        """
        self.max_rips_dimension=max_rips_dimension
        self.max_edge_length = max_edge_length
        self.collapse_edges = collapse_edges
        self.max_persistence_dimension = max_persistence_dimension
        self.only_this_dim = only_this_dim
        self.homology_coeff_field = homology_coeff_field
        self.min_persistence = min_persistence
        self.n_jobs = n_jobs

    def fit(self, X, Y=None):
        """
        Nothing to be done, but useful when included in a scikit-learn Pipeline.
        """
        return self

    def __transform(self, points):
        rips = RipsComplex(points=points, max_edge_length=self.max_edge_length)
        if self.collapse_edges:
            stree = rips.create_simplex_tree(max_dimension=1)
            stree.collapse_edges()
            stree.expansion(self.max_rips_dimension)
        else:
            stree = rips.create_simplex_tree(max_dimension=self.max_rips_dimension)
        stree.compute_persistence(
            homology_coeff_field=self.homology_coeff_field, min_persistence=self.min_persistence
        )
        return [
            stree.persistence_intervals_in_dimension(dim) for dim in range(self.max_persistence_dimension + 1)
        ]

    def __transform_only_this_dim(self, points):
        rips = RipsComplex(points=points, max_edge_length=self.max_edge_length)
        if self.collapse_edges:
            stree = rips.create_simplex_tree(max_dimension=1)
            stree.collapse_edges()
            stree.expansion(self.max_rips_dimension)
        else:
            stree = rips.create_simplex_tree(max_dimension=self.max_rips_dimension)
        stree.compute_persistence(
            homology_coeff_field=self.homology_coeff_field, min_persistence=self.min_persistence
        )
        return stree.persistence_intervals_in_dimension(self.only_this_dim)

    def transform(self, X, Y=None):
        """
        Compute all the Rips complexes and their associated persistence diagrams.

        Parameters:
            X (list of list of double OR list of numpy.ndarray): List of point clouds.

        Returns:
            Persistence diagrams in the format:
            - If `only_this_dim` was set to `n`: `[array( Hn(X[0]) ), array( Hn(X[1]) ), ...]` 
            - else: `[[array( H0(X[0]) ), array( H1(X[0]) ), ...], [array( H0(X[1]) ), array( H1(X[1]) ), ...], ...]` 
        """

        if self.only_this_dim == -1:
            # threads is preferred as rips construction and persistence computation releases the GIL
            return Parallel(n_jobs=self.n_jobs, prefer="threads")(delayed(self.__transform)(cells) for cells in X)
        else:
            # threads is preferred as rips construction and persistence computation releases the GIL
            return Parallel(n_jobs=self.n_jobs, prefer="threads")(
                delayed(self.__transform_only_this_dim)(cells) for cells in X
            )
