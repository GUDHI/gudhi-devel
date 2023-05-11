# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2022 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

from .. import RipsComplex
from sklearn.base import BaseEstimator, TransformerMixin

import numpy as np
# joblib is required by scikit-learn
from joblib import Parallel, delayed

# Mermaid sequence diagram - https://mermaid-js.github.io/mermaid-live-editor/
# sequenceDiagram
#     USER->>RipsPersistence: fit_transform(X)
#     Note right of RipsPersistence: homology_dimensions=[i,j]
#     RipsPersistence->>thread1: _tranform(X[0])
#     RipsPersistence->>thread2: _tranform(X[1])
#     Note right of RipsPersistence: ...
#     thread1->>RipsPersistence: [array( Hi(X[0]) ), array( Hj(X[0]) )]
#     thread2->>RipsPersistence: [array( Hi(X[1]) ), array( Hj(X[1]) )]
#     Note right of RipsPersistence: ...
#     RipsPersistence->>USER: [[array( Hi(X[0]) ), array( Hj(X[0]) )],<br/> [array( Hi(X[1]) ), array( Hj(X[1]) )],<br/> ...]


class RipsPersistence(BaseEstimator, TransformerMixin):
    """
    This is a class for computing the persistence diagrams from a Vietoris-Rips complex.
    """

    def __init__(
        self,
        homology_dimensions,
        max_edge_length=float('inf'),
        nb_collapse=-1,
        homology_coeff_field=11,
        min_persistence=0.0,
        n_jobs=None,
    ):
        """
        Constructor for the RipsPersistence class.

        Parameters:
            homology_dimensions (int or list of int): The returned persistence diagrams dimension(s).
                Short circuit the use of :class:`~gudhi.representations.preprocessing.DimensionSelector` when only one
                dimension matters (in other words, when `homology_dimensions` is an int).
            max_edge_length (float): Rips value. Default is +Inf.
            nb_collapse (int): The number of :func:`~gudhi.SimplexTree.collapse_edges` iterations to perform on the
                SimplexTree. Default is -1, which means "automatic" (a relatively good enough number of iterations is
                choosen).
            homology_coeff_field (int): The homology coefficient field. Must be a prime number. Default value is 11.
            min_persistence (float): The minimum persistence value to take into account (strictly greater than
                `min_persistence`). Default value is `0.0`. Set `min_persistence` to `-1.0` to see all values.
            n_jobs (int): cf. https://joblib.readthedocs.io/en/latest/generated/joblib.Parallel.html
        """
        self.homology_dimensions = homology_dimensions
        self.max_edge_length = max_edge_length
        self.nb_collapse = nb_collapse
        self.homology_coeff_field = homology_coeff_field
        self.min_persistence = min_persistence
        self.n_jobs = n_jobs

    def fit(self, X, Y=None):
        """
        Nothing to be done, but useful when included in a scikit-learn Pipeline.
        """
        return self

    def __get_stree_from_points(self, points, max_dimension):
        # nb_collapse "automatic" case management
        if self.nb_collapse < 0:
            nb_collapse = 1
        else:
            nb_collapse = self.nb_collapse
        
        rips = RipsComplex(points=points, max_edge_length = self.max_edge_length)
        if max_dimension > 1:
            stree = rips.create_simplex_tree(max_dimension=1)
            stree.collapse_edges(nb_iterations = nb_collapse)
            stree.expansion(max_dimension)
        else:
            stree = rips.create_simplex_tree(max_dimension=max_dimension)
        
        stree.compute_persistence(
            homology_coeff_field=self.homology_coeff_field, min_persistence=self.min_persistence
        )
        return stree

    def __transform(self, points):
        stree = (self.__get_stree_from_points)(points, (max(self.homology_dimensions) + 1))
        return [
            stree.persistence_intervals_in_dimension(dim) for dim in self.homology_dimensions
        ]

    def __transform_only_this_dim(self, points):
        stree = (self.__get_stree_from_points)(points, (self.homology_dimensions + 1))
        return stree.persistence_intervals_in_dimension(self.homology_dimensions)

    def transform(self, X, Y=None):
        """Compute all the Vietoris-Rips complexes and their associated persistence diagrams.

        :param X: List of point clouds.
        :type X: list of list of float OR list of numpy.ndarray

        :return: Persistence diagrams in the format:

              - If `homology_dimensions` was set to `n`: `[array( Hn(X[0]) ), array( Hn(X[1]) ), ...]` 
              - If `homology_dimensions` was set to `[i, j]`: `[[array( Hi(X[0]) ), array( Hj(X[0]) )], [array( Hi(X[1]) ), array( Hj(X[1]) )], ...]`
        :rtype: list of (,2) array_like or list of list of (,2) array_like
        """
        # Depends on homology_dimensions is an integer or a list of integer (else case)
        if isinstance(self.homology_dimensions, int):
            # threads is preferred as Rips construction and persistence computation releases the GIL
            return Parallel(n_jobs=self.n_jobs, prefer="threads")(
                delayed(self.__transform_only_this_dim)(points) for points in X
            )
        else:
            # threads is preferred as Rips construction and persistence computation releases the GIL
            return Parallel(n_jobs=self.n_jobs, prefer="threads")(delayed(self.__transform)(points) for points in X)

