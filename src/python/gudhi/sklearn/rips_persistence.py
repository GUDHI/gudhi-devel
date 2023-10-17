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
    This is a class for constructing Vietoris-Rips complexes and computing the persistence diagrams from them.
    """

    def __init__(
        self,
        homology_dimensions,
        threshold=float('inf'),
        input_type='point cloud',
        num_collapses=1,
        homology_coeff_field=11,
        expand_extra_dimension = True,
        n_jobs=None,
    ):
        """
        Constructor for the RipsPersistence class.

        Parameters:
            homology_dimensions (int or list of int): The returned persistence diagrams dimension(s).
                Short circuit the use of :class:`~gudhi.representations.preprocessing.DimensionSelector` when only one
                dimension matters (in other words, when `homology_dimensions` is an int).
            threshold (float): Rips maximal edge length value. Default is +Inf.
            input_type (str): Can be 'point cloud' when inputs are point clouds, or 'lower distance matrix', when
                inputs are lower triangular distance matrix (can be full square, but the upper part of the distance
                matrix will not be considered). Default is 'point cloud'.
            num_collapses (int): Specify the number of :func:`~gudhi.SimplexTree.collapse_edges` iterations to perform
                on the SimplexTree. Default value is 1 (a relatively good enough number of iterations).
            homology_coeff_field (int): The homology coefficient field. Must be a prime number. Default value is 11.
            expand_extra_dimension (bool): :func:`~gudhi.SimplexTree.expansion` is performed at
                `max(homology_dimensions) + 1` if true, and `max(homology_dimensions)` otherwise. Default is true.
            n_jobs (int): Number of jobs to run in parallel. `None` (default value) means `n_jobs = 1` unless in a
                joblib.parallel_backend context. `-1` means using all processors. cf.
                https://joblib.readthedocs.io/en/latest/generated/joblib.Parallel.html for more details.
        """
        self.homology_dimensions = homology_dimensions
        self.threshold = threshold
        self.input_type = input_type
        self.num_collapses = num_collapses
        self.homology_coeff_field = homology_coeff_field
        self.expand_extra_dimension = expand_extra_dimension
        self.n_jobs = n_jobs

    def fit(self, X, Y=None):
        """
        Nothing to be done, but useful when included in a scikit-learn Pipeline.
        """
        return self

    def __transform(self, inputs):
        # persistence_dim_max is False by default, but True when expand_extra_dimension is False
        persistence_dim_max = self.expand_extra_dimension
        max_dimension = max(self.dim_list_)
        if self.expand_extra_dimension:
            max_dimension += 1

        if self.input_type == 'point cloud':
            rips = RipsComplex(points=inputs, max_edge_length = self.threshold)
        elif self.input_type == 'lower distance matrix':
            rips = RipsComplex(distance_matrix=inputs, max_edge_length = self.threshold)
        else:
            raise ValueError("Only 'point cloud' and  'lower distance matrix' are valid input_type")
        
        if max_dimension > 1:
            stree = rips.create_simplex_tree(max_dimension=1)
            stree.collapse_edges(nb_iterations = self.num_collapses)
            stree.expansion(max_dimension)
        else:
            stree = rips.create_simplex_tree(max_dimension=max_dimension)

        # Specific case where, despite expansion(max_dimension), stree has a lower dimension
        if max_dimension > stree.dimension():
            persistence_dim_max = True

        stree.compute_persistence(
            homology_coeff_field=self.homology_coeff_field,
            persistence_dim_max=persistence_dim_max
        )

        return [
            stree.persistence_intervals_in_dimension(dim) for dim in self.dim_list_
        ]

    def transform(self, X, Y=None):
        """Compute all the Vietoris-Rips complexes and their associated persistence diagrams.

        :param X: list of point clouds as Euclidean coordinates or distance matrices.
        :type X: list of list of float OR list of numpy.ndarray

        :return: Persistence diagrams in the format:

              - If `homology_dimensions` was set to `n`: `[array( Hn(X[0]) ), array( Hn(X[1]) ), ...]` 
              - If `homology_dimensions` was set to `[i, j]`: `[[array( Hi(X[0]) ), array( Hj(X[0]) )], [array( Hi(X[1]) ), array( Hj(X[1]) )], ...]`
        :rtype: list of numpy ndarray of shape (,2) or list of list of numpy ndarray of shape (,2)
        """
        # Depends on homology_dimensions is an integer or a list of integer (else case)
        if isinstance(self.homology_dimensions, int):
            unwrap = True
            self.dim_list_ = [ self.homology_dimensions ]
        else:
            unwrap = False
            self.dim_list_ = self.homology_dimensions

        # threads is preferred as Rips construction and persistence computation releases the GIL
        res = Parallel(n_jobs=self.n_jobs, prefer="threads")(delayed(self.__transform)(inputs) for inputs in X)

        if unwrap:
            res = [d[0] for d in res]
        return res
