# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2022 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

from .._ripser import _lower, _full, _sparse, _lower_to_coo, _lower_cone_radius
from ..flag_filtration.edge_collapse import reduce_graph
from .. import SimplexTree
import math
import numpy as np
from typing import Union, Iterable, Literal, Optional
from sklearn.base import BaseEstimator, TransformerMixin
from scipy.sparse import coo_matrix
from scipy.spatial import cKDTree
from scipy.spatial.distance import pdist, squareform

# joblib is required by scikit-learn
from joblib import Parallel, delayed

# Mermaid sequence diagram - https://mermaid-js.github.io/mermaid-live-editor/
# sequenceDiagram
#   participant USER
#   participant R as RipsPersistence
#   USER->>R: fit_transform(X)
#   Note right of R: homology_dimensions=[i,j]
#   R->>thread1: _tranform(X[0])
#   R->>thread2: _tranform(X[1])
#   Note right of R: ...
#   thread1->>R: [array( Hi(X[0]) ), array( Hj(X[0]) )]
#   thread2->>R: [array( Hi(X[1]) ), array( Hj(X[1]) )]
#   Note right of R: ...
#   R->>USER: [[array( Hi(X[0]) ), array( Hj(X[0]) )],<br/> [array( Hi(X[1]) ), array( Hj(X[1]) )],<br/>...]


class RipsPersistence(BaseEstimator, TransformerMixin):
    """
    This is a class for constructing Vietoris-Rips complexes and computing the persistence diagrams from them.
    """

    def __init__(
        self,
        homology_dimensions: Union[int, Iterable[int]],
        threshold: float = float('inf'),
        input_type: Literal[
            "point cloud", "full distance matrix", "lower distance matrix", "distance coo_matrix"
        ] = 'point cloud',
        num_collapses: Union[int, Literal["auto"]] = 'auto',
        homology_coeff_field: int = 11,
        n_jobs: Optional[int] = None,
    ):
        """Constructor for the RipsPersistence class.

        Parameters:
            homology_dimensions: The returned persistence diagrams dimension(s).
                Short circuit the use of :class:`~gudhi.representations.preprocessing.DimensionSelector` when only one
                dimension matters (in other words, when `homology_dimensions` is an int).
            threshold: Rips maximal edge length value. Default is +Inf. Ignored if input_type is
                'distance coo_matrix'.
            input_type: Can be 'point cloud' when inputs are point clouds, 'full distance matrix',
                'lower distance matrix' when inputs are lower triangular distance matrix (can be full square,
                but the upper part will not be considered), or 'distance coo_matrix' for a distance matrix in SciPy's
                sparse format, which should contain each edge at most once (avoid the symmetric) and no diagonal entry.
                Default is 'point cloud'.
            num_collapses: Specify the number of iterations of
                :func:`~gudhi.flag_filtration.edge_collapse.reduce_graph` (edge collapse) to perform on the graph.
                Default value is 'auto'.
            homology_coeff_field: The homology coefficient field. Must be a prime number. Default value is 11.
            n_jobs: Number of jobs to run in parallel. `None` (default value) means `n_jobs = 1` unless in a
                joblib.parallel_backend context. `-1` means using all processors. cf.
                https://joblib.readthedocs.io/en/latest/generated/joblib.Parallel.html for more details.
        """
        self.homology_dimensions = homology_dimensions
        self.threshold = threshold
        self.input_type = input_type
        self.num_collapses = num_collapses
        self.homology_coeff_field = homology_coeff_field
        self.n_jobs = n_jobs

    def fit(self, X, Y=None):
        """
        Nothing to be done, but useful when included in a scikit-learn Pipeline.
        """
        return self

    def __transform(self, inp):
        # TODO: give the user more control over the strategy
        # Should we use threshold in the sparse case?
        num_collapses = self.num_collapses
        input_type = self.input_type
        threshold = self.threshold
        n = inp.shape[0] if input_type == 'distance coo_matrix' else len(inp)
        max_dimension = min(max(self.dim_list_), max(0, n-3))
        # Ripser needs to encode simplices and coefficients in 128 bits, which may not always fit
        # Instead of a 256 bit version which may not always suffice either, fall back to SimplexTree
        use_simplex_tree = math.comb(n, min(n // 2, max_dimension + 2)) >= (1 << (128 - (self.homology_coeff_field - 2).bit_length()))
        if num_collapses == 'auto':
            num_collapses = 1 if max_dimension > (not use_simplex_tree) else 0
            # or num_collapses=max_dimension-1 maybe?
        elif max_dimension == 0:
            num_collapses = 0

        # Points -> distance matrix
        if input_type == 'point cloud':
            if threshold < float('inf'):
                # Hope that the user gave a useful threshold
                tree = cKDTree(inp)

                ## V1: Returns self-loops and every edge twice (symmetry)
                # inp = tree.sparse_distance_matrix(tree, max_distance=threshold, output_type="coo_matrix")
                # mask = inp.row < inp.col
                # inp = coo_matrix((inp.data[mask], (inp.row[mask], inp.col[mask])), shape=inp.shape)

                # V2: Gets the right edges, but forgets the distances
                pairs = tree.query_pairs(r=threshold, output_type='ndarray')
                data = np.ravel(np.linalg.norm(np.diff(inp[pairs], axis=1), axis=-1))
                inp = coo_matrix((data, (pairs[:,0], pairs[:,1])), shape=(n,)*2)

                input_type = 'distance coo_matrix'
            else:
                inp = squareform(pdist(inp))
                input_type = 'full distance matrix'

        # Dense -> sparse
        if input_type in ('full distance matrix', 'lower distance matrix'):
            # After this filtration value, all complexes are cones, nothing happens
            if input_type == 'full distance matrix':
                inp = np.asarray(inp)
                cone_radius = inp.max(-1).min()
            else:
                cone_radius = _lower_cone_radius(inp)
            sparsify = use_simplex_tree or num_collapses > 0 or threshold < cone_radius # heuristic
            threshold = min(threshold, cone_radius)
            if sparsify:
                # For 'full' we could use i, j = np.triu_indices_from(inp, k=1), etc
                i, j, f = _lower_to_coo(inp, threshold)
                inp = coo_matrix((f, (i, j)), shape=(n,) * 2)
                input_type = 'distance coo_matrix'

        if num_collapses > 0:
            assert input_type == 'distance coo_matrix'
            inp = reduce_graph(inp, num_collapses)

        if use_simplex_tree:
            st = SimplexTree()
            # Use create_from_array in case of full matrix?
            # (not important since this fallback mostly matters in high dimension, where we use edge-collapse anyway)
            st.insert_batch(np.arange(n).reshape(1,-1), np.zeros(n))
            st.insert_edges_from_coo_matrix(inp)
            st.expansion(max_dimension + 1)
            st.compute_persistence(homology_coeff_field=self.homology_coeff_field, persistence_dim_max=max_dimension>=st.dimension())
            return [ st.persistence_intervals_in_dimension(dim) for dim in self.dim_list_ ]

        if input_type == 'full distance matrix':
            ## Possibly transpose for performance?
            # if inp.strides[0] > inp.strides[1]: # or the reverse?
            #     inp = inp.T
            dgm = _full(inp, max_dimension=max_dimension, max_edge_length=threshold, homology_coeff_field=self.homology_coeff_field)
        elif input_type == 'lower distance matrix':
            dgm = _lower(inp, max_dimension=max_dimension, max_edge_length=threshold, homology_coeff_field=self.homology_coeff_field)
        elif input_type == 'distance coo_matrix':
            # Switch to coo_array (danger: row/col seem deprecated)?
            dgm = _sparse(inp.row, inp.col, inp.data, inp.shape[0], max_dimension=max_dimension, max_edge_length=threshold, homology_coeff_field=self.homology_coeff_field)
        else:
            raise ValueError("Only 'point cloud', 'lower distance matrix', 'full distance matrix' and 'distance coo_matrix' are valid input_type") # move to __init__?

        # dgm stops at n-2
        return [dgm[dim] if dim < len(dgm) else np.empty((0,2)) for dim in self.dim_list_]

    def transform(self, X, Y=None):
        """Compute all the Vietoris-Rips complexes and their associated persistence diagrams.

        :param X: list of point clouds as Euclidean coordinates or distance matrices.
        :type X: list of list of float OR list of numpy.ndarray

        :return: Persistence diagrams in the format:

              - If `homology_dimensions` was set to `n`: `[array( Hn(X[0]) ), array( Hn(X[1]) ), ...]`
              - If `homology_dimensions` was set to `[i, j]`:
                `[[array( Hi(X[0]) ), array( Hj(X[0]) )], [array( Hi(X[1]) ), array( Hj(X[1]) )], ...]`
        :rtype: list of numpy ndarray of shape (,2) or list of list of numpy ndarray of shape (,2)
        """
        # Depends if homology_dimensions is an integer or a list of integers (else case)
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

    def get_feature_names_out(self):
        """Provide column names for implementing sklearn's set_output API."""
        return [f"H{i}" for i in self.dim_list_]
