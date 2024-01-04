# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Martin Royer

import copy

import numpy as np
import pandas as pd

from sklearn.base import BaseEstimator, TransformerMixin

from gudhi.representations.vector_methods import Atol, TopologicalVector


def _diag_to_dict_by_dim_format(pdiagram, max_dimension=None):
    """ transforms persistence diagrams in gudhi format (so list of Tuple(dimension, Tuple(x,y)))
     into dictionary of ndarray[[x_0, y_0], [x_1, y_1], ...] for each dimension, e.g.:
         [(1, (1.0577792405537423, 1.1003878733068035)),
          (0, (0.0, inf)),
          (0, (0.0, 1.0556057201636535)),
          (0, (0.0, 0.8756047102452433))]
     transforms into
         {'(D0)': array([[0.        ,        inf],
                         [0.        , 1.05560572],
                         [0.        , 0.87560471]]),
          '(D1)': array([[1.05777924, 1.10038787]])}
    """
    max_dimension = max(dim for (dim, _) in pdiagram) if max_dimension is None else max_dimension
    by_dim_pdiagram = {f"(D{i})": np.array([_ for (dim, _) in pdiagram if dim == i]) for i in range(0, max_dimension + 1)}
    return by_dim_pdiagram


class Archipelago(BaseEstimator, TransformerMixin):
    """
    Wrapper class for gudhi.representations.vector_methods in pandas format that is sklearn-API consistent.
    One provides persistence diagram vectorizers (by way of either `island` or `island_list`) and the target homology
    dimensions (`homology_dimensions`), and the Archipelago object will |fit on| and |transform = vectorize| dataframes
    of persistence diagrams.

    Parameters:
        homology_dimensions (int or list of int): The targeted persistence diagrams dimension(s).
            Short circuit the use of :class:`~gudhi.representations.preprocessing.DimensionSelector` when only one
            dimension matters (in other words, when `homology_dimensions` is an int).
        island: island for populating archipelago, i.e. object to vectorize the target in each homology
            dimensions. Must be `copy.deepcopy`-able. Will be ignored if island_list is given.
        island_list: island list for populating archipelago, i.e. list of object to vectorize the target
            in order of passed homology_dimensions.

    Examples
    >>> pdiagram1 = [(0, (0.0, 2.34)), (0, (0.0, 0.956)), (1, (0.536, 0.856)), (2, (1.202, 1.734))]
    >>> pdiagram2 = [(0, (0.0, 3.34)), (0, (0.0, 2.956)), (1, (0.536, 1.856)), (2, (1.202, 2.734))]
    >>> pdiagram3 = [(0, (1.0, 4.34)), (0, (2.0, 3.956)), (1, (1.536, 2.856)), (2, (3.202, 4.734))]
    >>> dict_pdiag1 = _diag_to_dict_by_dim_format(pdiagram1)
    >>> dict_pdiag2 = _diag_to_dict_by_dim_format(pdiagram2)
    >>> dict_pdiag3 = _diag_to_dict_by_dim_format(pdiagram3)
    >>> df_pdiags = pd.DataFrame([dict_pdiag1, dict_pdiag2, dict_pdiag3])
    >>> archipelago = Archipelago(homology_dimensions=range(3), island=Atol())
    >>> archipelago.fit(X=df_pdiags)
    >>> archipelago.transform(X=df_pdiags)
    >>> archipelago = Archipelago(homology_dimensions=range(2), island_list=[Atol(), TopologicalVector()])
    >>> archipelago.fit(X=df_pdiags)
    >>> archipelago.transform(X=df_pdiags)
    """

    def __init__(
            self,
            homology_dimensions=None,
            island=None,
            island_list=None
    ):
        if homology_dimensions is None:
            homology_dimensions = [0]
        self.homology_dimensions = homology_dimensions
        if island is None:
            island = Atol()
        self.island = island
        self.island_list = island_list

        if isinstance(self.homology_dimensions, int):
            self.dim_list_ = [ self.homology_dimensions ]
        else:
            self.dim_list_ = self.homology_dimensions

        if island_list is not None:
            if len(island_list) != len(self.dim_list_):
                raise ValueError("Archipelago initialization:"
                                 "`island_list` must be of the same length as number of targeted homology dimensions.")
            else:
                self.island_list_ = self.island_list
        else:
            self.island_list_ = [copy.deepcopy(self.island) for _ in range(len(self.dim_list_))]
        self.archipelago_ = {}

    def fit(self, X, y=None):
        """
        Calibration step: create and fit `island` vectorizer to each matching diagram element

        Args:
            X (pandas.DataFrame of diagrams): input sets of diagrams to be fitted on.
            y: possibly labels for each diagram

        Returns:
            self
        """

        for key, dgms in X.items():
            if not sum(len(dgm) for dgm in dgms):
                print(f"[Archipelago] No point to fit for {key}, dropping this key.")
                continue

            for (homology_dimension, island) in zip(self.dim_list_, self.island_list_):
                if f"(D{homology_dimension})" not in key:
                    continue
                print(f"[Archipelago] Fitting key {key} matching homology dimension {homology_dimension}.")
                try:
                    island.fit(X=dgms, y=y)
                except ValueError as ve:
                    print(f'[Archipelago] Fit of {key} returned "{ve}", archipelago will ignore this key.')
                    continue
                self.archipelago_[key] = island
        return self

    def transform(self, X, y=None):
        """
        Apply measure vectorisation on a dictionary of list of measures.

        Args:
            X (pandas.DataFrame of diagrams): input sets of diagrams to vectorize.
            y: Ignored, present for API consistency by convention.

        Returns:
            dict of numpy array with shape (number of measures) x (n_clusters_by_atol).
        """
        archipelago_vectorized = pd.DataFrame(index=X.index)
        for key, dgms in X.items():
            if key not in self.archipelago_.keys():
                continue
            vectorized_dgms = self.archipelago_[key].transform(dgms)
            col_keys = [f"{key} Center {i + 1}" for i in range(vectorized_dgms.shape[1])]
            archipelago_vectorized.loc[:, col_keys] = vectorized_dgms
        return archipelago_vectorized

    def get_feature_names_out(self):
        # hack to get access to skl `set_output`
        pass
