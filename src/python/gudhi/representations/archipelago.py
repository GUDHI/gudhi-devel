# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Martin Royer

import copy

import numpy as np

from sklearn.base import BaseEstimator, TransformerMixin

from gudhi.representations.vector_methods import Atol, TopologicalVector, BettiCurve


class Archipelago(BaseEstimator, TransformerMixin):
    """
    Transformer that dictionary-wraps persistence diagram vectorizers, i.e. objects from gudhi.representations.vector_methods.
    One provides persistence diagram vectorizers (by way of either `island` or `island_dict`), and the Archipelago object will
    |fit on| and |transform = vectorize| lists or series of persistence diagrams.
    The object is sklearn-API consistent.

    Parameters:
        island: island for populating archipelago, i.e. object to vectorize the target in each homology
            dimensions. Must be `copy.deepcopy`-able. *Will be ignored if island_list is given*.
        island_dict: island dict for populating archipelago, i.e. dictionary of objects to vectorize persistence
            diagrams according to homology_dimensions.

    Examples
    >>> pdiagram1 = [(0, (0.0, 2.34)), (0, (0.0, 0.956)), (1, (0.536, 0.856)), (2, (1.202, 1.734))]
    >>> pdiagram2 = [(0, (0.0, 3.34)), (0, (0.0, 2.956)), (1, (0.536, 1.856)), (2, (1.202, 2.734))]
    >>> pdiagram3 = [(0, (1.0, 4.34)), (0, (2.0, 3.956)), (1, (1.536, 2.856)), (2, (3.202, 4.734))]
    >>> list_pdiags = [pdiagram1, pdiagram2, pdiagram3]
    >>> archipelago = Archipelago(island=Atol())
    >>> archipelago.fit(X=list_pdiags)
    >>> archipelago.transform(X=list_pdiags)
    >>> archipelago = Archipelago(island_dict={2: BettiCurve(resolution=4), 0:Atol()})
    >>> import pandas as pd
    >>> series_pdiags = pd.Series(list_pdiags)
    >>> archipelago.set_output(transform="pandas")
    >>> archipelago.fit(X=series_pdiags)
    >>> archipelago.transform(X=series_pdiags)
    """

    def __init__(
            self,
            island=None,
            island_dict=None
    ):
        if island is None and island_dict is None:
            island = Atol()
        self.island = island
        self.island_dict = island_dict
        self.archipelago_ = {}
        self._running_transform_names = ""

    def fit(self, X, y=None):
        """
        Calibration step: create and fit `island` vectorizer to each matching diagram element

        Args:
            X (list of diagrams): input persistence diagrams to fit vectorizers on.
            y: possibly labels for each diagram

        Returns:
            self
        """

        max_dimension = max(dim for pdiagram in X for (dim, _) in pdiagram)
        by_dim_list_pdiags = [[
            np.array([_ for (dim, _) in pdiagram if dim == dimension]) for dimension in range(0, max_dimension + 1)
        ] for pdiagram in X]
        for dimension in range(0, max_dimension + 1):
            this_dim_list_pdiags = [pdiags[dimension] for pdiags in by_dim_list_pdiags]
            if not len(this_dim_list_pdiags):
                continue
            if self.island_dict is not None and dimension in self.island_dict.keys():
                island = self.island_dict[dimension]
            elif self.island_dict is not None:
                continue
            else:
                island = copy.deepcopy(self.island)
            try:
                island.fit(X=this_dim_list_pdiags, y=y)
            except ValueError as ve:
                print(f"[Archipelago] Fit of homology dimension {dimension} returned {ve}. Will ignore this key.")
                continue
            print(f"[Archipelago] Fit of homology dimension {dimension} with object {island.__class__} succeeded.")
            self.archipelago_[dimension] = island
        return self

    def transform(self, X, y=None):
        """
        Apply measure vectorisation on a dictionary of list of measures.

        Args:
            X (list of diagrams): input persistence diagrams to vectorize.
            y: Ignored, present for API consistency by convention.

        Returns:
            vectors : array of shape (len(X), n_features) where the columns features are vectorized homology dimension
                in increasing order.
        """

        max_dimension = max(dim for pdiagram in X for (dim, _) in pdiagram)
        by_dim_list_pdiags = [[
            np.array([_ for (dim, _) in pdiagram if dim == dimension]) for dimension in range(0, max_dimension + 1)
        ] for pdiagram in X]

        archipelago_vectorized = []
        running_transform_names = []
        for dimension in range(0, max_dimension + 1):
            if dimension not in self.archipelago_.keys():
                print(f"[Archipelago] Encounters homology dimension {dimension} that has not been fitted on. Will ignore this key")
                continue
            this_dim_list_pdiags = [pdiags[dimension] for pdiags in by_dim_list_pdiags]
            vectorized_dgms = self.archipelago_[dimension].transform(this_dim_list_pdiags)
            running_transform_names += [f"{dimension} Center {i + 1}" for i in range(vectorized_dgms.shape[1])]
            archipelago_vectorized.append(vectorized_dgms)
        self._running_transform_names = running_transform_names
        return np.concatenate(archipelago_vectorized, axis=1)

    def get_feature_names_out(self):
        return self._running_transform_names
