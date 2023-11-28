# coding: utf-8
"""
@author: Martin Royer
"""

import copy

import numpy as np
import pandas as pd

from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.cluster import KMeans

from .vector_methods import Atol

# @martin@todo: doc, exemple d'utilisation


class Archipelago(BaseEstimator, TransformerMixin):
    """
    Wrapper class for Atol dictionaries that is sklearn-API consistent.

    Example
    --------
    >>> archipelago = Archipelago(settler_class=Atol, settler_params={"quantiser": KMeans(n_init='auto'), "weighting_method": "iidproba"}),
    >>> archipelago.fit(X=train_diags)
    >>> archipelago.transform(X=diag_data)
    """

    def __init__(self, settler_class: type = Atol, settler_params=None):
        """
        Constructor for the Archipelago class that sets the future number of cluster per Atol created.

        @param settler_class: elementary class for populating the archipelago.
        @param settler_params: parameters dict with which to instantiate the populating class.
        """
        if settler_params is None:
            settler_params = {"quantiser": KMeans(n_init='auto'), "weighting_method": "iidproba"}
        self.settler_class = settler_class
        self.settler_params = settler_params
        self.archipelago_: dict[str: Atol] = {}

    def fit(self, X, y=None):
        """
        Calibration step: create and fit Atol elements to the corresponding diagram element

        Parameters:
            X (dictionary of diagrams): input sets of diagrams to be fit by the Atol algorithm.
            y: possibly labels for each diagram

        Returns:
            self
        """
        for key, dgms in X.items():
            if not (dgms.apply(type) == np.ndarray).all():
                print(f"[Archipelago] {key} does not yield only ndarrays.")
                continue
            if not np.sum(len(dgm) for dgm in dgms):
                print(f"[Archipelago] No point to fit for {key}, dropping this key.")
                continue

            atol = self.settler_class(**copy.deepcopy(self.settler_params))

            n_points = np.concatenate(dgms).shape[0]
            if atol.quantiser.n_clusters > n_points:
                print(f"[Archipelago] {key} has only {n_points} points for `fit`, reducing n_clusters.")
                atol.quantiser.n_clusters = np.min([atol.quantiser.n_clusters, n_points])
            try:
                atol.fit(X=dgms, y=y)
            except ValueError as ve:
                print(f'[Archipelago] Fit of {key} returned "{ve}", archipelago will ignore this key.')
                continue
            # atol.inertias *= 10
            self.archipelago_[key] = atol
        return self

    def transform(self, X, y=None):
        """
        Apply measure vectorisation on a dictionary of list of measures.

        Parameters:
            X (dict of list N x d numpy arrays): dict of input measures in R^d to vectorize according to Atol centers
                (measures can have different N).

        Returns:
            dict of numpy array with shape (number of measures) x (n_clusters_by_atol).
        """
        archipelago_vectorize = {key: pd.DataFrame(index=X.index, dtype=float,
                                                   columns=[f"{key} Center {i + 1}" for i in
                                                            range(self.archipelago_[key].quantiser.n_clusters)])
                                 for key in self.archipelago_.keys()}
        for key, dgms in X.items():
            if key not in self.archipelago_.keys():
                continue
            archipelago_vectorize[key].loc[X.index] = self.archipelago_[key].transform(dgms)
        return pd.concat(archipelago_vectorize.values(), axis=1)

    def get_feature_names_out(self):
        # hack to get access to skl `set_output`
        pass
