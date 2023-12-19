# coding: utf-8
"""
@author: Martin Royer
"""

import copy

import numpy as np
import pandas as pd

from sklearn.base import BaseEstimator, TransformerMixin

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

    # def __init__(self, homology_dimensions=None, settler_class: type = Atol, settler_params=None):
    def __init__(
            self,
            homology_dimensions=None,
            settler=None,
            settler_list=None
    ):
        """
        Constructor for the Archipelago class that sets the future number of cluster per Atol created.

        Parameters:
            homology_dimensions (int or list of int): The targeted persistence diagrams dimension(s).
                Short circuit the use of :class:`~gudhi.representations.preprocessing.DimensionSelector` when only one
                dimension matters (in other words, when `homology_dimensions` is an int).
            settler: settler for populating atol, i.e. object to vectorize persistence diagrams in each homology
                dimensions. Must be copy.deepcopy-able. Will be ignored if settler_list is given.
            settler_list: settler list for populating atols, i.e. list of object to vectorize persistence diagrams
                in order of passed homology_dimensions.
        """
        if homology_dimensions is None:
            homology_dimensions = [0]
        self.homology_dimensions = homology_dimensions
        # @todo: rename this `island` or something
        if settler is None:
            settler = Atol()
        self.settler = settler
        self.settler_list = settler_list

        if isinstance(self.homology_dimensions, int):
            self.dim_list_ = [ self.homology_dimensions ]
        else:
            self.dim_list_ = self.homology_dimensions

        if settler_list is not None:
            if len(settler_list) != len(self.dim_list_):
                raise ValueError("settler_list must be of the same length as number of targeted homology dimensions.")
            else:
                self.settler_list_ = self.settler_list
        else:
            self.settler_list_ = [copy.deepcopy(self.settler) for i in range(len(self.dim_list_))]
        self.archipelago_ = {}

    def fit(self, X, y=None):
        """
        Calibration step: create and fit settler elements to the corresponding diagram element

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

            n_points = np.concatenate(dgms).shape[0]
            for (homology_dimension, settler) in zip(self.dim_list_, self.settler_list_):
                if f"(D{homology_dimension})" not in key:
                    continue
                print(f"[Archipelago] Fitting key {key} matching homology dimension {homology_dimension}.")
                if settler.quantiser.n_clusters > n_points:
                    print(f"[Archipelago] {key} has only {n_points} points for `fit`, reducing n_clusters.")
                    settler.quantiser.n_clusters = np.min([settler.quantiser.n_clusters, n_points])
                try:
                    settler.fit(X=dgms, y=y)
                except ValueError as ve:
                    print(f'[Archipelago] Fit of {key} returned "{ve}", archipelago will ignore this key.')
                    continue
                self.archipelago_[key] = settler
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
