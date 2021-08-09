# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2021 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

from sklearn.base import BaseEstimator, TransformerMixin

# joblib is required by scikit-learn
from joblib import Parallel, delayed

# Mermaid sequence diagram - https://mermaid-js.github.io/mermaid-live-editor/
# sequenceDiagram
#     USER->>DimensionSelector: fit_transform(<br/>[[array( H0(X0) ), array( H1(X0) ), ...],<br/> [array( H0(X1) ), array( H1(X1) ), ...],<br/> ...])
#     DimensionSelector->>thread1: _transform([array( H0(X0) ), array( H1(X0) )], ...)
#     DimensionSelector->>thread2: _transform([array( H0(X1) ), array( H1(X1) )], ...)
#     Note right of DimensionSelector: ...
#     thread1->>DimensionSelector: array( Hn(X0) )
#     thread2->>DimensionSelector: array( Hn(X1) )
#     Note right of DimensionSelector: ...
#     DimensionSelector->>USER: [array( Hn(X0) ), <br/> array( Hn(X1) ), <br/> ...]


class DimensionSelector(BaseEstimator, TransformerMixin):
    """
    This is a class to select persistence diagrams in a specific dimension.
    """

    def __init__(self, persistence_dimension=0, n_jobs=None):
        """
        Constructor for the DimensionSelector class.

        Parameters:
            persistence_dimension (int): The returned persistence diagrams dimension. Default value is `0`.
        """
        self.persistence_dimension = persistence_dimension
        self.n_jobs = n_jobs

    def fit(self, X, Y=None):
        """
        Nothing to be done, but useful when included in a scikit-learn Pipeline.
        """
        return self

    def transform(self, X, Y=None):
        """
        Select persistence diagrams from its dimension.

        Parameters:
            X (list of list of pairs): List of list of persistence pairs, i.e.
            `[[array( H0(X0) ), array( H1(X0) ), ...], [array( H0(X1) ), array( H1(X1) ), ...], ...]` 

        Returns:
            Persistence diagrams in a specific dimension, i.e.
            `[array( Hn(X0) ), array( Hn(X1), ...]`
        """

        return [persistence[self.persistence_dimension] for persistence in X]
