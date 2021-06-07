# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2021 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

from .. import CubicalComplex
from sklearn.base import BaseEstimator, TransformerMixin

# joblib is required by scikit-learn
from joblib import Parallel, delayed


class CubicalPersistence(BaseEstimator, TransformerMixin):
    # Fast way to find primes and should be enough
    _available_primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]
    """
    This is a class for computing the persistence diagrams from a cubical complex.
    """

    def __init__(self, dimensions=None, persistence_dim=0, min_persistence=0, n_jobs=None):
        """
        Constructor for the CubicalPersistence class.

        Parameters:
            dimensions (list of int): A list of number of top dimensional cells if cells filtration values will require
                to be reshaped (cf. :func:`~gudhi.sklearn.cubical_persistence.CubicalPersistence.transform`)
            persistence_dim (int): The returned persistence diagrams dimension. Default value is `0`.
            min_persistence (float): The minimum persistence value to take into account (strictly greater than
                `min_persistence`). Default value is `0.0`. Sets `min_persistence` to `-1.0` to see all values.
            n_jobs (int): cf. https://joblib.readthedocs.io/en/latest/generated/joblib.Parallel.html
        """
        self.dimensions = dimensions
        self.persistence_dim = persistence_dim

        self.homology_coeff_field_ = None
        for dim in self._available_primes:
            if dim > persistence_dim + 1:
                self.homology_coeff_field_ = dim
                break
        if self.homology_coeff_field_ == None:
            raise ValueError("persistence_dim must be less than 96")

        self.min_persistence = min_persistence
        self.n_jobs = n_jobs

    def fit(self, X, Y=None):
        """
        Nothing to be done, but useful when included in a scikit-learn Pipeline.
        """
        return self

    def __transform(self, cells):
        cubical_complex = CubicalComplex(top_dimensional_cells=cells, dimensions=self.dimensions)
        cubical_complex.compute_persistence(
            homology_coeff_field=self.homology_coeff_field_, min_persistence=self.min_persistence
        )
        diagrams = cubical_complex.persistence_intervals_in_dimension(self.persistence_dim)
        return diagrams

    def transform(self, X, Y=None):
        """
        Compute all the cubical complexes and their associated persistence diagrams.

        Parameters:
            X (list of list of double OR list of numpy.ndarray): List of cells filtration values that can be flatten if
                dimensions is set in the constructor, or already with the correct shape in a numpy.ndarray (and
                dimensions must not be set).

        Returns:
            Persistence diagrams
        """

        # threads is preferred as cubical construction and persistence computation releases the GIL
        return Parallel(n_jobs=self.n_jobs, prefer="threads")(delayed(self.__transform)(cells) for cells in X)
