from .. import CubicalComplex
from sklearn.base import TransformerMixin

class CubicalPersistence(TransformerMixin):
    # Fast way to find primes and should be enough
    available_primes_ = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]
    """
    This is a class for computing the persistence diagrams from a cubical complex.
    """
    def __init__(self, dimensions=None, persistence_dim=0, min_persistence=0):
        """
        Constructor for the CubicalPersistence class.

        Parameters:
            dimensions (list of int): A list of number of top dimensional cells.
            persistence_dim (int): The returned persistence diagrams dimension. Default value is `0`.
            min_persistence (float): The minimum persistence value to take into account (strictly greater than
                `min_persistence`). Default value is `0.0`. Sets `min_persistence` to `-1.0` to see all values.
        """
        self.dimensions_ = dimensions
        self.persistence_dim_ = persistence_dim

        self.homology_coeff_field_ = None
        for dim in self.available_primes_:
            if dim > persistence_dim + 1:
                self.homology_coeff_field_ = dim
                break
        if self.homology_coeff_field_ == None:
            raise ValueError("persistence_dim must be less than 96")

        self.min_persistence_ = min_persistence

    def transform(self, X):
        """
        Compute all the cubical complexes and their persistence diagrams.

        Parameters:
            X (list of double OR numpy.ndarray): Cells filtration values.

        Returns:
            Persistence diagrams
        """
        cubical_complex = CubicalComplex(top_dimensional_cells = X,
                                         dimensions = self.dimensions_)
        cubical_complex.compute_persistence(homology_coeff_field = self.homology_coeff_field_,
                                            min_persistence = self.min_persistence_)
        self.diagrams_ = cubical_complex.persistence_intervals_in_dimension(self.persistence_dim_)
        if self.persistence_dim_ == 0:
            # return all but the last, always [ 0., inf]
            self.diagrams_ = self.diagrams_[:-1]
        return self.diagrams_

    def fit_transform(self, X):
        """
        Compute all the cubical complexes and their persistence diagrams.

        Parameters:
            X (list of double OR numpy.ndarray): Cells filtration values.

        Returns:
            Persistence diagrams
        """
        self.transform(X)
        return self.diagrams_
