import numpy               as np
import tensorflow          as tf
from ..cubical_complex  import CubicalComplex

######################
# Cubical filtration #
######################

# The parameters of the model are the pixel values.

def _Cubical(Xflat, Xdim, dimensions, homology_coeff_field):
    # Parameters: Xflat (flattened image),
    #             Xdim (shape of non-flattened image)
    #             dimensions (homology dimensions)

    # Compute the persistence pairs with Gudhi
    # We reverse the dimensions because CubicalComplex uses Fortran ordering
    cc = CubicalComplex(dimensions=Xdim[::-1], top_dimensional_cells=Xflat)
    cc.compute_persistence(homology_coeff_field=homology_coeff_field)

    # Retrieve and output image indices/pixels corresponding to positive and negative simplices
    cof_pp = cc.cofaces_of_persistence_pairs()
    
    L_cofs = []
    for dim in dimensions:

        try:
            cof = cof_pp[0][dim]
        except IndexError:
            cof = np.array([])

        L_cofs.append(np.array(cof, dtype=np.int32))

    return L_cofs

class CubicalLayer(tf.keras.layers.Layer):
    """
    TensorFlow layer for computing the persistent homology of a cubical complex
    """
    def __init__(self, homology_dimensions, min_persistence=None, homology_coeff_field=11, **kwargs):
        """
        Constructor for the CubicalLayer class

        Parameters:
            homology_dimensions (List[int]): list of homology dimensions
            min_persistence (List[float]): minimum distance-to-diagonal of the points in the output persistence diagrams (default None, in which case 0. is used for all dimensions)
            homology_coeff_field (int): homology field coefficient. Must be a prime number. Default value is 11. Max is 46337.
        """
        super().__init__(**kwargs)
        self.dimensions = homology_dimensions
        self.min_persistence = min_persistence if min_persistence is not None else [0.] * len(self.dimensions)
        self.hcf = homology_coeff_field
        assert len(self.min_persistence) == len(self.dimensions)

    def call(self, X):
        """
        Compute persistence diagram associated to a cubical complex filtered by some pixel values 

        Parameters:
            X (TensorFlow variable): pixel values of the cubical complex

        Returns:
            List[Tuple[tf.Tensor,tf.Tensor]]: List of cubical persistence diagrams. The length of this list is the same than that of dimensions, i.e., there is one persistence diagram per homology dimension provided in the input list dimensions. Moreover, the finite and essential parts of the persistence diagrams are provided separately: each element of this list is a tuple of size two that contains the finite and essential parts of the corresponding persistence diagram, of shapes [num_finite_points, 2] and [num_essential_points, 1] respectively. Note that the essential part is always empty in cubical persistence diagrams, except in homology dimension zero, where the essential part always contains a single point, with abscissa equal to the smallest value in the complex, and infinite ordinate
        """
        # Compute pixels associated to positive and negative simplices 
        # Don't compute gradient for this operation
        Xflat = tf.reshape(X, [-1])
        Xdim, Xflat_numpy = X.shape, Xflat.numpy()
        indices_list = _Cubical(Xflat_numpy, Xdim, self.dimensions, self.hcf)
        index_essential = np.argmin(Xflat_numpy) # index of minimum pixel value for essential persistence diagram
        # Get persistence diagram by simply picking the corresponding entries in the image
        self.dgms = []
        for idx_dim, dimension in enumerate(self.dimensions):
            finite_dgm = tf.reshape(tf.gather(Xflat, indices_list[idx_dim]), [-1,2])
            essential_dgm = tf.reshape(tf.gather(Xflat, index_essential), [-1,1]) if dimension == 0 else tf.zeros([0, 1])
            min_pers = self.min_persistence[idx_dim]
            if min_pers >= 0:
                persistent_indices = tf.where(tf.math.abs(finite_dgm[:,1]-finite_dgm[:,0]) > min_pers)
                self.dgms.append((tf.reshape(tf.gather(finite_dgm, indices=persistent_indices), [-1,2]), essential_dgm))
            else:
                self.dgms.append((finite_dgm, essential_dgm))
        return self.dgms
