import numpy               as np
import tensorflow          as tf
from ..cubical_complex  import CubicalComplex

######################
# Cubical filtration #
######################

# The parameters of the model are the pixel values.

def _Cubical(Xflat, Xdim, dimensions):
    # Parameters: Xflat (flattened image),
    #             Xdim (shape of non-flattened image)
    #             dimensions (homology dimensions)

    # Compute the persistence pairs with Gudhi
    # We reverse the dimensions because CubicalComplex uses Fortran ordering
    cc = CubicalComplex(dimensions=Xdim[::-1], top_dimensional_cells=Xflat)
    cc.compute_persistence()

    # Retrieve and ouput image indices/pixels corresponding to positive and negative simplices    
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
    def __init__(self, dimensions, min_persistence=None, **kwargs):
        """
        Constructor for the CubicalLayer class

        Parameters:
            dimensions (List[int]): homology dimensions
            min_persistence (List[float]): minimum distance-to-diagonal of the points in the output persistence diagrams (default None, in which case 0. is used for all dimensions)
        """
        super().__init__(dynamic=True, **kwargs)
        self.dimensions = dimensions
        self.min_persistence = min_persistence if min_persistence != None else [0.] * len(self.dimensions)
        assert len(self.min_persistence) == len(self.dimensions)

    def call(self, X):
        """
        Compute persistence diagram associated to a cubical complex filtered by some pixel values 

        Parameters:
            X (TensorFlow variable): pixel values of the cubical complex

        Returns:
            dgms (list of TensorFlow variables): list of cubical persistence diagrams of length self.dimensions, where each element contains a finite persistence diagram of shape [num_finite_points, 2]
        """
        # Compute pixels associated to positive and negative simplices 
        # Don't compute gradient for this operation
        Xflat = tf.reshape(X, [-1])
        Xdim = X.shape
        indices_list = _Cubical(Xflat.numpy(), Xdim, self.dimensions)
        # Get persistence diagram by simply picking the corresponding entries in the image
        self.dgms = [tf.reshape(tf.gather(Xflat, indices), [-1,2]) for indices in indices_list]
        for idx_dim in range(len(self.min_persistence)):
            min_pers = self.min_persistence[idx_dim]
            if min_pers >= 0:
                finite_dgm = self.dgms[idx_dim]
                persistent_indices = tf.where(tf.math.abs(finite_dgm[:,1]-finite_dgm[:,0]) > min_pers)
                self.dgms[idx_dim] = tf.reshape(tf.gather(finite_dgm, indices=persistent_indices), [-1,2])
        return self.dgms
