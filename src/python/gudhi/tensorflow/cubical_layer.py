import numpy               as np
import tensorflow          as tf
from ..cubical_complex  import CubicalComplex

######################
# Cubical filtration #
######################

# The parameters of the model are the pixel values.

def _Cubical(X, dimensions):
    # Parameters: X (image),
    #             dimensions (homology dimensions)

    # Compute the persistence pairs with Gudhi
    Xs = X.shape
    cc = CubicalComplex(top_dimensional_cells=X)
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
    TensorFlow layer for computing cubical persistence out of a cubical complex
    """
    def __init__(self, dimensions, **kwargs):
        """
        Constructor for the CubicalLayer class

        Parameters:
            dimensions (list of int): homology dimensions
        """
        super().__init__(dynamic=True, **kwargs)
        self.dimensions = dimensions

    def build(self):
        super.build()
        
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
        indices = _Cubical(X.numpy(), self.dimensions)
        # Get persistence diagram by simply picking the corresponding entries in the image
        self.dgms = [tf.reshape(tf.gather(  tf.reshape(tf.transpose(X), [-1]), indice  ), [-1,2]) for indice in indices]
        return self.dgms
