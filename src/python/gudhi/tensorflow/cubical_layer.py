import numpy               as np
import tensorflow          as tf
from ..cubical_complex  import CubicalComplex

######################
# Cubical filtration #
######################

# The parameters of the model are the pixel values.

def _Cubical(X, dimension):
    # Parameters: X (image),
    #             dimension (homology dimension)

    # Compute the persistence pairs with Gudhi
    cc = CubicalComplex(dimensions=X.shape, top_dimensional_cells=X.flatten())
    cc.persistence()
    try:
        cof = cc.cofaces_of_persistence_pairs()[0][dimension]
    except IndexError:
        cof = np.array([])

    # Retrieve and ouput image indices/pixels corresponding to positive and negative simplices
    D = len(Xs) if len(cof) > 0 else 1
    ocof = np.array([0 for _ in range(D*2*cof.shape[0])])
    count = 0
    for idx in range(0,2*cof.shape[0],2):
        ocof[D*idx:D*(idx+1)]     = np.unravel_index(cof[count,0], Xs)
        ocof[D*(idx+1):D*(idx+2)] = np.unravel_index(cof[count,1], Xs)
        count += 1
    return np.array(ocof, dtype=np.int32)

class CubicalLayer(tf.keras.layers.Layer):
    """
    TensorFlow layer for computing cubical persistence out of a cubical complex
    """
    def __init__(self, dimension=1, **kwargs):
        """
        Constructor for the CubicalLayer class

        Parameters:
            dimension (int): homology dimension
        """
        super().__init__(dynamic=True, **kwargs)
        self.dimension = dimension

    def build(self):
        super.build()
        
    def call(self, X):
        """
        Compute persistence diagram associated to a cubical complex filtered by some pixel values 

        Parameters:
            X (TensorFlow variable): pixel values of the cubical complex

        Returns:
            dgm (TensorFlow variable): cubical persistence diagram with shape [num_points, 2]
        """
        # Compute pixels associated to positive and negative simplices 
        # Don't compute gradient for this operation
        indices = tf.stop_gradient(_Cubical(X.numpy(), self.dimension))
        # Get persistence diagram by simply picking the corresponding entries in the image
        dgm = tf.reshape(tf.gather_nd(X, tf.reshape(indices, [-1,len(X.shape)])), [-1,2])
        return dgm
