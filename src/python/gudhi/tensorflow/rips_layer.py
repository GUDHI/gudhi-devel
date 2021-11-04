import numpy               as np
import tensorflow          as tf
from ..rips_complex     import RipsComplex

############################
# Vietoris-Rips filtration #
############################

# The parameters of the model are the point coordinates.

def _Rips(DX, max_edge, dimension):
    # Parameters: DX (distance matrix), 
    #             max_edge (maximum edge length for Rips filtration), 
    #             dimension (homology dimension)

    # Compute the persistence pairs with Gudhi
    rc = RipsComplex(distance_matrix=DX, max_edge_length=max_edge)
    st = rc.create_simplex_tree(max_dimension=dimension+1)
    dgm = st.persistence()
    if dimension == 0:
        pairs = st.flag_persistence_generators()[0]
    else:
        pairs = st.flag_persistence_generators()[1][dimension-1]

    indices = pairs.flatten()
    return np.array(indices, dtype=np.int32)

class RipsLayer(tf.keras.layers.Layer):
    """
    TensorFlow layer for computing Rips persistence out of a point cloud
    """
    def __init__(self, maximum_edge_length=12, dimension=1, **kwargs):
        """
        Constructor for the RipsLayer class

        Parameters:
            maximum_edge_length (float): maximum edge length for the Rips complex 
            dimension (int): homology dimension
        """
        super().__init__(dynamic=True, **kwargs)
        self.max_edge = maximum_edge_length
        self.dimension = dimension

    def build(self):
        super.build()
        
    def call(self, X):
        """
        Compute Rips persistence diagram associated to a point cloud

        Parameters:   
            X (TensorFlow variable): point cloud of shape [number of points, number of dimensions]

        Returns:
            dgm (TensorFlow variable): Rips persistence diagram with shape [num_points, 2] with points sorted by 
        """    
        # Compute distance matrix
        DX = tf.math.sqrt(tf.reduce_sum((tf.expand_dims(X, 1)-tf.expand_dims(X, 0))**2, 2))
        # Compute vertices associated to positive and negative simplices 
        # Don't compute gradient for this operation
        indices = tf.stop_gradient(_Rips(DX.numpy(), self.max_edge, self.dimension))
        # Get persistence diagram by simply picking the corresponding entries in the distance matrix
        if self.dimension > 0:
            dgm = tf.reshape(tf.gather_nd(DX, tf.reshape(indices, [-1,2])), [-1,2])
        else:
            #indices = tf.reshape(indices, [-1,2])[1::2,:]
            indices = indices[:,1:]
            dgm = tf.concat([tf.zeros([indices.shape[0],1]), tf.reshape(tf.gather_nd(DX, indices), [-1,1])], axis=1)
        return dgm

