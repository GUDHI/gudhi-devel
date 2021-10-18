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
    pairs = st.persistence_pairs()

    # Retrieve vertices v_a and v_b by picking the ones achieving the maximal
    # distance among all pairwise distances between the simplex vertices
    indices, pers = [], []
    for s1, s2 in pairs:
        if len(s1) == dimension+1 and len(s2) > 0:
            l1, l2 = np.array(s1), np.array(s2)
            i1 = [l1[v] for v in np.unravel_index(np.argmax(DX[l1,:][:,l1]),[len(l1), len(l1)])]
            i2 = [l2[v] for v in np.unravel_index(np.argmax(DX[l2,:][:,l2]),[len(l2), len(l2)])]
            indices.append(i1)
            indices.append(i2)
            pers.append(st.filtration(s2)-st.filtration(s1))
    
    # Sort points with distance-to-diagonal
    perm = np.argsort(pers)
    indices = np.reshape(indices, [-1,4])[perm][::-1,:].flatten()

    return np.array(indices, dtype=np.int32)

class RipsLayer(tf.keras.layers.Layer):
    """
    TensorFlow layer for computing Rips persistence out of a point cloud

    Attributes:
        maximum_edge_length (float): maximum edge length for the Rips complex 
        dimension (int): homology dimension
    """
    def __init__(self, maximum_edge_length=12, dimension=1, **kwargs):
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
            indices = tf.reshape(indices, [-1,2])[1::2,:]
            dgm = tf.concat([tf.zeros([indices.shape[0],1]), tf.reshape(tf.gather_nd(DX, indices), [-1,1])], axis=1)
        return dgm

