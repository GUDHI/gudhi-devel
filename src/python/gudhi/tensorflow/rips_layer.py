import numpy               as np
import tensorflow          as tf
from ..rips_complex     import RipsComplex

############################
# Vietoris-Rips filtration #
############################

# The parameters of the model are the point coordinates.

def _Rips(DX, max_edge, dimensions):
    # Parameters: DX (distance matrix), 
    #             max_edge (maximum edge length for Rips filtration), 
    #             dimensions (homology dimensions)

    # Compute the persistence pairs with Gudhi
    rc = RipsComplex(distance_matrix=DX, max_edge_length=max_edge)
    st = rc.create_simplex_tree(max_dimension=max(dimensions)+1)
    st.compute_persistence()
    pairs = st.flag_persistence_generators()

    L_indices = []
    for dimension in dimensions:

        if dimension == 0:
            finite_pairs = pairs[0]
            essential_pairs = pairs[2]
        else:
            finite_pairs = pairs[1][dimension-1] if len(pairs[1]) >= dimension else np.empty(shape=[0,4])
            essential_pairs = pairs[3][dimension-1] if len(pairs[3]) >= dimension else np.empty(shape=[0,2])
        
        finite_indices = np.array(finite_pairs.flatten(), dtype=np.int32)
        essential_indices = np.array(essential_pairs.flatten(), dtype=np.int32)

        L_indices.append((finite_indices, essential_indices))

    return L_indices

class RipsLayer(tf.keras.layers.Layer):
    """
    TensorFlow layer for computing Rips persistence out of a point cloud
    """
    def __init__(self, dimensions, maximum_edge_length=np.inf, min_persistence=None, **kwargs):
        """
        Constructor for the RipsLayer class

        Parameters:
            maximum_edge_length (float): maximum edge length for the Rips complex 
            dimensions (List[int]): homology dimensions
            min_persistence (List[float]): minimum distance-to-diagonal of the points in the output persistence diagrams (default None, in which case 0. is used for all dimensions)
        """
        super().__init__(dynamic=True, **kwargs)
        self.max_edge = maximum_edge_length
        self.dimensions = dimensions
        self.min_persistence = min_persistence if min_persistence != None else [0. for _ in range(len(self.dimensions))]
        assert len(self.min_persistence) == len(self.dimensions)
        
    def call(self, X):
        """
        Compute Rips persistence diagram associated to a point cloud

        Parameters:   
            X (TensorFlow variable): point cloud of shape [number of points, number of dimensions]

        Returns:
            dgms (list of tuple of TensorFlow variables): list of Rips persistence diagrams of length self.dimensions, where each element of the list is a tuple that contains the finite and essential persistence diagrams of shapes [num_finite_points, 2] and [num_essential_points, 1] respectively
        """    
        # Compute distance matrix
        DX = tf.norm(tf.expand_dims(X, 1)-tf.expand_dims(X, 0), axis=2)
        # Compute vertices associated to positive and negative simplices 
        # Don't compute gradient for this operation
        indices = _Rips(DX.numpy(), self.max_edge, self.dimensions)
        # Get persistence diagrams by simply picking the corresponding entries in the distance matrix
        self.dgms = []
        for idx_dim, dimension in enumerate(self.dimensions):
            cur_idx = indices[idx_dim]
            if dimension > 0:
                finite_dgm = tf.reshape(tf.gather_nd(DX, tf.reshape(cur_idx[0], [-1,2])), [-1,2])
                essential_dgm = tf.reshape(tf.gather_nd(DX, tf.reshape(cur_idx[1], [-1,2])), [-1,1])
            else:
                reshaped_cur_idx = tf.reshape(cur_idx[0], [-1,3])
                finite_dgm = tf.concat([tf.zeros([reshaped_cur_idx.shape[0],1]), tf.reshape(tf.gather_nd(DX, reshaped_cur_idx[:,1:]), [-1,1])], axis=1)
                essential_dgm = tf.zeros([cur_idx[1].shape[0],1])
            min_pers = self.min_persistence[idx_dim]
            if min_pers >= 0:
                persistent_indices = np.argwhere(np.abs(finite_dgm[:,1]-finite_dgm[:,0]) > min_pers).ravel()
                self.dgms.append((tf.gather(finite_dgm, indices=persistent_indices), essential_dgm))
            else:
                self.dgms.append((finite_dgm, essential_dgm))
        return self.dgms
