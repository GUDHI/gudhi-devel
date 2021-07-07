import numpy               as np
import tensorflow          as tf
from ..rips_complex     import RipsComplex
from ..cubical_complex  import CubicalComplex

# In this file, we write functions based on the Gudhi library that compute persistence diagrams associated to 
# different filtrations (lower star, Rips, cubical), as well as the corresponding positive and negative 
# simplices. We also wrap these functions into Tensorflow models.



#########################################
# Lower star filtration on simplex tree #
#########################################

# The parameters of the model are the vertex function values of the simplex tree.

def _LowerStarSimplexTree(simplextree, filtration, dimension):
    # Parameters: simplextree (simplex tree on which to compute persistence)
    #             filtration (function values on the vertices of st),
    #             dimension (homology dimension),
    
    for s,_ in simplextree.get_filtration():
        simplextree.assign_filtration(s, -1e10)

    # Assign new filtration values
    for i in range(simplextree.num_vertices()):
        simplextree.assign_filtration([i], filtration[i])
    simplextree.make_filtration_non_decreasing()
    
    # Compute persistence diagram
    dgm = simplextree.persistence()
    
    # Get vertex pairs for optimization. First, get all simplex pairs
    pairs = simplextree.persistence_pairs()
    
    # Then, loop over all simplex pairs
    indices, pers = [], []
    for s1, s2 in pairs:
        # Select pairs with good homological dimension and finite lifetime
        if len(s1) == dimension+1 and len(s2) > 0:
            # Get IDs of the vertices corresponding to the filtration values of the simplices
            l1, l2 = np.array(s1), np.array(s2)
            i1 = l1[np.argmax(filtration[l1])]
            i2 = l2[np.argmax(filtration[l2])]
            indices.append(i1)
            indices.append(i2)
            # Compute lifetime
            pers.append(simplextree.filtration(s2)-simplextree.filtration(s1))
    
    # Sort vertex pairs wrt lifetime
    perm = np.argsort(pers)
    indices = np.reshape(indices, [-1,2])[perm][::-1,:].flatten()
    
    return np.array(indices, dtype=np.int32)

class LowerStarSimplexTreeLayer(tf.keras.layers.Layer):
    """
    TensorFlow layer for computing lower-star persistence out of a simplex tree

    Attributes:
        simplextree (gudhi.SimplexTree()): underlying simplex tree 
        dimension (int): homology dimension
    """
    def __init__(self, simplextree, dimension=0, **kwargs):
        super().__init__(dynamic=True, **kwargs)
        self.dimension   = dimension
        self.simplextree = simplextree
    
    def build(self):
        super.build()
    
    def call(self, filtration):
        """
        Compute lower-star persistence diagram associated to a function defined on the vertices of the simplex tree

        Parameters:
            F (TensorFlow variable): filter function values over the vertices of the simplex tree
        """
        # Don't try to compute gradients for the vertex pairs
        indices = tf.stop_gradient(_LowerStarSimplexTree(self.simplextree, filtration.numpy(), self.dimension))
        # Get persistence diagram
        self.dgm = tf.reshape(tf.gather(filtration, indices), [-1,2]) 
        return self.dgm









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

    if len(cof) > 0:
        # Sort points with distance-to-diagonal
        Xs = X.shape
        pers = [X[np.unravel_index(cof[idx,1], Xs)] - X[np.unravel_index(cof[idx,0], Xs)] for idx in range(len(cof))]
        perm = np.argsort(pers)
        cof = cof[perm[::-1]]
    
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

    Attributes:
        dimension (int): homology dimension
    """
    def __init__(self, dimension=1, **kwargs):
        super().__init__(dynamic=True, **kwargs)
        self.dimension = dimension

    def build(self):
        super.build()
        
    def call(self, X):
        """
        Compute persistence diagram associated to a cubical complex filtered by some pixel values 

        Parameters:
            X (TensorFlow variable): pixel values of the cubical complex
        """
        # Compute pixels associated to positive and negative simplices 
        # Don't compute gradient for this operation
        indices = tf.stop_gradient(_Cubical(X.numpy(), self.dimension))
        # Get persistence diagram by simply picking the corresponding entries in the image
        dgm = tf.reshape(tf.gather_nd(X, tf.reshape(indices, [-1,len(X.shape)])), [-1,2])
        return dgm
