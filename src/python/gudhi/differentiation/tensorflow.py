import numpy               as np
import tensorflow          as tf
import tensorflow_addons   as tfa
from ..simplex_tree     import SimplexTree
from ..rips_complex     import RipsComplex
from ..cubical_complex  import CubicalComplex

# In this file, we write functions based on the Gudhi library that compute persistence diagrams associated to 
# different filtrations (lower star, Rips, cubical), as well as the corresponding positive and negative 
# simplices. We also wrap these functions into Tensorflow models.



#########################################
# Lower star filtration on simplex tree #
#########################################

# The parameters of the model are the vertex function values of the simplex tree.

def _SimplexTree(stbase, fct, dim, card):
    # Parameters: stbase (array containing the name of the file where the simplex tree is located)
    #             fct (function values on the vertices of stbase),
    #             dim (homological dimension),
    #             card (number of persistence diagram points, sorted by distance-to-diagonal)
    
    # Copy stbase in another simplex tree st
    st = SimplexTree()
    f = open(stbase[0], "r")
    for line in f:
        ints = line.split(" ")
        s = [int(v) for v in ints[:-1]]
        st.insert(s, -1e10)
    f.close()
        
    # Assign new filtration values
    for i in range(st.num_vertices()):
        st.assign_filtration([i], fct[i])
    st.make_filtration_non_decreasing()
    
    # Compute persistence diagram
    dgm = st.persistence()
    
    # Get vertex pairs for optimization. First, get all simplex pairs
    pairs = st.persistence_pairs()
    
    # Then, loop over all simplex pairs
    indices, pers = [], []
    for s1, s2 in pairs:
        # Select pairs with good homological dimension and finite lifetime
        if len(s1) == dim+1 and len(s2) > 0:
            # Get IDs of the vertices corresponding to the filtration values of the simplices
            l1, l2 = np.array(s1), np.array(s2)
            i1 = l1[np.argmax(fct[l1])]
            i2 = l2[np.argmax(fct[l2])]
            indices.append(i1)
            indices.append(i2)
            # Compute lifetime
            pers.append(st.filtration(s2) - st.filtration(s1))
    
    # Sort vertex pairs wrt lifetime
    perm = np.argsort(pers)
    indices = list(np.reshape(indices, [-1,2])[perm][::-1,:].flatten())
    
    # Pad vertex pairs
    indices = indices[:2*card] + [0 for _ in range(0,max(0,2*card-len(indices)))]
    return list(np.array(indices, dtype=np.int32))

class LowerStarSimplexTreeModel(tf.keras.Model):
    """
    TensorFlow model for computing lower-star persistence out of a simplex tree. Since simplex trees cannot be easily encoded as TensorFlow variables, the model takes as input a path to a file containing the simplex tree simplices, and read it each time the simplex tree is required for computations.

    Attributes:
        F (TensorFlow variable): filter function values over the vertices of the simplex tree
        stbase (string): path to the file containing the simplex tree. Each line of the file should represent a simplex as a sequence of integers separated by spaces
        card (int): maximum number of points in the persistence diagram
        dim (int): homology dimension
    """
    def __init__(self, F, stbase="simplextree.txt", dim=0, card=50):
        super(SimplexTreeModel, self).__init__()
        self.F = F
        self.dim = dim
        self.card = card
        self.st = stbase
        
    def call(self):
        d, c = self.dim, self.card
        st, fct = self.st, self.F

        # Turn STPers into a numpy function
        SimplexTreeTF = lambda fct: tf.numpy_function(_SimplexTree, [np.array([st], dtype=str), fct, d, c], [tf.int32 for _ in range(2*c)])
        
        # Don't try to compute gradients for the vertex pairs
        fcts = tf.reshape(fct, [1, self.F.shape[0]])
        inds = tf.nest.map_structure(tf.stop_gradient, tf.map_fn(SimplexTreeTF, 
                                                                 fcts, dtype=[tf.int32 for _ in range(2*c)]))
        
        # Get persistence diagram
        self.dgm = tf.reshape(tf.gather_nd(self.F, inds), [c,2]) 
        return self.dgm










############################
# Vietoris-Rips filtration #
############################

# The parameters of the model are the point coordinates.

def _Rips(DX, mel, dim, card):
    # Parameters: DX (distance matrix), 
    #             mel (maximum edge length for Rips filtration), 
    #             dim (homological dimension), 
    #             card (number of persistence diagram points, sorted by distance-to-diagonal)

    # Compute the persistence pairs with Gudhi
    rc = RipsComplex(distance_matrix=DX, max_edge_length=mel)
    st = rc.create_simplex_tree(max_dimension=dim+1)
    dgm = st.persistence()
    pairs = st.persistence_pairs()

    # Retrieve vertices v_a and v_b by picking the ones achieving the maximal
    # distance among all pairwise distances between the simplex vertices
    indices, pers = [], []
    for s1, s2 in pairs:
        if len(s1) == dim+1 and len(s2) > 0:
            l1, l2 = np.array(s1), np.array(s2)
            i1 = [s1[v] for v in np.unravel_index(np.argmax(DX[l1,:][:,l1]),[len(s1), len(s1)])]
            i2 = [s2[v] for v in np.unravel_index(np.argmax(DX[l2,:][:,l2]),[len(s2), len(s2)])]
            indices += i1
            indices += i2
            pers.append(st.filtration(s2) - st.filtration(s1))
    
    # Sort points with distance-to-diagonal
    perm = np.argsort(pers)
    indices = list(np.reshape(indices, [-1,4])[perm][::-1,:].flatten())
    
    # Output indices
    indices = indices[:4*card] + [0 for _ in range(0,max(0,4*card-len(indices)))]
    return list(np.array(indices, dtype=np.int32))

class RipsModel(tf.keras.Model):
    """
    TensorFlow model for computing Rips persistence out of a point cloud.

    Attributes:
        X (TensorFlow variable): point cloud of shape [number of points, number of dimensions]
        mel (float): maximum edge length for the Rips complex 
        card (int): maximum number of points in the persistence diagram
        dim (int): homology dimension
    """
    def __init__(self, X, mel=12, dim=1, card=50):
        super(RipsModel, self).__init__()
        self.X = X
        self.mel = mel
        self.dim = dim
        self.card = card
        
    def call(self):
        m, d, c = self.mel, self.dim, self.card
        
        # Compute distance matrix
        DX = tfa.losses.metric_learning.pairwise_distance(self.X)
        DXX = tf.reshape(DX, [1, DX.shape[0], DX.shape[1]])
        
        # Turn numpy function into tensorflow function
        RipsTF = lambda DX: tf.numpy_function(_Rips, [DX, m, d, c], [tf.int32 for _ in range(4*c)])
        
        # Compute vertices associated to positive and negative simplices 
        # Don't compute gradient for this operation
        ids = tf.nest.map_structure(tf.stop_gradient, tf.map_fn(RipsTF,DXX,dtype=[tf.int32 for _ in range(4*c)]))
        
        # Get persistence diagram by simply picking the corresponding entries in the distance matrix
        dgm = tf.reshape(tf.gather_nd(DX, tf.reshape(ids, [2*c,2])), [c,2])
        return dgm









######################
# Cubical filtration #
######################

# The parameters of the model are the pixel values.

def _Cubical(X, dim, card):
    # Parameters: X (image),
    #             dim (homological dimension), 
    #             card (number of persistence diagram points, sorted by distance-to-diagonal)

    # Compute the persistence pairs with Gudhi
    cc = CubicalComplex(dimensions=X.shape, top_dimensional_cells=X.flatten())
    cc.persistence()
    cof = cc.cofaces_of_persistence_pairs()[0][dim]

    # Sort points with distance-to-diagonal
    Xs = X.shape
    pers = [X[np.unravel_index(cof[idx,1], Xs)] - X[np.unravel_index(cof[idx,0], Xs)] for idx in range(len(cof))]
    perm = np.argsort(pers)
    cof = cof[perm[::-1]]
    
    # Retrieve and ouput image indices/pixels corresponding to positive and negative simplices
    D = len(Xs)
    ocof = np.array([0 for _ in range(D*card*2)])
    count = 0
    for idx in range(0,min(2*card, 2*cof.shape[0]),2):
        ocof[D*idx:D*(idx+1)]     = np.unravel_index(cof[count,0], Xs)
        ocof[D*(idx+1):D*(idx+2)] = np.unravel_index(cof[count,1], Xs)
        count += 1
    return list(np.array(ocof, dtype=np.int32))

class CubicalModel(tf.keras.Model):
    """
    TensorFlow model for computing cubical persistence out of a cubical complex.

    Attributes:
        X (TensorFlow variable): pixel values of the cubical complex
        card (int): maximum number of points in the persistence diagram
        dim (int): homology dimension
    """
    def __init__(self, X, dim=1, card=50):
        super(CubicalModel, self).__init__()
        self.X = X
        self.dim = dim
        self.card = card
        
    def call(self):
        d, c, D = self.dim, self.card, len(self.X.shape)
        XX = tf.reshape(self.X, [1, self.X.shape[0], self.X.shape[1]])
        
        # Turn numpy function into tensorflow function
        CbTF = lambda X: tf.numpy_function(_Cubical, [X, d, c], [tf.int32 for _ in range(2*D*c)])
        
        # Compute pixels associated to positive and negative simplices 
        # Don't compute gradient for this operation
        inds = tf.nest.map_structure(tf.stop_gradient, tf.map_fn(CbTF,XX,dtype=[tf.int32 for _ in range(2*D*c)]))
        
        # Get persistence diagram by simply picking the corresponding entries in the image
        dgm = tf.reshape(tf.gather_nd(self.X, tf.reshape(inds, [-1,D])), [-1,2])
        return dgm
