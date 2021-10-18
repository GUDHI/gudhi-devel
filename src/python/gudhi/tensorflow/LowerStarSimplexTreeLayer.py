import numpy               as np
import tensorflow          as tf

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

