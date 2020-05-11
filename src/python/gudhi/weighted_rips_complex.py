# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Raphaël Tinarrage, Yuichi Ike, Masatoshi Takenouchi
#
# Copyright (C) 2020 Inria, Copyright (C) 2020 FUjitsu Laboratories Ltd.
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

from gudhi import SimplexTree

class WeightedRipsComplex:
    """
    Class to generate a weighted Rips complex from a distance matrix and weights on vertices, 
    in the way described in the paper 'DTM-based filtrations' https://arxiv.org/abs/1811.04757.
    Remark that the filtration value of a vertex is twice of its weight for the consistency with 
    RipsComplex, which is different from the definition in the paper.
    """
    def __init__(self, 
                distance_matrix, 
                weights=None,
                max_filtration=float('inf')):
        """
        Args:
            distance_matrix (Sequence[Sequence[float]]): distance matrix (full square or lower triangular).
            weights (Sequence[float]): (one half of) weight for each vertex.
            max_filtration (float): specifies the maximal filtration value to be considered.      
        """
        self.distance_matrix = distance_matrix
        if weights is not None:
            self.weights = weights
        else:
            self.weights = [0] * len(distance_matrix)
        self.max_filtration = max_filtration
            
    def create_simplex_tree(self, max_dimension):
        """
        Args:
            max_dimension (int): graph expansion until this given dimension.
        """
        dist = self.distance_matrix
        F = self.weights
        num_pts = len(dist)
        
        st = SimplexTree()
        
        for i in range(num_pts):
            if 2*F[i] <= self.max_filtration:
                st.insert([i], 2*F[i])
        for i in range(num_pts):
            for j in range(i):
                value = max(2*F[i], 2*F[j], dist[i][j] + F[i] + F[j])
                # max is needed when F is not 1-Lipschitz
                if value <= self.max_filtration:
                    st.insert([i,j], filtration=value)
                    
        st.expansion(max_dimension) 
        return st
        
