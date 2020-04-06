# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Raphaël Tinarrage and Yuichi Ike
#
# Copyright (C) 2020 Inria, Copyright (C) 2020 FUjitsu Laboratories Ltd.
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

from gudhi import SimplexTree

class WeightedRipsComplex:
    """
    class to generate a weighted Rips complex 
    from a distance matrix and filtration value
    """
    def __init__(self, 
                distance_matrix=None, 
                filtration_values=None,
                max_filtration=float('inf'), sparse=None):
        self.distance_matrix = distance_matrix
        self.filtration_values = filtration_values
        self.max_filtration = max_filtration
            
    def create_simplex_tree(self, max_dimension):
        dist = self.distance_matrix
        F = self.filtration_values
        num_pts = len(dist)
        
        st = SimplexTree()
        
        for i in range(num_pts):
            if F[i] < self.max_filtration:
                st.insert([i], F[i])
        for i in range(num_pts):
            for j in range(num_pts):
                value = (dist[i][j] + F[i] + F[j]) / 2
                if value < self.max_filtration:
                    st.insert([i,j], filtration=value)
        return st
        
