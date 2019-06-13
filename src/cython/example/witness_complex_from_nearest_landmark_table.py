#!/usr/bin/env python

from gudhi import StrongWitnessComplex, SimplexTree

"""This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

   Copyright (C) 2016 Inria

   Modification(s):
     - YYYY/MM Author: Description of the modification
"""

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "GPL v3"

print("#####################################################################")
print("WitnessComplex creation from nearest landmark table")
nearest_landmark_table = [[[0, 0.0], [1, 0.1], [2, 0.2], [3, 0.3], [4, 0.4]],
                          [[1, 0.0], [2, 0.1], [3, 0.2], [4, 0.3], [0, 0.4]],
                          [[2, 0.0], [3, 0.1], [4, 0.2], [0, 0.3], [1, 0.4]],
                          [[3, 0.0], [4, 0.1], [0, 0.2], [1, 0.3], [2, 0.4]],
                          [[4, 0.0], [0, 0.1], [1, 0.2], [2, 0.3], [3, 0.4]]]

witness_complex = StrongWitnessComplex(nearest_landmark_table=nearest_landmark_table)
simplex_tree = witness_complex.create_simplex_tree(max_alpha_square=0.41)

message = "Number of simplices: " + repr(simplex_tree.num_simplices())
print(message)

diag = simplex_tree.persistence(min_persistence=-0.1, homology_coeff_field=11)
print(diag)
