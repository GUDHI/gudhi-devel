#!/usr/bin/env python

""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2016 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

__license__ = "MIT"


from gudhi import StrongWitnessComplex


print("#####################################################################")
print("WitnessComplex creation from nearest landmark table")
nearest_landmark_table = [
    [[0, 0.0], [1, 0.1], [2, 0.2], [3, 0.3], [4, 0.4]],
    [[1, 0.0], [2, 0.1], [3, 0.2], [4, 0.3], [0, 0.4]],
    [[2, 0.0], [3, 0.1], [4, 0.2], [0, 0.3], [1, 0.4]],
    [[3, 0.0], [4, 0.1], [0, 0.2], [1, 0.3], [2, 0.4]],
    [[4, 0.0], [0, 0.1], [1, 0.2], [2, 0.3], [3, 0.4]],
]

witness_complex = StrongWitnessComplex(nearest_landmark_table=nearest_landmark_table)
simplex_tree = witness_complex.create_simplex_tree(max_alpha_square=0.41)

print(f"Number of simplices: {simplex_tree.num_simplices()}")

diag = simplex_tree.persistence(min_persistence=-0.1, homology_coeff_field=11)
print(diag)
