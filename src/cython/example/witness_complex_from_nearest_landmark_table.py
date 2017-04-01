#!/usr/bin/env python

from gudhi import StrongWitnessComplex, SimplexTree

"""This file is part of the Gudhi Library. The Gudhi library
   (Geometric Understanding in Higher Dimensions) is a generic C++
   library for computational topology.

   Author(s):       Vincent Rouvreau

   Copyright (C) 2016 INRIA

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 INRIA"
__license__ = "GPL v3"

print("#####################################################################")
print("WitnessComplex creation from nearest landmark table")
nearest_landmark_table = [[[0, 0], [1, 1], [2, 2], [3, 3], [4, 4]],
                          [[1, 0], [2, 1], [3, 2], [4, 3], [0, 4]],
                          [[2, 0], [3, 1], [4, 2], [0, 3], [1, 4]],
                          [[3, 0], [4, 1], [0, 2], [1, 3], [2, 4]],
                          [[4, 0], [0, 1], [1, 2], [2, 3], [3, 4]]]

witness_complex = StrongWitnessComplex(nearest_landmark_table=nearest_landmark_table)
simplex_tree = witness_complex.create_simplex_tree(max_alpha_square=4.1)

message = "Number of simplices: " + repr(simplex_tree.num_simplices())
print(message)

diag = simplex_tree.persistence(min_persistence=-0.1, homology_coeff_field=11)
print(diag)
