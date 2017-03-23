#!/usr/bin/env python

from gudhi import AlphaComplex, SimplexTree

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
print("AlphaComplex creation from points")
alpha_complex = AlphaComplex(points=[[0, 0], [1, 0], [0, 1], [1, 1]])
simplex_tree = alpha_complex.create_simplex_tree(max_alpha_square=60.0)

if simplex_tree.find([0, 1]):
    print("[0, 1] Found !!")
else:
    print("[0, 1] Not found...")

if simplex_tree.find([4]):
    print("[4] Found !!")
else:
    print("[4] Not found...")

if simplex_tree.insert([0, 1, 2], filtration=4.0):
    print("[0, 1, 2] Inserted !!")
else:
    print("[0, 1, 2] Not inserted...")

if simplex_tree.insert([0, 1, 4], filtration=4.0):
    print("[0, 1, 4] Inserted !!")
else:
    print("[0, 1, 4] Not inserted...")

if simplex_tree.find([4]):
    print("[4] Found !!")
else:
    print("[4] Not found...")

print("dimension=", simplex_tree.dimension())
print("filtered_tree=", simplex_tree.get_filtered_tree())
print("star([0])=", simplex_tree.get_star([0]))
print("coface([0], 1)=", simplex_tree.get_cofaces([0], 1))

print("point[0]=", alpha_complex.get_point(0))
print("point[5]=", alpha_complex.get_point(5))
