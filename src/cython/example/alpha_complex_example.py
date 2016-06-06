#!/usr/bin/env python

import gudhi

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
alpha_complex = gudhi.AlphaComplex(points=[[0, 0], [1, 0], [0, 1], [1, 1]],
                                   max_alpha_square=60.0)

if alpha_complex.find([0, 1]):
    print("[0, 1] Found !!")
else:
    print("[0, 1] Not found...")

if alpha_complex.find([4]):
    print("[4] Found !!")
else:
    print("[4] Not found...")

if alpha_complex.insert([0, 1, 2], filtration=4.0):
    print("[0, 1, 2] Inserted !!")
else:
    print("[0, 1, 2] Not inserted...")

if alpha_complex.insert([0, 1, 4], filtration=4.0):
    print("[0, 1, 4] Inserted !!")
else:
    print("[0, 1, 4] Not inserted...")

if alpha_complex.find([4]):
    print("[4] Found !!")
else:
    print("[4] Not found...")

print("dimension=", alpha_complex.dimension())
print("filtered_tree=", alpha_complex.get_filtered_tree())
print("star([0])=", alpha_complex.get_star_tree([0]))
print("coface([0], 1)=", alpha_complex.get_coface_tree([0], 1))

print("point[0]=", alpha_complex.get_point(0))
print("point[5]=", alpha_complex.get_point(5))

print("betti_numbers()=")
print(alpha_complex.betti_numbers())

alpha_complex.initialize_filtration()
print("persistence(homology_coeff_field=2, min_persistence=0)=")
print(alpha_complex.persistence(homology_coeff_field=2, min_persistence=0))

print("betti_numbers()=")
print(alpha_complex.betti_numbers())
