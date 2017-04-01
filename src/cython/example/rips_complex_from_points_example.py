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
print("RipsComplex creation from points")
rips = gudhi.RipsComplex(points=[[0, 0], [1, 0], [0, 1], [1, 1]],
                         max_edge_length=42)

simplex_tree = rips.create_simplex_tree(max_dimension=1)

print("filtered_tree=", simplex_tree.get_filtered_tree())
print("star([0])=", simplex_tree.get_star([0]))
print("coface([0], 1)=", simplex_tree.get_cofaces([0], 1))
