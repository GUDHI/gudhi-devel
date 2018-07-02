#!/usr/bin/env python

import gudhi

"""This file is part of the Gudhi Library. The Gudhi library
   (Geometric Understanding in Higher Dimensions) is a generic C++
   library for computational topology.

   Author(s):       Vincent Rouvreau

   Copyright (C) 2016 Inria

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
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "GPL v3"

print("#####################################################################")
print("SimplexTree creation from insertion")

st = gudhi.SimplexTree()

if st.insert([0, 1]):
    print("Inserted !!")
else:
    print("Not inserted...")

if st.find([0, 1]):
    print("Found !!")
else:
    print("Not found...")

if st.insert([0, 1, 2], filtration=4.0):
    print("Inserted !!")
else:
    print("Not inserted...")

print("dimension=", st.dimension())

st.initialize_filtration()
print("filtration=", st.get_filtration())
print("filtration[1, 2]=", st.filtration([1, 2]))
print("filtration[4, 2]=", st.filtration([4, 2]))

print("num_simplices=", st.num_simplices())
print("num_vertices=", st.num_vertices())

print("skeleton[2]=", st.get_skeleton(2))
print("skeleton[1]=", st.get_skeleton(1))
print("skeleton[0]=", st.get_skeleton(0))
