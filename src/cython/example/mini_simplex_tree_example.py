#!/usr/bin/env python

import gudhi

"""This file is part of the Gudhi Library. The Gudhi library
   (Geometric Understanding in Higher Dimensions) is a generic C++
   library for computational topology.

   Author(s):       Vincent Rouvreau

   Copyright (C) 2016  INRIA Saclay (France)

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
__copyright__ = "Copyright (C) 2016  INRIA Saclay (France)"
__license__ = "GPL v3"

print("#####################################################################")
print("MiniSimplexTree creation from insertion")

""" Complex to build.
     1   3
     o---o
    /X\ /
   o---o   o
   2   0   4
"""

triangle012 = [0, 1, 2]
edge03 = [0, 3]
edge13 = [1, 3]
vertex4 = [4]
mini_st = gudhi.MiniSimplexTree()
mini_st.insert(triangle012)
mini_st.insert(edge03)
mini_st.insert(edge13)
mini_st.insert(vertex4)

# FIXME: Remove this line
mini_st.set_dimension(2)

# initialize_filtration required before plain_homology
mini_st.initialize_filtration()

print("persistence(homology_coeff_field=2)=")
print(mini_st.persistence(homology_coeff_field=2))

edge02 = [0, 2]
if mini_st.find(edge02):
    # Only coface is 012
    print("coface(edge02,1)=", mini_st.get_coface_tree(edge02, 1))

if mini_st.get_coface_tree(triangle012, 1) == []:
    # Precondition: Check the simplex has no coface before removing it.
    mini_st.remove_maximal_simplex(triangle012)

# initialize_filtration required after removing
mini_st.initialize_filtration()

print("filtered_tree after triangle012 removal =", mini_st.get_filtered_tree())
