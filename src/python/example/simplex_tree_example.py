#!/usr/bin/env python

import gudhi

""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2016 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "MIT"

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

print("simplices=")
for simplex_with_filtration in st.get_simplices():
    print("(%s, %.2f)" % tuple(simplex_with_filtration))

st.initialize_filtration()
print("filtration=")
for simplex_with_filtration in st.get_filtration():
    print("(%s, %.2f)" % tuple(simplex_with_filtration))

print("filtration[1, 2]=", st.filtration([1, 2]))
print("filtration[4, 2]=", st.filtration([4, 2]))

print("num_simplices=", st.num_simplices())
print("num_vertices=", st.num_vertices())

print("skeleton[2]=", st.get_skeleton(2))
print("skeleton[1]=", st.get_skeleton(1))
print("skeleton[0]=", st.get_skeleton(0))
