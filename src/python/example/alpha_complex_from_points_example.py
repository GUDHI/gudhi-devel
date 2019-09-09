#!/usr/bin/env python

from gudhi import AlphaComplex, SimplexTree

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
print("filtrations=", simplex_tree.get_filtration())
print("star([0])=", simplex_tree.get_star([0]))
print("coface([0], 1)=", simplex_tree.get_cofaces([0], 1))

print("point[0]=", alpha_complex.get_point(0))
print("point[5]=", alpha_complex.get_point(5))
