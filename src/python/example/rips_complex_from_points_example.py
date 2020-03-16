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
print("RipsComplex creation from points")
rips = gudhi.RipsComplex(points=[[0, 0], [1, 0], [0, 1], [1, 1]], max_edge_length=42)

simplex_tree = rips.create_simplex_tree(max_dimension=1)

print("filtrations=")
for simplex_with_filtration in simplex_tree.get_filtration():
    print("(%s, %.2f)" % tuple(simplex_with_filtration))

print("star([0])=", simplex_tree.get_star([0]))
print("coface([0], 1)=", simplex_tree.get_cofaces([0], 1))
