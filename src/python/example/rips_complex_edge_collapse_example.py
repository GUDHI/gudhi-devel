#!/usr/bin/env python

""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2016 Inria

    Modification(s):
      - 2025/06 Vincent Rouvreau: Use datasets generators
      - YYYY/MM Author: Description of the modification
"""

__license__ = "MIT"


import matplotlib.pyplot as plt
import time
import gudhi as gd
from gudhi.datasets.generators import points


print("################################################################################")
print("RipsComplex (only the one-skeleton) creation from 300 random points on a 2-torus")

# Set the seed for reproductibility
points.seed(0)
pts = points.c3dtorus(n_samples = 300, R = 2, r = 1)
rips_complex = gd.RipsComplex(points=pts, max_edge_length=12.)
simplex_tree = rips_complex.create_simplex_tree(max_dimension=1)
print(f"1. Rips complex has {simplex_tree.num_simplices()} simplices - {simplex_tree.num_vertices()} vertices.")

# Expansion of this one-skeleton would require a lot of memory. Let's collapse it
start = time.process_time()
simplex_tree.collapse_edges()
print(f"2. Rips complex has {simplex_tree.num_simplices()} simplices - {simplex_tree.num_vertices()} vertices.")

simplex_tree.expansion(3)
diag = simplex_tree.persistence()
print(f"Collapse, expansion and persistence computation took {time.process_time() - start} sec.")

# Use subplots to display diagram and density side by side
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))
gd.plot_persistence_diagram(diag, axes=axes[0])
axes[0].set_title("Persistence after 1 collapse")

# Collapse can be performed several times. Let's collapse it 3 times
start = time.process_time()
simplex_tree.collapse_edges(nb_iterations=3)
print(f"3. Rips complex has {simplex_tree.num_simplices()} simplices - {simplex_tree.num_vertices()} vertices.")

simplex_tree.expansion(3)
diag = simplex_tree.persistence()
print(f"Collapse, expansion and persistence computation took {time.process_time() - start} sec.")

gd.plot_persistence_diagram(diag, axes=axes[1])
axes[1].set_title("Persistence after 3 more collapses")

# Plot the 2 persistence diagrams side to side to check the persistence is the same
plt.show()
