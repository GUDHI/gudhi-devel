#!/usr/bin/env python

""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Marc Glisse

    Copyright (C) 2019 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

__license__ = "GPL v3"  # Because of AlphaComplex


import numpy as np
import gudhi as gd


points = gd.read_points_from_off_file(off_file="../../data/points/tore3D_1307.off")
ac = gd.AlphaComplex(points=points)
st = ac.create_simplex_tree()
points = np.array([ac.get_point(i) for i in range(st.num_vertices())])
# We want to plot the alpha-complex with alpha=0.1.
# We are only going to plot the triangles
triangles = np.array([s[0] for s in st.get_skeleton(2) if len(s[0]) == 3 and s[1] <= 0.1])

# First possibility: plotly
import plotly.graph_objects as go

fig = go.Figure(
    data=[
        go.Mesh3d(
            x=points[:, 0],
            y=points[:, 1],
            z=points[:, 2],
            i=triangles[:, 0],
            j=triangles[:, 1],
            k=triangles[:, 2],
        )
    ]
)
fig.show()

# Second possibility: matplotlib
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.gca(projection="3d")
ax.plot_trisurf(points[:, 0], points[:, 1], points[:, 2], triangles=triangles)
plt.show()

# Third possibility: mayavi
from mayavi import mlab

mlab.triangular_mesh(points[:, 0], points[:, 1], points[:, 2], triangles)
mlab.show()
