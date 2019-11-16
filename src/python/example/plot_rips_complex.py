#!/usr/bin/env python

import numpy as np
import gudhi
points = np.array(gudhi.read_off('../../data/points/Kl.off'))
rc = gudhi.RipsComplex(points=points, max_edge_length=.2)
st = rc.create_simplex_tree(max_dimension=2)
# We are only going to plot the triangles
triangles = np.array([s[0] for s in st.get_skeleton(2) if len(s[0])==3])

# First possibility: plotly
import plotly.graph_objects as go
fig = go.Figure(data=[
    go.Mesh3d(
        # Use the first 3 coordinates, but we could as easily pick others
        x=points[:,0],
        y=points[:,1],
        z=points[:,2],
        i = triangles[:,0],
        j = triangles[:,1],
        k = triangles[:,2],
    )
])
fig.show()

# Second possibility: matplotlib
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_trisurf(points[:,0], points[:,1], points[:,2], triangles=triangles)
plt.show()

# Third possibility: mayavi.mlab.triangular_mesh
