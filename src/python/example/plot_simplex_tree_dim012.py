#!/usr/bin/env python
import numpy as np
import gudhi

# Coordinates of the points
points=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,1],[1,1,0],[0,1,1]])
# Build the simplicial complex with a tetrahedon, an edge and an isolated vertex
cplx=gudhi.SimplexTree()
cplx.insert([1,2,3,5])
cplx.insert([4,6])
cplx.insert([0])
# List of triangles (point indices)
triangles = np.array([s[0] for s in cplx.get_skeleton(2) if len(s[0])==3])
# List of edges (point coordinates)
edges = []
for s in cplx.get_skeleton(1):
    e = s[0]
    if len(e) == 2:
        edges.append(points[[e[0],e[1]]])

## With plotly
import plotly.graph_objects as go
# Plot triangles
f2 = go.Mesh3d(
        x=points[:,0],
        y=points[:,1],
        z=points[:,2],
        i = triangles[:,0],
        j = triangles[:,1],
        k = triangles[:,2],
    )
# Plot points
f0 = go.Scatter3d(x=points[:,0], y=points[:,1], z=points[:,2], mode="markers")
data = [f2, f0]
# Plot edges
for pts in edges:
    seg = go.Scatter3d(x=pts[:,0],y=pts[:,1],z=pts[:,2],mode="lines",line=dict(color='green'))
    data.append(seg)
fig = go.Figure(data=data,layout=dict(showlegend=False))
# By default plotly would give each edge its own color and legend, that's too much
fig.show()

## With matplotlib
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.gca(projection='3d')
# Plot triangles
ax.plot_trisurf(points[:,0], points[:,1], points[:,2], triangles=triangles)
# Plot points
ax.scatter3D(points[:,0], points[:,1], points[:,2])
# Plot edges
ax.add_collection3d(Line3DCollection(segments=edges))
plt.show()

## With mayavi
from mayavi import mlab
# Plot triangles
mlab.triangular_mesh(points[:,0], points[:,1], points[:,2], triangles);
# Plot points
mlab.points3d(points[:,0], points[:,1], points[:,2])
# Plot edges
for pts in edges:
    mlab.plot3d(pts[:,0],pts[:,1],pts[:,2],tube_radius=None)
mlab.show()
