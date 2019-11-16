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

from mayavi import mlab
# Plot triangles
triangles = np.array([s[0] for s in cplx.get_skeleton(2) if len(s[0])==3])
mlab.triangular_mesh(points[:,0], points[:,1], points[:,2], triangles);
# Plot edges
for s in cplx.get_skeleton(1):
    e = s[0]
    if len(e) == 2:
        pts = points[[e[0],e[1]]]
        mlab.plot3d(pts[:,0],pts[:,1],pts[:,2],tube_radius=None)
# Plot points
mlab.points3d(points[:,0], points[:,1], points[:,2])
mlab.show()
