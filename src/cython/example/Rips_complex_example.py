#!/usr/bin/env python

import gudhi

print("#######################################################################")
print("RipsComplex creation from points")
rips = gudhi.RipsComplex(points=[[0,0],[1,0],[0,1],[1,1]],max_dimension=1,max_edge_length=42)

print("filtered_tree=", rips.get_filtered_tree())
print("star([0])=", rips.get_star_tree([0]))
print("coface([0],1)=", rips.get_coface_tree([0], 1))
