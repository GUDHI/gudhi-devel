import numpy as np
from sklearn.metrics import pairwise_distances
import os
import gudhi as gd
import sys
sys.path.append("../sktda/")
from clustering import *

X = np.loadtxt("human")

print("Mapper computation with point cloud")
mapper = MapperComplex(inp="point cloud", 
                       filters=X[:,[2,0]], 
                       filter_bnds=np.array([[np.nan,np.nan],[np.nan,np.nan]]), 
                       resolutions=np.array([np.nan,np.nan]), gains=np.array([0.33,0.33]), 
                       colors=X[:,2:3],
                       ).fit(X)

f = open("mapper_pc", "w")
f.write("%s\n%s\n%s\n%f %f\n%d %d\n" % ("human", "coord2-0", "coord2", 10, 0.33, len(mapper.mapper_.get_skeleton(0)), len([edge for (edge,f) in mapper.mapper_.get_skeleton(1) if len(edge)==2])))
for (vertex,_) in mapper.mapper_.get_skeleton(0):
    f.write(str(vertex[0]) + " " + str(mapper.node_info_[vertex[0]]["colors"][0]) + " " + str(mapper.node_info_[vertex[0]]["size"]) + "\n")
for (edge,_) in mapper.mapper_.get_skeleton(1):
    if len(edge) == 2:
        f.write(str(edge[0]) + " " + str(edge[1]) + "\n")
f.close()

os.system("python3 ~/Git/gudhi-devel/src/Nerve_GIC/utilities/KeplerMapperVisuFromTxtFile.py -f ~/Git/gudhi-devel/src/cython/example/mapper_pc")
os.system("rm ~/Git/gudhi-devel/src/cython/example/mapper_pc")

dgms = mapper.compute_persistence_diagrams()
plot = gd.plot_persistence_diagram(dgms[0])
plot.show()

distrib = mapper.compute_distribution(X, N=10)
print("Distance threshold associated to confidence 90 percent is " + str(distrib[int(np.floor(0.9 * len(distrib)))]))

print("Mapper computation with pairwise distances only")
X = pairwise_distances(X)
mapper = MapperComplex(inp="distance matrix", 
                       filters=X[:,[2,0]], 
                       filter_bnds=np.array([[np.nan,np.nan],[np.nan,np.nan]]), 
                       resolutions=np.array([np.nan,np.nan]), gains=np.array([0.33,0.33]), 
                       colors=np.max(X, axis=1)[:,np.newaxis],
                       ).fit(X)

f = open("mapper_dm", "w")
f.write("%s\n%s\n%s\n%f %f\n%d %d\n" % ("human", "coord2-0", "coord2", 10, 0.33, len(mapper.mapper_.get_skeleton(0)), len([edge for (edge,f) in mapper.mapper_.get_skeleton(1) if len(edge)==2])))
for (vertex,_) in mapper.mapper_.get_skeleton(0):
    f.write(str(vertex[0]) + " " + str(mapper.node_info_[vertex[0]]["colors"][0]) + " " + str(mapper.node_info_[vertex[0]]["size"]) + "\n")
for (edge,_) in mapper.mapper_.get_skeleton(1):
    if len(edge) == 2:
        f.write(str(edge[0]) + " " + str(edge[1]) + "\n")
f.close()

os.system("python3 ~/Git/gudhi-devel/src/Nerve_GIC/utilities/KeplerMapperVisuFromTxtFile.py -f ~/Git/gudhi-devel/src/cython/example/mapper_dm")
os.system("rm ~/Git/gudhi-devel/src/cython/example/mapper_dm")

dgms = mapper.compute_persistence_diagrams()
plot = gd.plot_persistence_diagram(dgms[0])
plot.show()

distrib = mapper.compute_distribution(X, N=10)
print("Distance threshold associated to confidence 90 percent is " + str(distrib[int(np.floor(0.9 * len(distrib)))]))
