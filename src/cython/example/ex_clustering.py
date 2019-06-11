import numpy as np
import sklearn_tda as tda
from sklearn.cluster import DBSCAN
import os
from clustering import *

X = np.loadtxt("human")

map = tda.MapperComplex(filters=X[:,[2,0]], resolutions=[10,3], gains=[0.33,0.33], colors=X[:,2:3], clustering=DBSCAN(eps=0.05, min_samples=5)).fit(X)

num_pts, num_edges = 0, 0
min_col, min_sz, max_col, max_sz = 1e10, 1e10, -1e10, -1e10
for value in map.graph_:
    if len(value[0]) == 1:
        min_col = min(min_col, value[1])
        max_col = max(max_col, value[1])
        min_sz = min(min_sz, value[2])
        max_sz = max(max_sz, value[2])
        num_pts += 1
    if len(value[0]) == 2:
        num_edges += 1

f = open("mapper", "w")
f.write("%s\n%s\n%s\n%f %f\n%d %d\n" % ("human", "coord2-0", "coord2", 10, 0.33, num_pts, num_edges))
for value in map.graph_:
    if len(value[0]) == 1:
        f.write(str(value[0][0]) + " " + str(value[1][0]) + " " + str(value[2]) + "\n")
for value in map.graph_:
    if len(value[0]) == 2:
        f.write(str(value[0][0]) + " " + str(value[0][1]) + "\n")
f.close()

os.system("python3 ~/Git/gudhi-devel/src/sklearn_tda/KeplerMapperVisuFromTxtFile.py -f ~/Git/gudhi-devel/src/sklearn_tda/mapper")
os.system("firefox ~/Git/gudhi-devel/src/sklearn_tda/mapper.html")
