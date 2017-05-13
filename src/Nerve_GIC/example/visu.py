import km
import numpy as np
from collections import defaultdict

network = {}
mapper = km.KeplerMapper(verbose=0)
data = np.zeros((3,3))
projected_data = mapper.fit_transform( data, projection="sum", scaler=None )

f = open('SC.txt','r')
nodes = defaultdict(list)
links = defaultdict(list)
custom = defaultdict(list)

dat = f.readline()
lens = f.readline()
color = f.readline();
param = [float(i) for i in f.readline().split(" ")]

nums = [int(i) for i in f.readline().split(" ")]
num_nodes = nums[0]
num_edges = nums[1]

for i in range(0,num_nodes):
	point = [float(j) for j in f.readline().split(" ")]
	nodes[  str(int(point[0]))  ] = [  int(point[0]), point[1], int(point[2])  ]
	links[  str(int(point[0]))  ] = []
	custom[  int(point[0])  ] = point[1]

m = min([custom[i] for i in range(0,num_nodes)])
M = max([custom[i] for i in range(0,num_nodes)])

for i in range(0,num_edges):
	edge = [int(j) for j in f.readline().split(" ")]	
	links[  str(edge[0])  ].append(  str(edge[1])  )	
	links[  str(edge[1])  ].append(  str(edge[0])  )

network["nodes"] = nodes
network["links"] = links
network["meta"] = lens

mapper.visualize(network, color_function = color, path_html="SC.html", title=dat,
graph_link_distance=30, graph_gravity=0.1, graph_charge=-120, custom_tooltips=custom, width_html=0, 
height_html=0, show_tooltips=True, show_title=True, show_meta=True, res=param[0],gain=param[1], minimum=m,maximum=M)