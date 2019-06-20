import sys
sys.path.append("../")
from clustering import *
import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity

# ToMATo with custom parameters
X = np.loadtxt("../../../data/points/spiral_w_o_density.txt")[::10,:]
X[:,0]/=max(X[:,0])
X[:,1]/=max(X[:,1])

tom = ToMATo(tau=.3, density_estimator=KernelDensity(bandwidth=.02), n_neighbors=10, verbose=True)
tom.fit(X)
plt.scatter(X[:,0], X[:,1], s=5., c=tom.labels_)
plt.show()

# ToMATo with all parameters estimated from data 
X = np.loadtxt("../../../data/points/toy_example_w_o_density.txt")
X[:,0]/=max(X[:,0])
X[:,1]/=max(X[:,1])

tom = ToMATo(verbose=True)
tom.fit(X)
plt.scatter(X[:,0], X[:,1], s=5., c=tom.labels_)
plt.show()
