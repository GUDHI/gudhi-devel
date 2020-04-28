#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
from sklearn.kernel_approximation import RBFSampler
from sklearn.preprocessing import MinMaxScaler

from gudhi.representations import DiagramSelector, Clamping, Landscape, Silhouette, BettiCurve, ComplexPolynomial,\
  TopologicalVector, DiagramScaler, BirthPersistenceTransform,\
  PersistenceImage, PersistenceWeightedGaussianKernel, Entropy, \
  PersistenceScaleSpaceKernel, SlicedWassersteinDistance,\
  SlicedWassersteinKernel, BottleneckDistance, PersistenceFisherKernel

D = np.array([[0.,4.],[1.,2.],[3.,8.],[6.,8.], [0., np.inf], [5., np.inf]])
diags = [D]

diags = DiagramSelector(use=True, point_type="finite").fit_transform(diags)
diags = DiagramScaler(use=True, scalers=[([0,1], MinMaxScaler())]).fit_transform(diags)
diags = DiagramScaler(use=True, scalers=[([1], Clamping(maximum=.9))]).fit_transform(diags)

D = diags[0]
plt.scatter(D[:,0],D[:,1])
plt.plot([0.,1.],[0.,1.])
plt.title("Test Persistence Diagram for vector methods")
plt.show()

LS = Landscape(resolution=1000)
L = LS.fit_transform(diags)
plt.plot(L[0][:1000])
plt.plot(L[0][1000:2000])
plt.plot(L[0][2000:3000])
plt.title("Landscape")
plt.show()

def pow(n):
  return lambda x: np.power(x[1]-x[0],n)

SH = Silhouette(resolution=1000, weight=pow(2))
sh = SH.fit_transform(diags)
plt.plot(sh[0])
plt.title("Silhouette")
plt.show()

BC = BettiCurve(resolution=1000)
bc = BC.fit_transform(diags)
plt.plot(bc[0])
plt.title("Betti Curve")
plt.show()

CP = ComplexPolynomial(threshold=-1, polynomial_type="T")
cp = CP.fit_transform(diags)
print("Complex polynomial is " + str(cp[0,:]))

TV = TopologicalVector(threshold=-1)
tv = TV.fit_transform(diags)
print("Topological vector is " + str(tv[0,:]))

PI = PersistenceImage(bandwidth=.1, weight=lambda x: x[1], im_range=[0,1,0,1], resolution=[100,100])
pi = PI.fit_transform(diags)
plt.imshow(np.flip(np.reshape(pi[0], [100,100]), 0))
plt.title("Persistence Image")
plt.show()

ET = Entropy(mode="scalar")
et = ET.fit_transform(diags)
print("Entropy statistic is " + str(et[0,:]))

ET = Entropy(mode="vector", normalized=False)
et = ET.fit_transform(diags)
plt.plot(et[0])
plt.title("Entropy function")
plt.show()

D = np.array([[1.,5.],[3.,6.],[2.,7.]])
diags2 = [D]

diags2 = DiagramScaler(use=True, scalers=[([0,1], MinMaxScaler())]).fit_transform(diags2)

D = diags[0]
plt.scatter(D[:,0],D[:,1])
D = diags2[0]
plt.scatter(D[:,0],D[:,1])
plt.plot([0.,1.],[0.,1.])
plt.title("Test Persistence Diagrams for kernel methods")
plt.show()

def arctan(C,p):
  return lambda x: C*np.arctan(np.power(x[1], p))

PWG = PersistenceWeightedGaussianKernel(bandwidth=1., kernel_approx=None, weight=arctan(1.,1.))
X = PWG.fit(diags)
Y = PWG.transform(diags2)
print("PWG kernel is " + str(Y[0][0]))

PWG = PersistenceWeightedGaussianKernel(kernel_approx=RBFSampler(gamma=1./2, n_components=100000).fit(np.ones([1,2])), weight=arctan(1.,1.))
X = PWG.fit(diags)
Y = PWG.transform(diags2)
print("Approximate PWG kernel is " + str(Y[0][0]))

PSS = PersistenceScaleSpaceKernel(bandwidth=1.)
X = PSS.fit(diags)
Y = PSS.transform(diags2)
print("PSS kernel is " + str(Y[0][0]))

PSS = PersistenceScaleSpaceKernel(kernel_approx=RBFSampler(gamma=1./2, n_components=100000).fit(np.ones([1,2])))
X = PSS.fit(diags)
Y = PSS.transform(diags2)
print("Approximate PSS kernel is " + str(Y[0][0]))

sW = SlicedWassersteinDistance(num_directions=100)
X = sW.fit(diags)
Y = sW.transform(diags2)
print("SW distance is " + str(Y[0][0]))

SW = SlicedWassersteinKernel(num_directions=100, bandwidth=1.)
X = SW.fit(diags)
Y = SW.transform(diags2)
print("SW kernel is " + str(Y[0][0]))

W = WassersteinDistance(order=2, internal_p=2, mode="pot")
X = W.fit(diags)
Y = W.transform(diags2)
print("Wasserstein distance (POT) is " + str(Y[0][0]))

W = WassersteinDistance(order=2, internal_p=2, mode="hera", delta=0.0001)
X = W.fit(diags)
Y = W.transform(diags2)
print("Wasserstein distance (hera) is " + str(Y[0][0]))

W = BottleneckDistance(epsilon=.001)
X = W.fit(diags)
Y = W.transform(diags2)
print("Bottleneck distance is " + str(Y[0][0]))

PF = PersistenceFisherKernel(bandwidth_fisher=1., bandwidth=1.)
X = PF.fit(diags)
Y = PF.transform(diags2)
print("PF kernel is " + str(Y[0][0]))

PF = PersistenceFisherKernel(bandwidth_fisher=1., bandwidth=1., kernel_approx=RBFSampler(gamma=1./2, n_components=100000).fit(np.ones([1,2])))
X = PF.fit(diags)
Y = PF.transform(diags2)
print("Approximate PF kernel is " + str(Y[0][0]))
