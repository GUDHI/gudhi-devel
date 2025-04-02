#!/usr/bin/env python

""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       ???

    Copyright (C) 20?? Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

__author__ = "???"
__maintainer__ = ""
__copyright__ = "Copyright (C) 20?? Inria"
__license__ = "MIT"


import matplotlib.pyplot as plt
import numpy as np
from sklearn.kernel_approximation import RBFSampler
from sklearn.preprocessing import MinMaxScaler

from gudhi.representations import (
    DiagramSelector,
    Clamping,
    Landscape,
    Silhouette,
    BettiCurve,
    ComplexPolynomial,
    TopologicalVector,
    DiagramScaler,
    BirthPersistenceTransform,
    PersistenceImage,
    PersistenceWeightedGaussianKernel,
    Entropy,
    PersistenceScaleSpaceKernel,
    SlicedWassersteinDistance,
    SlicedWassersteinKernel,
    PersistenceFisherKernel,
    WassersteinDistance,
)


D1 = np.array([[0.0, 4.0], [1.0, 2.0], [3.0, 8.0], [6.0, 8.0], [0.0, np.inf], [5.0, np.inf]])

proc1 = DiagramSelector(use=True, point_type="finite")
proc2 = DiagramScaler(use=True, scalers=[([0, 1], MinMaxScaler())])
proc3 = DiagramScaler(use=True, scalers=[([1], Clamping(maximum=0.9))])
D1 = proc3(proc2(proc1(D1)))

plt.scatter(D1[:, 0], D1[:, 1])
plt.plot([0.0, 1.0], [0.0, 1.0])
plt.title("Test Persistence Diagram for vector methods")
plt.show()

LS = Landscape(resolution=1000)
L = LS(D1)
plt.plot(L[:1000])
plt.plot(L[1000:2000])
plt.plot(L[2000:3000])
plt.title("Landscape")
plt.show()


def pow(n):
    return lambda x: np.power(x[1] - x[0], n)


SH = Silhouette(resolution=1000, weight=pow(2))
plt.plot(SH(D1))
plt.title("Silhouette")
plt.show()

BC = BettiCurve(resolution=1000)
plt.plot(BC(D1))
plt.title("Betti Curve")
plt.show()

CP = ComplexPolynomial(threshold=-1, polynomial_type="T")
print("Complex polynomial is " + str(CP(D1)))

TV = TopologicalVector(threshold=-1)
print("Topological vector is " + str(TV(D1)))

PI = PersistenceImage(
    bandwidth=0.1, weight=lambda x: x[1], im_range=[0, 1, 0, 1], resolution=[100, 100]
)
plt.imshow(np.flip(np.reshape(PI(D1), [100, 100]), 0))
plt.title("Persistence Image")
plt.show()

ET = Entropy(mode="scalar")
print("Entropy statistic is " + str(ET(D1)))

ET = Entropy(mode="vector", normalized=False)
plt.plot(ET(D1))
plt.title("Entropy function")
plt.show()

D2 = np.array([[1.0, 5.0], [3.0, 6.0], [2.0, 7.0]])
D2 = proc3(proc2(proc1(D2)))

plt.scatter(D1[:, 0], D1[:, 1])
plt.scatter(D2[:, 0], D2[:, 1])
plt.plot([0.0, 1.0], [0.0, 1.0])
plt.title("Test Persistence Diagrams for kernel methods")
plt.show()


def arctan(C, p):
    return lambda x: C * np.arctan(np.power(x[1], p))


PWG = PersistenceWeightedGaussianKernel(
    bandwidth=1.0, kernel_approx=None, weight=arctan(1.0, 1.0)
)
print("PWG kernel is " + str(PWG(D1, D2)))

PWG = PersistenceWeightedGaussianKernel(
    kernel_approx=RBFSampler(gamma=1.0 / 2, n_components=100000).fit(np.ones([1, 2])),
    weight=arctan(1.0, 1.0),
)
print("Approximate PWG kernel is " + str(PWG(D1, D2)))

PSS = PersistenceScaleSpaceKernel(bandwidth=1.0)
print("PSS kernel is " + str(PSS(D1, D2)))

PSS = PersistenceScaleSpaceKernel(
    kernel_approx=RBFSampler(gamma=1.0 / 2, n_components=100000).fit(np.ones([1, 2]))
)
print("Approximate PSS kernel is " + str(PSS(D1, D2)))

sW = SlicedWassersteinDistance(num_directions=100)
print("SW distance is " + str(sW(D1, D2)))

SW = SlicedWassersteinKernel(num_directions=100, bandwidth=1.0)
print("SW kernel is " + str(SW(D1, D2)))

try:
    W = WassersteinDistance(order=2, internal_p=2, mode="pot")
    print("Wasserstein distance (POT) is " + str(W(D1, D2)))
except ImportError:
    print("WassersteinDistance (POT) is not available, you may be missing pot.")

W = WassersteinDistance(order=2, internal_p=2, mode="hera", delta=0.0001)
print("Wasserstein distance (hera) is " + str(W(D1, D2)))

try:
    from gudhi.representations import BottleneckDistance

    W = BottleneckDistance(epsilon=0.001)
    print("Bottleneck distance is " + str(W(D1, D2)))
except ImportError:
    print("BottleneckDistance is not available, you may be missing CGAL.")

PF = PersistenceFisherKernel(bandwidth_fisher=1.0, bandwidth=1.0)
print("PF kernel is " + str(PF(D1, D2)))

PF = PersistenceFisherKernel(
    bandwidth_fisher=1.0,
    bandwidth=1.0,
    kernel_approx=RBFSampler(gamma=1.0 / 2, n_components=100000).fit(np.ones([1, 2])),
)
print("Approximate PF kernel is " + str(PF(D1, D2)))
