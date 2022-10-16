import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pytest
import random

from sklearn.cluster import KMeans

# Vectorization
from gudhi.representations import (Landscape, Silhouette, BettiCurve, ComplexPolynomial,\
  TopologicalVector, PersistenceImage, Entropy)

# Preprocessing
from gudhi.representations import (BirthPersistenceTransform, Clamping, DiagramScaler, Padding, ProminentPoints, \
  DiagramSelector)

# Kernel
from gudhi.representations import (PersistenceWeightedGaussianKernel, \
  PersistenceScaleSpaceKernel, SlicedWassersteinDistance,\
  SlicedWassersteinKernel, PersistenceFisherKernel, WassersteinDistance)


def test_representations_examples():
    # Disable graphics for testing purposes
    plt.show = lambda: None
    here = os.path.dirname(os.path.realpath(__file__))
    sys.path.append(here + "/../example")
    import diagram_vectorizations_distances_kernels

    return None


from gudhi.representations.vector_methods import Atol
from gudhi.representations.metrics import *
from gudhi.representations.kernel_methods import *


def _n_diags(n):
    l = []
    for _ in range(n):
        a = np.random.rand(50, 2)
        a[:, 1] += a[:, 0]  # So that y >= x
        l.append(a)
    return l


def test_multiple():
    l1 = _n_diags(9)
    l2 = _n_diags(11)
    l1b = l1.copy()
    d1 = pairwise_persistence_diagram_distances(l1, e=0.00001, n_jobs=4)
    d2 = BottleneckDistance(epsilon=0.00001).fit_transform(l1)
    d3 = pairwise_persistence_diagram_distances(l1, l1b, e=0.00001, n_jobs=4)
    assert d1 == pytest.approx(d2)
    assert d3 == pytest.approx(d2, abs=1e-5)  # Because of 0 entries (on the diagonal)
    d1 = pairwise_persistence_diagram_distances(l1, l2, metric="wasserstein", order=2, internal_p=2)
    d2 = WassersteinDistance(order=2, internal_p=2, n_jobs=4).fit(l2).transform(l1)
    print(d1.shape, d2.shape)
    assert d1 == pytest.approx(d2, rel=0.02)


# Test sorted values as points order can be inverted, and sorted test is not documentation-friendly
# Note the test below must be up to date with the Atol class documentation
def test_atol_doc():
    a = np.array([[1, 2, 4], [1, 4, 0], [1, 0, 4]])
    b = np.array([[4, 2, 0], [4, 4, 0], [4, 0, 2]])
    c = np.array([[3, 2, -1], [1, 2, -1]])

    atol_vectoriser = Atol(quantiser=KMeans(n_clusters=2, random_state=202006))
    # Atol will do
    # X = np.concatenate([a,b,c])
    # kmeans = KMeans(n_clusters=2, random_state=202006).fit(X) 
    # kmeans.labels_ will be : array([1, 0, 1, 0, 0, 1, 0, 0])
    first_cluster = np.asarray([a[0], a[2], b[2]])
    second_cluster = np.asarray([a[1], b[0], b[2], c[0], c[1]])

    # Check the center of the first_cluster and second_cluster are in Atol centers
    centers = atol_vectoriser.fit(X=[a, b, c]).centers
    np.isclose(centers, first_cluster.mean(axis=0)).all(1).any() 
    np.isclose(centers, second_cluster.mean(axis=0)).all(1).any() 

    vectorization = atol_vectoriser.transform(X=[a, b, c])
    assert np.allclose(vectorization[0], atol_vectoriser(a))
    assert np.allclose(vectorization[1], atol_vectoriser(b))
    assert np.allclose(vectorization[2], atol_vectoriser(c))


def test_dummy_atol():
    a = np.array([[1, 2, 4], [1, 4, 0], [1, 0, 4]])
    b = np.array([[4, 2, 0], [4, 4, 0], [4, 0, 2]])
    c = np.array([[3, 2, -1], [1, 2, -1]])

    for weighting_method in ["cloud", "iidproba"]:
        for contrast in ["gaussian", "laplacian", "indicator"]:
            atol_vectoriser = Atol(
                quantiser=KMeans(n_clusters=1, random_state=202006),
                weighting_method=weighting_method,
                contrast=contrast,
            )
            atol_vectoriser.fit([a, b, c])
            atol_vectoriser(a)
            atol_vectoriser.transform(X=[a, b, c])


from gudhi.representations.vector_methods import BettiCurve

def test_infinity():
    a = np.array([[1.0, 8.0], [2.0, np.inf], [3.0, 4.0]])
    c = BettiCurve(20, [0.0, 10.0])(a)
    assert c[1] == 0
    assert c[7] == 3
    assert c[9] == 2

def test_preprocessing_empty_diagrams():
    empty_diag = np.empty(shape = [0, 2])
    assert not np.any(BirthPersistenceTransform()(empty_diag))
    assert not np.any(Clamping().fit_transform(empty_diag))
    assert not np.any(DiagramScaler()(empty_diag))
    assert not np.any(Padding()(empty_diag))
    assert not np.any(ProminentPoints()(empty_diag))
    assert not np.any(DiagramSelector()(empty_diag))

def pow(n):
  return lambda x: np.power(x[1]-x[0],n)

def test_vectorization_empty_diagrams():
    empty_diag = np.empty(shape = [0, 2])
    random_resolution = random.randint(50,100)*10 # between 500 and 1000
    print("resolution = ", random_resolution)
    lsc = Landscape(resolution=random_resolution)(empty_diag)
    assert not np.any(lsc)
    assert lsc.shape[0]%random_resolution == 0
    slt = Silhouette(resolution=random_resolution, weight=pow(2))(empty_diag)
    assert not np.any(slt)
    assert slt.shape[0] == random_resolution
    btc = BettiCurve(resolution=random_resolution)(empty_diag)
    assert not np.any(btc)
    assert btc.shape[0] == random_resolution
    cpp = ComplexPolynomial(threshold=random_resolution, polynomial_type="T")(empty_diag)
    assert not np.any(cpp)
    assert cpp.shape[0] == random_resolution
    tpv = TopologicalVector(threshold=random_resolution)(empty_diag)
    assert tpv.shape[0] == random_resolution
    assert not np.any(tpv)
    prmg = PersistenceImage(resolution=[random_resolution,random_resolution])(empty_diag)
    assert not np.any(prmg)
    assert prmg.shape[0] == random_resolution * random_resolution
    sce = Entropy(mode="scalar", resolution=random_resolution)(empty_diag)
    assert not np.any(sce)
    assert sce.shape[0] == 1
    scv = Entropy(mode="vector", normalized=False, resolution=random_resolution)(empty_diag)
    assert not np.any(scv)
    assert scv.shape[0] == random_resolution
    
def test_entropy_miscalculation():
    diag_ex = np.array([[0.0,1.0], [0.0,1.0], [0.0,2.0]])
    def pe(pd):
        l = pd[:,1] - pd[:,0]
        l = l/sum(l)
        return -np.dot(l, np.log(l))
    sce = Entropy(mode="scalar")
    assert [[pe(diag_ex)]] == sce.fit_transform([diag_ex])
    sce = Entropy(mode="vector", resolution=4, normalized=False)
    pef = [-1/4*np.log(1/4)-1/4*np.log(1/4)-1/2*np.log(1/2),
           -1/4*np.log(1/4)-1/4*np.log(1/4)-1/2*np.log(1/2),
           -1/2*np.log(1/2), 
           0.0]
    assert all(([pef] == sce.fit_transform([diag_ex]))[0])
    sce = Entropy(mode="vector", resolution=4, normalized=True)
    pefN = (sce.fit_transform([diag_ex]))[0]
    area = np.linalg.norm(pefN, ord=1)
    assert area==1
        
def test_kernel_empty_diagrams():
    empty_diag = np.empty(shape = [0, 2])
    assert SlicedWassersteinDistance(num_directions=100)(empty_diag, empty_diag) == 0.
    assert SlicedWassersteinKernel(num_directions=100, bandwidth=1.)(empty_diag, empty_diag) == 1.
    assert WassersteinDistance(mode="hera", delta=0.0001)(empty_diag, empty_diag) == 0.
    assert WassersteinDistance(mode="pot")(empty_diag, empty_diag) == 0.
    assert BottleneckDistance(epsilon=.001)(empty_diag, empty_diag) == 0.
    assert BottleneckDistance()(empty_diag, empty_diag) == 0.
#    PersistenceWeightedGaussianKernel(bandwidth=1., kernel_approx=None, weight=arctan(1.,1.))(empty_diag, empty_diag)
#    PersistenceWeightedGaussianKernel(kernel_approx=RBFSampler(gamma=1./2, n_components=100000).fit(np.ones([1,2])), weight=arctan(1.,1.))(empty_diag, empty_diag)
#    PersistenceScaleSpaceKernel(bandwidth=1.)(empty_diag, empty_diag)
#    PersistenceScaleSpaceKernel(kernel_approx=RBFSampler(gamma=1./2, n_components=100000).fit(np.ones([1,2])))(empty_diag, empty_diag)
#    PersistenceFisherKernel(bandwidth_fisher=1., bandwidth=1.)(empty_diag, empty_diag)
#    PersistenceFisherKernel(bandwidth_fisher=1., bandwidth=1., kernel_approx=RBFSampler(gamma=1./2, n_components=100000).fit(np.ones([1,2])))(empty_diag, empty_diag)

