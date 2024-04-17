# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Mathieu Carrière
#
# Copyright (C) 2018-2019 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.metrics import pairwise_distances, pairwise_kernels
from .metrics import SlicedWassersteinDistance, PersistenceFisherDistance, _sklearn_wrapper, _pairwise, pairwise_persistence_diagram_distances, _sliced_wasserstein_distance, _persistence_fisher_distance
from .preprocessing import Padding

#############################################
# Kernel methods ############################
#############################################

def _persistence_weighted_gaussian_kernel(D1, D2, weight=lambda x: 1, kernel_approx=None, bandwidth=1.):
    """
    This is a function for computing the persistence weighted Gaussian kernel value from two persistence diagrams. The persistence weighted Gaussian kernel is computed by convolving the persistence diagram points with weighted Gaussian kernels. See http://proceedings.mlr.press/v48/kusano16.html for more details.

    Parameters:
        D1: (n x 2) numpy.array encoding the (finite points of the) first diagram. Must not contain essential points (i.e. with infinite coordinate).
        D2: (m x 2) numpy.array encoding the second diagram.
        bandwidth (double): bandwidth of the Gaussian kernel with which persistence diagrams will be convolved
        weight: weight function for the persistence diagram points (default constant function, ie lambda x: 1). This function must be defined on 2D points, ie lists or numpy arrays of the form [p_x,p_y].
        kernel_approx: kernel approximation class used to speed up computation. Common kernel approximations classes can be found in the scikit-learn library (such as RBFSampler for instance).

    Returns:
        float: the persistence weighted Gaussian kernel value between persistence diagrams. 
    """
    ws1 = np.array([weight(D1[j,:]) for j in range(len(D1))])
    ws2 = np.array([weight(D2[j,:]) for j in range(len(D2))])
    if kernel_approx is not None:
        approx1 = np.sum(np.multiply(ws1[:,np.newaxis], kernel_approx.transform(D1)), axis=0)
        approx2 = np.sum(np.multiply(ws2[:,np.newaxis], kernel_approx.transform(D2)), axis=0)
        return (1./(np.sqrt(2*np.pi)*bandwidth)) * np.matmul(approx1, approx2.T)
    else:
        W = np.matmul(ws1[:,np.newaxis], ws2[np.newaxis,:])
        E = (1./(np.sqrt(2*np.pi)*bandwidth)) * np.exp(-np.square(pairwise_distances(D1,D2))/(2*bandwidth*bandwidth))
        return np.sum(np.multiply(W, E))

def _persistence_scale_space_kernel(D1, D2, kernel_approx=None, bandwidth=1.):
    """
    This is a function for computing the persistence scale space kernel value from two persistence diagrams. The persistence scale space kernel is computed by adding the symmetric to the diagonal of each point in each persistence diagram, with negative weight, and then convolving the points with a Gaussian kernel. See https://www.cv-foundation.org/openaccess/content_cvpr_2015/papers/Reininghaus_A_Stable_Multi-Scale_2015_CVPR_paper.pdf for more details.
    
    Parameters:
        D1: (n x 2) numpy.array encoding the (finite points of the) first diagram. Must not contain essential points (i.e. with infinite coordinate).
        D2: (m x 2) numpy.array encoding the second diagram.
        bandwidth (double): bandwidth of the Gaussian kernel with which persistence diagrams will be convolved
        kernel_approx: kernel approximation class used to speed up computation. Common kernel approximations classes can be found in the scikit-learn library (such as RBFSampler for instance).
    
    Returns:
        float: the persistence scale space kernel value between persistence diagrams. 
    """
    DD1 = np.concatenate([D1, D1[:,[1,0]]], axis=0)
    DD2 = np.concatenate([D2, D2[:,[1,0]]], axis=0)
    weight_pss = lambda x: 1 if x[1] >= x[0] else -1
    return 0.5 * _persistence_weighted_gaussian_kernel(DD1, DD2, weight=weight_pss, kernel_approx=kernel_approx, bandwidth=bandwidth)


def pairwise_persistence_diagram_kernels(X, Y=None, kernel="sliced_wasserstein", n_jobs=None, **kwargs):
    """
    This function computes the kernel matrix between two lists of persistence diagrams given as numpy arrays of shape (nx2).

    Parameters:    
        X (list of n numpy arrays of shape (numx2)): first list of persistence diagrams. 
        Y (list of m numpy arrays of shape (numx2)): second list of persistence diagrams (optional). If None, pairwise kernel values are computed from the first list only.
        kernel: kernel to use. It can be either a string ("sliced_wasserstein", "persistence_scale_space", "persistence_weighted_gaussian", "persistence_fisher") or a function taking two numpy arrays of shape (nx2) and (mx2) as inputs. If it is a function, make sure that it is symmetric.
        n_jobs (int): number of jobs to use for the computation. This uses joblib.Parallel(prefer="threads"), so kernels that do not release the GIL may not scale unless run inside a `joblib.parallel_backend <https://joblib.readthedocs.io/en/latest/parallel.html#joblib.parallel_backend>`_ block.
        **kwargs: optional keyword parameters. Any further parameters are passed directly to the kernel function. See the docs of the various kernel classes in this module.

    Returns: 
        numpy array of shape (nxm): kernel matrix.
    """    
    XX = np.reshape(np.arange(len(X)), [-1,1])
    YY = None if Y is None or Y is X else np.reshape(np.arange(len(Y)), [-1,1])
    if kernel == "sliced_wasserstein":
        return np.exp(-pairwise_persistence_diagram_distances(X, Y, metric="sliced_wasserstein", num_directions=kwargs["num_directions"], n_jobs=n_jobs) / kwargs["bandwidth"])
    elif kernel == "persistence_fisher":
        return np.exp(-pairwise_persistence_diagram_distances(X, Y, metric="persistence_fisher", kernel_approx=kwargs["kernel_approx"], bandwidth=kwargs["bandwidth_fisher"], n_jobs=n_jobs) / kwargs["bandwidth"])
    elif kernel == "persistence_scale_space":
        return _pairwise(pairwise_kernels, False, XX, YY, metric=_sklearn_wrapper(_persistence_scale_space_kernel, X, Y, **kwargs), n_jobs=n_jobs)
    elif kernel == "persistence_weighted_gaussian":
        return _pairwise(pairwise_kernels, False, XX, YY, metric=_sklearn_wrapper(_persistence_weighted_gaussian_kernel, X, Y, **kwargs), n_jobs=n_jobs)
    else:
        return _pairwise(pairwise_kernels, False, XX, YY, metric=_sklearn_wrapper(kernel, **kwargs), n_jobs=n_jobs)

class SlicedWassersteinKernel(BaseEstimator, TransformerMixin):
    """
    This is a class for computing the sliced Wasserstein kernel matrix from a list of persistence diagrams. The sliced Wasserstein kernel is computed by exponentiating the corresponding sliced Wasserstein distance with a Gaussian kernel. See http://proceedings.mlr.press/v70/carriere17a.html for more details. 
    """
    def __init__(self, num_directions=10, bandwidth=1.0, n_jobs=None):
        """
        Constructor for the SlicedWassersteinKernel class.

        Parameters:
            bandwidth (double): bandwidth of the Gaussian kernel applied to the sliced Wasserstein distance (default 1.).
            num_directions (int): number of lines evenly sampled from [-pi/2,pi/2] in order to approximate and speed up the kernel computation (default 10).
            n_jobs (int): number of jobs to use for the computation. See :func:`pairwise_persistence_diagram_kernels` for details.
        """
        self.bandwidth = bandwidth
        self.num_directions = num_directions
        self.n_jobs = n_jobs

    def fit(self, X, y=None):
        """
        Fit the SlicedWassersteinKernel class on a list of persistence diagrams: an instance of the SlicedWassersteinDistance class is fitted on the diagrams and then stored. 

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
            y (n x 1 array): persistence diagram labels (unused).
        """
        self.diagrams_ = X
        return self

    def transform(self, X):
        """
        Compute all sliced Wasserstein kernel values between the persistence diagrams that were stored after calling the fit() method, and a given list of (possibly different) persistence diagrams.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.

        Returns:
            numpy array of shape (number of diagrams in X) x (number of diagrams in **diagrams**): matrix of pairwise sliced Wasserstein kernel values.
        """
        return pairwise_persistence_diagram_kernels(X, self.diagrams_, kernel="sliced_wasserstein", bandwidth=self.bandwidth, num_directions=self.num_directions, n_jobs=self.n_jobs)

    def __call__(self, diag1, diag2):
        """
        Apply SlicedWassersteinKernel on a single pair of persistence diagrams and outputs the result.

        Parameters:
            diag1 (n x 2 numpy array): first input persistence diagram.
            diag2 (n x 2 numpy array): second input persistence diagram.

        Returns:
            float: sliced Wasserstein kernel value.
        """
        return np.exp(-_sliced_wasserstein_distance(diag1, diag2, num_directions=self.num_directions) / self.bandwidth)

class PersistenceWeightedGaussianKernel(BaseEstimator, TransformerMixin):
    """
    This is a class for computing the persistence weighted Gaussian kernel matrix from a list of persistence diagrams. The persistence weighted Gaussian kernel is computed by convolving the persistence diagram points with weighted Gaussian kernels. See http://proceedings.mlr.press/v48/kusano16.html for more details. 
    """
    def __init__(self, bandwidth=1., weight=lambda x: 1, kernel_approx=None, n_jobs=None):
        """
        Constructor for the PersistenceWeightedGaussianKernel class.
  
        Parameters:
            bandwidth (double): bandwidth of the Gaussian kernel with which persistence diagrams will be convolved (default 1.)
            weight (function): weight function for the persistence diagram points (default constant function, ie lambda x: 1). This function must be defined on 2D points, ie lists or numpy arrays of the form [p_x,p_y].
            kernel_approx (class): kernel approximation class used to speed up computation (default None). Common kernel approximations classes can be found in the scikit-learn library (such as RBFSampler for instance).
            n_jobs (int): number of jobs to use for the computation. See :func:`pairwise_persistence_diagram_kernels` for details.
        """
        self.bandwidth, self.weight = bandwidth, weight
        self.kernel_approx = kernel_approx
        self.n_jobs = n_jobs

    def fit(self, X, y=None):
        """
        Fit the PersistenceWeightedGaussianKernel class on a list of persistence diagrams: persistence diagrams are stored in a numpy array called **diagrams** and the kernel approximation class (if not None) is applied on them. 

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
            y (n x 1 array): persistence diagram labels (unused).
        """
        self.diagrams_ = X
        return self

    def transform(self, X):
        """
        Compute all persistence weighted Gaussian kernel values between the persistence diagrams that were stored after calling the fit() method, and a given list of (possibly different) persistence diagrams.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.

        Returns:
            numpy array of shape (number of diagrams in X) x (number of diagrams in **diagrams**): matrix of pairwise persistence weighted Gaussian kernel values.
        """
        return pairwise_persistence_diagram_kernels(X, self.diagrams_, kernel="persistence_weighted_gaussian", bandwidth=self.bandwidth, weight=self.weight, kernel_approx=self.kernel_approx, n_jobs=self.n_jobs)

    def __call__(self, diag1, diag2):
        """
        Apply PersistenceWeightedGaussianKernel on a single pair of persistence diagrams and outputs the result.

        Parameters:
            diag1 (n x 2 numpy array): first input persistence diagram.
            diag2 (n x 2 numpy array): second input persistence diagram.

        Returns:
            float: persistence weighted Gaussian kernel value.
        """
        return _persistence_weighted_gaussian_kernel(diag1, diag2, weight=self.weight, kernel_approx=self.kernel_approx, bandwidth=self.bandwidth)

class PersistenceScaleSpaceKernel(BaseEstimator, TransformerMixin):
    """
    This is a class for computing the persistence scale space kernel matrix from a list of persistence diagrams. The persistence scale space kernel is computed by adding the symmetric to the diagonal of each point in each persistence diagram, with negative weight, and then convolving the points with a Gaussian kernel. See https://www.cv-foundation.org/openaccess/content_cvpr_2015/papers/Reininghaus_A_Stable_Multi-Scale_2015_CVPR_paper.pdf for more details. 
    """
    def __init__(self, bandwidth=1., kernel_approx=None, n_jobs=None):
        """
        Constructor for the PersistenceScaleSpaceKernel class.
  
        Parameters:
            bandwidth (double): bandwidth of the Gaussian kernel with which persistence diagrams will be convolved (default 1.)
            kernel_approx (class): kernel approximation class used to speed up computation (default None). Common kernel approximations classes can be found in the scikit-learn library (such as RBFSampler for instance).
            n_jobs (int): number of jobs to use for the computation. See :func:`pairwise_persistence_diagram_kernels` for details.
        """
        self.bandwidth, self.kernel_approx = bandwidth, kernel_approx
        self.n_jobs = n_jobs

    def fit(self, X, y=None):
        """
        Fit the PersistenceScaleSpaceKernel class on a list of persistence diagrams: symmetric to the diagonal of all points are computed and an instance of the PersistenceWeightedGaussianKernel class is fitted on the diagrams and then stored. 

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
            y (n x 1 array): persistence diagram labels (unused).
        """
        self.diagrams_ = X
        return self

    def transform(self, X):
        """
        Compute all persistence scale space kernel values between the persistence diagrams that were stored after calling the fit() method, and a given list of (possibly different) persistence diagrams.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.

        Returns:
            numpy array of shape (number of diagrams in X) x (number of diagrams in **diagrams**): matrix of pairwise persistence scale space kernel values.
        """
        return pairwise_persistence_diagram_kernels(X, self.diagrams_, kernel="persistence_scale_space", bandwidth=self.bandwidth, kernel_approx=self.kernel_approx, n_jobs=self.n_jobs)

    def __call__(self, diag1, diag2):
        """
        Apply PersistenceScaleSpaceKernel on a single pair of persistence diagrams and outputs the result.

        Parameters:
            diag1 (n x 2 numpy array): first input persistence diagram.
            diag2 (n x 2 numpy array): second input persistence diagram.

        Returns:
            float: persistence scale space kernel value.
        """
        return _persistence_scale_space_kernel(diag1, diag2, bandwidth=self.bandwidth, kernel_approx=self.kernel_approx)

class PersistenceFisherKernel(BaseEstimator, TransformerMixin):
    """
    This is a class for computing the persistence Fisher kernel matrix from a list of persistence diagrams. The persistence Fisher kernel is computed by exponentiating the corresponding persistence Fisher distance with a Gaussian kernel. See papers.nips.cc/paper/8205-persistence-fisher-kernel-a-riemannian-manifold-kernel-for-persistence-diagrams for more details. 
    """
    def __init__(self, bandwidth_fisher=1., bandwidth=1., kernel_approx=None, n_jobs=None):
        """
        Constructor for the PersistenceFisherKernel class.

        Parameters:
            bandwidth (double): bandwidth of the Gaussian kernel applied to the persistence Fisher distance (default 1.).
            bandwidth_fisher (double): bandwidth of the Gaussian kernel used to turn persistence diagrams into probability distributions by PersistenceFisherDistance class (default 1.).
            kernel_approx (class): kernel approximation class used to speed up computation (default None). Common kernel approximations classes can be found in the scikit-learn library (such as RBFSampler for instance).
            n_jobs (int): number of jobs to use for the computation. See :func:`pairwise_persistence_diagram_kernels` for details.
        """
        self.bandwidth = bandwidth
        self.bandwidth_fisher, self.kernel_approx = bandwidth_fisher, kernel_approx
        self.n_jobs = n_jobs

    def fit(self, X, y=None):
        """
        Fit the PersistenceFisherKernel class on a list of persistence diagrams: an instance of the PersistenceFisherDistance class is fitted on the diagrams and then stored. 

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
            y (n x 1 array): persistence diagram labels (unused).
        """
        self.diagrams_ = X
        return self

    def transform(self, X):
        """
        Compute all persistence Fisher kernel values between the persistence diagrams that were stored after calling the fit() method, and a given list of (possibly different) persistence diagrams.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.

        Returns:
            numpy array of shape (number of diagrams in X) x (number of diagrams in **diagrams**): matrix of pairwise persistence Fisher kernel values.
        """
        return pairwise_persistence_diagram_kernels(X, self.diagrams_, kernel="persistence_fisher", bandwidth=self.bandwidth, bandwidth_fisher=self.bandwidth_fisher, kernel_approx=self.kernel_approx, n_jobs=self.n_jobs)

    def __call__(self, diag1, diag2):
        """
        Apply PersistenceFisherKernel on a single pair of persistence diagrams and outputs the result.

        Parameters:
            diag1 (n x 2 numpy array): first input persistence diagram.
            diag2 (n x 2 numpy array): second input persistence diagram.

        Returns:
            float: persistence Fisher kernel value.
        """
        return np.exp(-_persistence_fisher_distance(diag1, diag2, bandwidth=self.bandwidth_fisher, kernel_approx=self.kernel_approx) / self.bandwidth)

