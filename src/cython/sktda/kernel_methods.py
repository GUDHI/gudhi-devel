"""
@author: Mathieu Carriere
All rights reserved
"""

import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.metrics import pairwise_distances
from .metrics import SlicedWassersteinDistance, PersistenceFisherDistance

#############################################
# Kernel methods ############################
#############################################

class SlicedWassersteinKernel(BaseEstimator, TransformerMixin):
    """
    This is a class for computing the sliced Wasserstein kernel matrix from a list of persistence diagrams. The sliced Wasserstein kernel is computed by exponentiating the corresponding sliced Wasserstein distance with a Gaussian kernel. See http://proceedings.mlr.press/v70/carriere17a.html for more details. 
    """
    def __init__(self, num_directions=10, bandwidth=1.0):
        """
        Constructor for the SlicedWassersteinDistance class.

        Attributes:
            bandwidth (double): bandwidth of the Gaussian kernel applied to the sliced Wasserstein distance (default 1.).
            num_directions (int): number of lines to sample uniformly from [-pi,pi] in order to approximate and speed up the kernel computation (default 10). If -1, the exact kernel is computed.
        """
        self.bandwidth = bandwidth
        self.sw_ = SlicedWassersteinDistance(num_directions=num_directions)

    def fit(self, X, y=None):
        """
        Fit the SlicedWassersteinKernel class on a list of persistence diagrams: an instance of the SlicedWassersteinDistance class is fitted on the diagrams and then stored. 

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
            y (n x 1 array): persistence diagram labels (unused).
        """
        self.sw_.fit(X, y)
        return self

    def transform(self, X):
        """
        Compute all sliced Wasserstein kernel values between the persistence diagrams that were stored after calling the fit() method, and a given list of (possibly different) persistence diagrams.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.

        Returns:
            Xfit (numpy array of shape (number of diagrams in **diagrams**) x (number of diagrams in X)): matrix of pairwise sliced Wasserstein kernel values.
        """
        return np.exp(-self.sw_.transform(X)/self.bandwidth)

class PersistenceWeightedGaussianKernel(BaseEstimator, TransformerMixin):
    """
    This is a class for computing the persistence weighted Gaussian kernel matrix from a list of persistence diagrams. The persistence weighted Gaussian kernel is computed by convolving the persistence diagram points with weighted Gaussian kernels. See http://proceedings.mlr.press/v48/kusano16.html for more details. 
    """
    def __init__(self, bandwidth=1., weight=lambda x: 1, kernel_approx=None):
        """
        Constructor for the PersistenceWeightedGaussianKernel class.
  
        Attributes:
            bandwidth (double): bandwidth of the Gaussian kernel with which persistence diagrams will be convolved (default 1.)
            weight (function): weight function for the persistence diagram points (default constant function, ie lambda x: 1). This function must be defined on 2D points, ie lists or numpy arrays of the form [p_x,p_y].
            kernel_approx (class): kernel approximation class used to speed up computation (default None). Common kernel approximations classes can be found in the scikit-learn library (such as RBFSampler for instance).
        """
        self.bandwidth, self.weight = bandwidth, weight
        self.kernel_approx = kernel_approx

    def fit(self, X, y=None):
        """
        Fit the PersistenceWeightedGaussianKernel class on a list of persistence diagrams: persistence diagrams are stored in a numpy array called **diagrams** and the kernel approximation class (if not None) is applied on them. 

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
            y (n x 1 array): persistence diagram labels (unused).
        """
        self.diagrams_ = list(X)
        self.ws_ = [ np.array([self.weight(self.diagrams_[i][j,:]) for j in range(self.diagrams_[i].shape[0])]) for i in range(len(self.diagrams_)) ]
        if self.kernel_approx is not None:
            self.approx_ = np.concatenate([np.sum(np.multiply(self.ws_[i][:,np.newaxis], self.kernel_approx.transform(self.diagrams_[i])), axis=0)[np.newaxis,:] for i in range(len(self.diagrams_))])
        return self

    def transform(self, X):
        """
        Compute all sliced Wasserstein kernel values between the persistence diagrams that were stored after calling the fit() method, and a given list of (possibly different) persistence diagrams.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.

        Returns:
            Xfit (numpy array of shape (number of diagrams in **diagrams**) x (number of diagrams in X)): matrix of pairwise persistence weighted Gaussian kernel values.
        """
        Xp = list(X)
        Xfit = np.zeros((len(Xp), len(self.diagrams_)))
        if len(self.diagrams_) == len(Xp) and np.all([np.array_equal(self.diagrams_[i], Xp[i]) for i in range(len(Xp))]):
            if self.kernel_approx is not None:
                Xfit = (1./(np.sqrt(2*np.pi)*self.bandwidth)) * np.matmul(self.approx_, self.approx_.T)
            else:
                for i in range(len(self.diagrams_)):
                    for j in range(i+1, len(self.diagrams_)):
                        W = np.matmul(self.ws_[i][:,np.newaxis], self.ws_[j][np.newaxis,:])
                        E = (1./(np.sqrt(2*np.pi)*self.bandwidth)) * np.exp(-np.square(pairwise_distances(self.diagrams_[i], self.diagrams_[j]))/(2*np.square(self.bandwidth)))
                        Xfit[i,j] = np.sum(np.multiply(W, E))
                        Xfit[j,i] = X[i,j]
        else:
            ws = [ np.array([self.weight(Xp[i][j,:]) for j in range(Xp[i].shape[0])]) for i in range(len(Xp)) ]
            if self.kernel_approx is not None:
                approx = np.concatenate([np.sum(np.multiply(ws[i][:,np.newaxis], self.kernel_approx.transform(Xp[i])), axis=0)[np.newaxis,:] for i in range(len(Xp))])
                Xfit = (1./(np.sqrt(2*np.pi)*self.bandwidth)) * np.matmul(approx, self.approx_.T)
            else:
                for i in range(len(Xp)):
                    for j in range(len(self.diagrams_)):
                        W = np.matmul(ws[i][:,np.newaxis], self.ws_[j][np.newaxis,:])
                        E = (1./(np.sqrt(2*np.pi)*self.bandwidth)) * np.exp(-np.square(pairwise_distances(Xp[i], self.diagrams_[j]))/(2*np.square(self.bandwidth)))
                        Xfit[i,j] = np.sum(np.multiply(W, E))
        
        return Xfit

class PersistenceScaleSpaceKernel(BaseEstimator, TransformerMixin):
    """
    This is a class for computing the persistence scale space kernel matrix from a list of persistence diagrams. The persistence scale space kernel is computed by adding the symmetric to the diagonal of each point in each persistence diagram, and then convolving the points with a Gaussian kernel. See https://www.cv-foundation.org/openaccess/content_cvpr_2015/papers/Reininghaus_A_Stable_Multi-Scale_2015_CVPR_paper.pdf for more details. 
    """
    def __init__(self, bandwidth=1., kernel_approx=None):
        """
        Constructor for the PersistenceScaleSpaceKernel class.
  
        Attributes:
            bandwidth (double): bandwidth of the Gaussian kernel with which persistence diagrams will be convolved (default 1.)
            kernel_approx (class): kernel approximation class used to speed up computation (default None). Common kernel approximations classes can be found in the scikit-learn library (such as RBFSampler for instance).
        """
        self.pwg_ = PersistenceWeightedGaussianKernel(bandwidth=bandwidth, weight=lambda x: 1 if x[1] >= x[0] else -1, kernel_approx=kernel_approx)

    def fit(self, X, y=None):
        """
        Fit the PersistenceScaleSpaceKernel class on a list of persistence diagrams: symmetric to the diagonal of all points are computed and an instance of the PersistenceWeightedGaussianKernel class is fitted on the diagrams and then stored. 

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
            y (n x 1 array): persistence diagram labels (unused).
        """
        self.diagrams_ = list(X)
        for i in range(len(self.diagrams_)):
            op_D = np.matmul(self.diagrams_[i], np.array([[0.,1.], [1.,0.]]))
            self.diagrams_[i] = np.concatenate([self.diagrams_[i], op_D], axis=0)
        self.pwg_.fit(X)
        return self

    def transform(self, X):
        """
        Compute all persistence scale space kernel values between the persistence diagrams that were stored after calling the fit() method, and a given list of (possibly different) persistence diagrams.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.

        Returns:
            Xfit (numpy array of shape (number of diagrams in **diagrams**) x (number of diagrams in X)): matrix of pairwise persistence scale space kernel values.
        """
        Xp = list(X)
        for i in range(len(Xp)):
            op_X = np.matmul(Xp[i], np.array([[0.,1.], [1.,0.]]))
            Xp[i] = np.concatenate([Xp[i], op_X], axis=0)
        return self.pwg_.transform(Xp)

class PersistenceFisherKernel(BaseEstimator, TransformerMixin):
    """
    This is a class for computing the persistence Fisher kernel matrix from a list of persistence diagrams. The persistence Fisher kernel is computed by exponentiating the corresponding persistence Fisher distance with a Gaussian kernel. See papers.nips.cc/paper/8205-persistence-fisher-kernel-a-riemannian-manifold-kernel-for-persistence-diagrams for more details. 
    """
    def __init__(self, bandwidth_fisher=1., bandwidth=1., kernel_approx=None):
        """
        Constructor for the PersistenceFisherKernel class.

        Attributes:
            bandwidth (double): bandwidth of the Gaussian kernel applied to the persistence Fisher distance (default 1.).
            bandwidth_fisher (double): bandwidth of the Gaussian kernel used to turn persistence diagrams into probability distributions by PersistenceFisherDistance class (default 1.).
            kernel_approx (class): kernel approximation class used to speed up computation (default None). Common kernel approximations classes can be found in the scikit-learn library (such as RBFSampler for instance).
        """
        self.bandwidth = bandwidth
        self.pf_ = PersistenceFisherDistance(bandwidth=bandwidth_fisher, kernel_approx=kernel_approx)

    def fit(self, X, y=None):
        """
        Fit the PersistenceFisherKernel class on a list of persistence diagrams: an instance of the PersistenceFisherDistance class is fitted on the diagrams and then stored. 

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
            y (n x 1 array): persistence diagram labels (unused).
        """
        self.pf_.fit(X, y)
        return self

    def transform(self, X):
        """
        Compute all persistence Fisher kernel values between the persistence diagrams that were stored after calling the fit() method, and a given list of (possibly different) persistence diagrams.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.

        Returns:
            Xfit (numpy array of shape (number of diagrams in **diagrams**) x (number of diagrams in X)): matrix of pairwise persistence Fisher kernel values.
        """
        return np.exp(-self.pf_.transform(X)/self.bandwidth)

