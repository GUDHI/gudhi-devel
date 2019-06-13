"""
@author: Mathieu Carriere
All rights reserved
"""

import numpy as np
from sklearn.base          import BaseEstimator, TransformerMixin
from sklearn.preprocessing import MinMaxScaler
from sklearn.neighbors     import DistanceMetric

from preprocessing import DiagramPreprocessor

#############################################
# Finite Vectorization methods ##############
#############################################

class PersistenceImage(BaseEstimator, TransformerMixin):
    """
    This is a class for computing persistence images from a list of persistence diagrams. A persistence image is a 2D function computed from a persistence diagram by convolving the diagram points with a weighted Gaussian kernel. The plane is then discretized into an image with pixels, which is flattened and returned as a vector. See http://jmlr.org/papers/v18/16-337.html for more details.
    """
    def __init__(self, bandwidth=1., weight=lambda x: 1, resolution=[20,20], im_range=[np.nan, np.nan, np.nan, np.nan]):
        """
        Constructor for the PersistenceImage class.

        Attributes:
            bandwidth (double): bandwidth of the Gaussian kernel (default 1.).
            weight (function): weight function for the persistence diagram points (default constant function, ie lambda x: 1). This function must be defined on 2D points, ie lists or numpy arrays of the form [p_x,p_y].
            resolution ([int,int]): size (in pixels) of the persistence image (default [20,20]).
            im_range ([double,double,double,double]): minimum and maximum of each axis of the persistence image, of the form [x_min, x_max, y_min, y_max] (default [numpy.nan, numpy.nan, numpy.nan, numpyp.nan]). If one of the values is numpy.nan, it can be computed from the persistence diagrams with the fit() method.
        """
        self.bandwidth, self.weight = bandwidth, weight
        self.resolution, self.im_range = resolution, im_range

    def fit(self, X, y=None):
        """
        Fit the PersistenceImage class on a list of persistence diagrams: if any of the values in **im_range** is numpy.nan, replace it with the corresponding value computed on the given list of persistence diagrams.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
            y (n x 1 array): persistence diagram labels (unused).
        """
        if np.isnan(np.array(self.im_range)).any():
            pre = DiagramPreprocessor(use=True, scalers=[([0,1], MinMaxScaler())]).fit(X,y)
            [mx,my],[Mx,My] = pre.scalers[0][1].data_min_, pre.scalers[0][1].data_max_
            self.im_range = np.where(np.isnan(np.array(self.im_range)), np.array([mx, Mx, my, My]), np.array(self.im_range))
        return self

    def transform(self, X):
        """
        Compute the persistence image for each persistence diagram individually and concatenate the results.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
    
        Returns:
            Xfit (numpy array with shape (number of diagrams) x (number of pixels = **resolution[0]** x **resolution[1]**)): output persistence images.
        """
        num_diag, Xfit = len(X), []
        for i in range(num_diag):

            diagram, num_pts_in_diag = X[i], X[i].shape[0]

            w = np.ones(num_pts_in_diag)
            for j in range(num_pts_in_diag):
                w[j] = self.weight(diagram[j,:])

            x_values, y_values = np.linspace(self.im_range[0], self.im_range[1], self.resolution[0]), np.linspace(self.im_range[2], self.im_range[3], self.resolution[1])
            Xs, Ys = np.tile((diagram[:,0][:,np.newaxis,np.newaxis]-x_values[np.newaxis,np.newaxis,:]),[1,self.resolution[1],1]), np.tile(diagram[:,1][:,np.newaxis,np.newaxis]-y_values[np.newaxis,:,np.newaxis],[1,1,self.resolution[0]])
            image = np.tensordot(w, np.exp((-np.square(Xs)-np.square(Ys))/(2*np.square(self.bandwidth)))/(self.bandwidth*np.sqrt(2*np.pi)), 1)

            Xfit.append(image.flatten()[np.newaxis,:])
            Xfit = np.concatenate(Xfit,0)

        return Xfit

class Landscape(BaseEstimator, TransformerMixin):
    """
    This is a class for computing persistence landscapes from a list of persistence diagrams. A persistence landscape is a collection of 1D piecewise-linear functions computed from the rank function associated to the persistence diagram. These piecewise-linear functions are then sampled uniformly on a given range and the corresponding vectors of samples are concatenated and returned. See http://jmlr.org/papers/v16/bubenik15a.html for more details.
    """
    def __init__(self, num_landscapes=5, resolution=100, ls_range=[np.nan, np.nan]):
        """
        Constructor for the Landscape class.

        Attributes:
            num_landscapes (int): number of piecewise-linear functions to output (default 5).
            resolution (int): number of sample for all piecewise-linear functions (default 100).
            ls_range ([double, double]): minimum and maximum of all piecewise-linear function domains, of the form [x_min, x_max] (default [numpy.nan, numpy.nan]). It is the interval on which samples will be drawn uniformly. If one of the values is numpy.nan, it can be computed from the persistence diagrams with the fit() method.
        """
        self.num_landscapes, self.resolution, self.ls_range = num_landscapes, resolution, ls_range

    def fit(self, X, y=None):
        """
        Fit the Landscape class on a list of persistence diagrams: if any of the values in **ls_range** is numpy.nan, replace it with the corresponding value computed on the given list of persistence diagrams.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
            y (n x 1 array): persistence diagram labels (unused).
        """
        if np.isnan(np.array(self.ls_range)).any():
            pre = DiagramPreprocessor(use=True, scalers=[([0,1], MinMaxScaler())]).fit(X,y)
            [mx,my],[Mx,My] = pre.scalers[0][1].data_min_, pre.scalers[0][1].data_max_
            self.ls_range = np.where(np.isnan(np.array(self.ls_range)), np.array([mx, My]), np.array(self.ls_range))
        return self

    def transform(self, X):
        """
        Compute the persistence landscape for each persistence diagram individually and concatenate the results.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
    
        Returns:
            Xfit (numpy array with shape (number of diagrams) x (number of samples = **num_landscapes** x **resolution**)): output persistence landscapes.
        """
        num_diag, Xfit = len(X), []
        x_values = np.linspace(self.ls_range[0], self.ls_range[1], self.resolution)
        step_x = x_values[1] - x_values[0]

        for i in range(num_diag):

            diagram, num_pts_in_diag = X[i], X[i].shape[0]

            ls = np.zeros([self.num_landscapes, self.resolution])

            events = []
            for j in range(self.resolution):
                events.append([])

            for j in range(num_pts_in_diag):
                [px,py] = diagram[j,:]
                min_idx = np.minimum(np.maximum(np.ceil((px          - self.ls_range[0]) / step_x).astype(int), 0), self.resolution)
                mid_idx = np.minimum(np.maximum(np.ceil((0.5*(py+px) - self.ls_range[0]) / step_x).astype(int), 0), self.resolution)
                max_idx = np.minimum(np.maximum(np.ceil((py          - self.ls_range[0]) / step_x).astype(int), 0), self.resolution)

                if min_idx < self.resolution and max_idx > 0:

                    landscape_value = self.ls_range[0] + min_idx * step_x - px
                    for k in range(min_idx, mid_idx):
                        events[k].append(landscape_value)
                        landscape_value += step_x

                    landscape_value = py - self.ls_range[0] - mid_idx * step_x
                    for k in range(mid_idx, max_idx):
                        events[k].append(landscape_value)
                        landscape_value -= step_x

            for j in range(self.resolution):
                events[j].sort(reverse=True)
                for k in range( min(self.num_landscapes, len(events[j])) ):
                    ls[k,j] = events[j][k]

            Xfit.append(np.sqrt(2)*np.reshape(ls,[1,-1]))
            Xfit = np.concatenate(Xfit,0)

        return Xfit

class Silhouette(BaseEstimator, TransformerMixin):
    """
    This is a class for computing persistence silhouettes from a list of persistence diagrams. A persistence silhouette is computed by taking a weighted average of the collection of 1D piecewise-linear functions given by the persistence landscapes, and then by uniformly sampling this average on a given range. Finally, the corresponding vector of samples is returned. See https://arxiv.org/abs/1312.0308 for more details.
    """
    def __init__(self, weight=lambda x: 1, resolution=100, sh_range=[np.nan, np.nan]):
        """
        Constructor for the Silhouette class.

        Attributes:
            weight (function): weight function for the persistence diagram points (default constant function, ie lambda x: 1). This function must be defined on 2D points, ie on lists or numpy arrays of the form [p_x,p_y].
            resolution (int): number of samples for the weighted average (default 100).
            sh_range ([double, double]): minimum and maximum for the weighted average domain, of the form [x_min, x_max] (default [numpy.nan, numpy.nan]). It is the interval on which samples will be drawn uniformly. If one of the values is numpy.nan, it can be computed from the persistence diagrams with the fit() method.
        """
        self.weight, self.resolution, self.sh_range = weight, resolution, sh_range

    def fit(self, X, y=None):
        """
        Fit the Silhouette class on a list of persistence diagrams: if any of the values in **sh_range** is numpy.nan, replace it with the corresponding value computed on the given list of persistence diagrams.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
            y (n x 1 array): persistence diagram labels (unused).
        """
        if np.isnan(np.array(self.sh_range)).any():
            pre = DiagramPreprocessor(use=True, scalers=[([0,1], MinMaxScaler())]).fit(X,y)
            [mx,my],[Mx,My] = pre.scalers[0][1].data_min_, pre.scalers[0][1].data_max_
            self.sh_range = np.where(np.isnan(np.array(self.sh_range)), np.array([mx, My]), np.array(self.sh_range))
        return self

    def transform(self, X):
        """
        Compute the persistence silhouette for each persistence diagram individually and concatenate the results.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
    
        Returns:
            Xfit (numpy array with shape (number of diagrams) x (**resolution**): output persistence silhouettes.
        """
        num_diag, Xfit = len(X), []
        x_values = np.linspace(self.sh_range[0], self.sh_range[1], self.resolution)
        step_x = x_values[1] - x_values[0]

        for i in range(num_diag):

            diagram, num_pts_in_diag = X[i], X[i].shape[0]

            sh, weights = np.zeros(self.resolution), np.zeros(num_pts_in_diag)
            for j in range(num_pts_in_diag):
                weights[j] = self.weight(diagram[j,:])
            total_weight = np.sum(weights)

            for j in range(num_pts_in_diag):

                [px,py] = diagram[j,:]
                weight  = weights[j] / total_weight
                min_idx = np.minimum(np.maximum(np.ceil((px          - self.sh_range[0]) / step_x).astype(int), 0), self.resolution)
                mid_idx = np.minimum(np.maximum(np.ceil((0.5*(py+px) - self.sh_range[0]) / step_x).astype(int), 0), self.resolution)
                max_idx = np.minimum(np.maximum(np.ceil((py          - self.sh_range[0]) / step_x).astype(int), 0), self.resolution)

                if min_idx < self.resolution and max_idx > 0:

                    silhouette_value = self.sh_range[0] + min_idx * step_x - px
                    for k in range(min_idx, mid_idx):
                        sh[k] += weight * silhouette_value
                        silhouette_value += step_x

                    silhouette_value = py - self.sh_range[0] - mid_idx * step_x
                    for k in range(mid_idx, max_idx):
                        sh[k] += weight * silhouette_value
                        silhouette_value -= step_x

            Xfit.append(np.reshape(np.sqrt(2) * sh, [1,-1]))
            Xfit = np.concatenate(Xfit, 0)

        return Xfit 

class BettiCurve(BaseEstimator, TransformerMixin):
    """
    This is a class for computing Betti curves from a list of persistence diagrams. A Betti curve is a 1D piecewise-constant function obtained from the rank function. It is sampled uniformly on a given range and the vector of samples is returned. See https://www.researchgate.net/publication/316604237_Time_Series_Classification_via_Topological_Data_Analysis for more details.
    """
    def __init__(self, resolution=100, bc_range=[np.nan, np.nan]):
        """
        Constructor for the BettiCurve class.

        Attributes:
            resolution (int): number of sample for the piecewise-constant function (default 100).
            bc_range ([double, double]): minimum and maximum of the piecewise-constant function domain, of the form [x_min, x_max] (default [numpy.nan, numpy.nan]). It is the interval on which samples will be drawn uniformly. If one of the values is numpy.nan, it can be computed from the persistence diagrams with the fit() method.
        """
        self.resolution, self.bc_range = resolution, bc_range

    def fit(self, X, y=None):
        """
        Fit the BettiCurve class on a list of persistence diagrams: if any of the values in **bc_range** is numpy.nan, replace it with the corresponding value computed on the given list of persistence diagrams.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
            y (n x 1 array): persistence diagram labels (unused).
        """
        if np.isnan(np.array(self.bc_range)).any():
            pre = DiagramPreprocessor(use=True, scalers=[([0,1], MinMaxScaler())]).fit(X,y)
            [mx,my],[Mx,My] = pre.scalers[0][1].data_min_, pre.scalers[0][1].data_max_
            self.bc_range = np.where(np.isnan(np.array(self.bc_range)), np.array([mx, My]), np.array(self.bc_range))
        return self

    def transform(self, X):
        """
        Compute the Betti curve for each persistence diagram individually and concatenate the results.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
    
        Returns:
            Xfit (numpy array with shape (number of diagrams) x (**resolution**): output Betti curves.
        """
        num_diag, Xfit = len(X), []
        x_values = np.linspace(self.bc_range[0], self.bc_range[1], self.resolution)
        step_x = x_values[1] - x_values[0]

        for i in range(num_diag):

            diagram, num_pts_in_diag = X[i], X[i].shape[0]

            bc =  np.zeros(self.resolution)
            for j in range(num_pts_in_diag):
                [px,py] = diagram[j,:]
                min_idx = np.minimum(np.maximum(np.ceil((px - self.bc_range[0]) / step_x).astype(int), 0), self.resolution)
                max_idx = np.minimum(np.maximum(np.ceil((py - self.bc_range[0]) / step_x).astype(int), 0), self.resolution)
                for k in range(min_idx, max_idx):
                    bc[k] += 1

            Xfit.append(np.reshape(bc,[1,-1]))
            Xfit = np.concatenate(Xfit, 0)

        return Xfit

class TopologicalVector(BaseEstimator, TransformerMixin):
    """
    This is a class for computing topological vectors from a list of persistence diagrams. The topological vector associated to a persistence diagram is the sorted vector of a slight modification of the pairwise distances between the persistence diagram points. See https://diglib.eg.org/handle/10.1111/cgf12692 for more details.
    """
    def __init__(self, threshold=10):
        """
        Constructor for the TopologicalVector class.

        Attributes:
            threshold (int): number of distances to keep (default 10). This is the dimension of the topological vector. If -1, this threshold is computed from the list of persistence diagrams by considering the one with the largest number of points and using the dimension of its corresponding topological vector as threshold. 
        """
        self.threshold = threshold

    def fit(self, X, y=None):
        """
        Fit the TopologicalVector class on a list of persistence diagrams (this function actually does nothing but is useful when TopologicalVector is included in a scikit-learn Pipeline).

        Parameters:
            X (list of n x 2 or n x 1 numpy arrays): input persistence diagrams.
            y (n x 1 array): persistence diagram labels (unused).
        """
        return self

    def transform(self, X):
        """
        Compute the topological vector for each persistence diagram individually and concatenate the results.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
    
        Returns:
            Xfit (numpy array with shape (number of diagrams) x (**threshold**): output topological vectors.
        """
        if self.threshold == -1:
            thresh = np.array([X[i].shape[0] for i in range(len(X))]).max()
        else:
            thresh = self.threshold

        num_diag = len(X)
        Xfit = np.zeros([num_diag, thresh])

        for i in range(num_diag):

            diagram, num_pts_in_diag = X[i], X[i].shape[0]
            pers = 0.5 * np.matmul(diagram, np.array([[-1.0],[1.0]]))
            min_pers = np.minimum(pers,np.transpose(pers))
            distances = DistanceMetric.get_metric("chebyshev").pairwise(diagram)
            vect = np.flip(np.sort(np.triu(np.minimum(distances, min_pers)), axis=None), 0)
            dim = min(len(vect), thresh)
            Xfit[i, :dim] = vect[:dim]

        return Xfit

class ComplexPolynomial(BaseEstimator, TransformerMixin):
    """
    This is a class for computing complex polynomials from a list of persistence diagrams. The persistence diagram points are seen as the roots of some complex polynomial, whose coefficients are returned in a complex vector. See https://link.springer.com/chapter/10.1007%2F978-3-319-23231-7_27 for more details.
    """
    def __init__(self, F="R", threshold=10):
        """
        Constructor for the ComplexPolynomial class.

        Attributes:
           F (char): either "R", "S" or "T" (default "R"). Type of complex polynomial that is going to be computed (explained in https://link.springer.com/chapter/10.1007%2F978-3-319-23231-7_27).
           threshold (int): number of coefficients (default 10). This is the dimension of the complex vector of coefficients. If -1, this threshold is computed from the list of persistence diagrams by considering the one with the largest number of points and using the dimension of its corresponding complex vector of coefficients as threshold. 
        """
        self.threshold, self.F = threshold, F

    def fit(self, X, y=None):
        """
        Fit the ComplexPolynomial class on a list of persistence diagrams (this function actually does nothing but is useful when ComplexPolynomial is included in a scikit-learn Pipeline).

        Parameters:
            X (list of n x 2 or n x 1 numpy arrays): input persistence diagrams.
            y (n x 1 array): persistence diagram labels (unused).
        """
        return self

    def transform(self, X):
        """
        Compute the complex vector of coefficients for each persistence diagram individually and concatenate the results.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
    
        Returns:
            Xfit (numpy array with shape (number of diagrams) x (**threshold**): output complex vectors of coefficients.
        """
        if self.threshold == -1:
            thresh = np.array([X[i].shape[0] for i in range(len(X))]).max()
        else:
            thresh = self.threshold

        Xfit = np.zeros([len(X), thresh]) + 1j * np.zeros([len(X), thresh])
        for d in range(len(X)):
            D, N = X[d], X[d].shape[0]
            if self.F == "R":
                roots = D[:,0] + 1j * D[:,1]
            elif self.F == "S":
                alpha = np.linalg.norm(D, axis=1)
                alpha = np.where(alpha==0, np.ones(N), alpha)
                roots = np.multiply( np.multiply(  (D[:,0]+1j*D[:,1]), (D[:,1]-D[:,0])  ), 1./(np.sqrt(2)*alpha) )
            elif self.F == "T":
                alpha = np.linalg.norm(D, axis=1)
                roots = np.multiply(  (D[:,1]-D[:,0])/2, np.cos(alpha) - np.sin(alpha) + 1j * (np.cos(alpha) + np.sin(alpha))  )
            coeff = [0] * (N+1)
            coeff[N] = 1
            for i in range(1, N+1): 
                for j in range(N-i-1, N): 
                    coeff[j] += ((-1) * roots[i-1] * coeff[j+1])   
            coeff = np.array(coeff[::-1])[1:]
            Xfit[d, :min(thresh, coeff.shape[0])] = coeff[:min(thresh, coeff.shape[0])]
        return Xfit
