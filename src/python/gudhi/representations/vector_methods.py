# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Mathieu Carrière, Martin Royer
#
# Copyright (C) 2018-2020 Inria
#
# Modification(s):
#   - 2020/06 Martin: ATOL integration

import numpy as np
from sklearn.base          import BaseEstimator, TransformerMixin
from sklearn.preprocessing import MinMaxScaler, MaxAbsScaler
from sklearn.neighbors     import DistanceMetric
from sklearn.metrics       import pairwise

from .preprocessing import DiagramScaler, BirthPersistenceTransform

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

        Parameters:
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
            new_X = BirthPersistenceTransform().fit_transform(X)
            pre = DiagramScaler(use=True, scalers=[([0], MinMaxScaler()), ([1], MinMaxScaler())]).fit(new_X,y)
            [mx,my],[Mx,My] = [pre.scalers[0][1].data_min_[0], pre.scalers[1][1].data_min_[0]], [pre.scalers[0][1].data_max_[0], pre.scalers[1][1].data_max_[0]]
            self.im_range = np.where(np.isnan(np.array(self.im_range)), np.array([mx, Mx, my, My]), np.array(self.im_range))
        return self

    def transform(self, X):
        """
        Compute the persistence image for each persistence diagram individually and store the results in a single numpy array.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
    
        Returns:
            numpy array with shape (number of diagrams) x (number of pixels = **resolution[0]** x **resolution[1]**): output persistence images.
        """
        num_diag, Xfit = len(X), []
        new_X = BirthPersistenceTransform().fit_transform(X)

        for i in range(num_diag):

            diagram, num_pts_in_diag = new_X[i], X[i].shape[0]

            w = np.empty(num_pts_in_diag)
            for j in range(num_pts_in_diag):
                w[j] = self.weight(diagram[j,:])

            x_values, y_values = np.linspace(self.im_range[0], self.im_range[1], self.resolution[0]), np.linspace(self.im_range[2], self.im_range[3], self.resolution[1])
            Xs, Ys = np.tile((diagram[:,0][:,np.newaxis,np.newaxis]-x_values[np.newaxis,np.newaxis,:]),[1,self.resolution[1],1]), np.tile(diagram[:,1][:,np.newaxis,np.newaxis]-y_values[np.newaxis,:,np.newaxis],[1,1,self.resolution[0]])
            image = np.tensordot(w, np.exp((-np.square(Xs)-np.square(Ys))/(2*np.square(self.bandwidth)))/(np.square(self.bandwidth)*2*np.pi), 1)

            Xfit.append(image.flatten()[np.newaxis,:])

        Xfit = np.concatenate(Xfit,0)

        return Xfit

    def __call__(self, diag):
        """
        Apply PersistenceImage on a single persistence diagram and outputs the result.

        Parameters:
            diag (n x 2 numpy array): input persistence diagram.

        Returns:
            numpy array with shape (number of pixels = **resolution[0]** x **resolution[1]**):: output persistence image.
        """
        return self.fit_transform([diag])[0,:]

class Landscape(BaseEstimator, TransformerMixin):
    """
    This is a class for computing persistence landscapes from a list of persistence diagrams. A persistence landscape is a collection of 1D piecewise-linear functions computed from the rank function associated to the persistence diagram. These piecewise-linear functions are then sampled evenly on a given range and the corresponding vectors of samples are concatenated and returned. See http://jmlr.org/papers/v16/bubenik15a.html for more details.
    """
    def __init__(self, num_landscapes=5, resolution=100, sample_range=[np.nan, np.nan]):
        """
        Constructor for the Landscape class.

        Parameters:
            num_landscapes (int): number of piecewise-linear functions to output (default 5).
            resolution (int): number of sample for all piecewise-linear functions (default 100).
            sample_range ([double, double]): minimum and maximum of all piecewise-linear function domains, of the form [x_min, x_max] (default [numpy.nan, numpy.nan]). It is the interval on which samples will be drawn evenly. If one of the values is numpy.nan, it can be computed from the persistence diagrams with the fit() method.
        """
        self.num_landscapes, self.resolution, self.sample_range = num_landscapes, resolution, sample_range
        self.nan_in_range = np.isnan(np.array(self.sample_range))
        self.new_resolution = self.resolution + self.nan_in_range.sum()

    def fit(self, X, y=None):
        """
        Fit the Landscape class on a list of persistence diagrams: if any of the values in **sample_range** is numpy.nan, replace it with the corresponding value computed on the given list of persistence diagrams.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
            y (n x 1 array): persistence diagram labels (unused).
        """
        if self.nan_in_range.any():
            pre = DiagramScaler(use=True, scalers=[([0], MinMaxScaler()), ([1], MinMaxScaler())]).fit(X,y)
            [mx,my],[Mx,My] = [pre.scalers[0][1].data_min_[0], pre.scalers[1][1].data_min_[0]], [pre.scalers[0][1].data_max_[0], pre.scalers[1][1].data_max_[0]]
            self.sample_range = np.where(self.nan_in_range, np.array([mx, My]), np.array(self.sample_range))
        return self

    def transform(self, X):
        """
        Compute the persistence landscape for each persistence diagram individually and concatenate the results.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
    
        Returns:
            numpy array with shape (number of diagrams) x (number of samples = **num_landscapes** x **resolution**): output persistence landscapes.
        """
        num_diag, Xfit = len(X), []
        x_values = np.linspace(self.sample_range[0], self.sample_range[1], self.new_resolution)
        step_x = x_values[1] - x_values[0]

        for i in range(num_diag):

            diagram, num_pts_in_diag = X[i], X[i].shape[0]

            ls = np.zeros([self.num_landscapes, self.new_resolution])

            events = []
            for j in range(self.new_resolution):
                events.append([])

            for j in range(num_pts_in_diag):
                [px,py] = diagram[j,:2]
                min_idx = np.clip(np.ceil((px          - self.sample_range[0]) / step_x).astype(int), 0, self.new_resolution)
                mid_idx = np.clip(np.ceil((0.5*(py+px) - self.sample_range[0]) / step_x).astype(int), 0, self.new_resolution)
                max_idx = np.clip(np.ceil((py          - self.sample_range[0]) / step_x).astype(int), 0, self.new_resolution)

                if min_idx < self.new_resolution and max_idx > 0:

                    landscape_value = self.sample_range[0] + min_idx * step_x - px
                    for k in range(min_idx, mid_idx):
                        events[k].append(landscape_value)
                        landscape_value += step_x

                    landscape_value = py - self.sample_range[0] - mid_idx * step_x
                    for k in range(mid_idx, max_idx):
                        events[k].append(landscape_value)
                        landscape_value -= step_x

            for j in range(self.new_resolution):
                events[j].sort(reverse=True)
                for k in range( min(self.num_landscapes, len(events[j])) ):
                    ls[k,j] = events[j][k]

            if self.nan_in_range[0]:
                ls = ls[:,1:]
            if self.nan_in_range[1]:
                ls = ls[:,:-1]
            ls = np.sqrt(2)*np.reshape(ls,[1,-1])
            Xfit.append(ls)

        Xfit = np.concatenate(Xfit,0)

        return Xfit

    def __call__(self, diag):
        """
        Apply Landscape on a single persistence diagram and outputs the result.

        Parameters:
            diag (n x 2 numpy array): input persistence diagram.

        Returns:
            numpy array with shape (number of samples = **num_landscapes** x **resolution**): output persistence landscape.
        """
        return self.fit_transform([diag])[0,:]

class Silhouette(BaseEstimator, TransformerMixin):
    """
    This is a class for computing persistence silhouettes from a list of persistence diagrams. A persistence silhouette is computed by taking a weighted average of the collection of 1D piecewise-linear functions given by the persistence landscapes, and then by evenly sampling this average on a given range. Finally, the corresponding vector of samples is returned. See https://arxiv.org/abs/1312.0308 for more details.
    """
    def __init__(self, weight=lambda x: 1, resolution=100, sample_range=[np.nan, np.nan]):
        """
        Constructor for the Silhouette class.

        Parameters:
            weight (function): weight function for the persistence diagram points (default constant function, ie lambda x: 1). This function must be defined on 2D points, ie on lists or numpy arrays of the form [p_x,p_y].
            resolution (int): number of samples for the weighted average (default 100).
            sample_range ([double, double]): minimum and maximum for the weighted average domain, of the form [x_min, x_max] (default [numpy.nan, numpy.nan]). It is the interval on which samples will be drawn evenly. If one of the values is numpy.nan, it can be computed from the persistence diagrams with the fit() method.
        """
        self.weight, self.resolution, self.sample_range = weight, resolution, sample_range

    def fit(self, X, y=None):
        """
        Fit the Silhouette class on a list of persistence diagrams: if any of the values in **sample_range** is numpy.nan, replace it with the corresponding value computed on the given list of persistence diagrams.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
            y (n x 1 array): persistence diagram labels (unused).
        """
        if np.isnan(np.array(self.sample_range)).any():
            pre = DiagramScaler(use=True, scalers=[([0], MinMaxScaler()), ([1], MinMaxScaler())]).fit(X,y)
            [mx,my],[Mx,My] = [pre.scalers[0][1].data_min_[0], pre.scalers[1][1].data_min_[0]], [pre.scalers[0][1].data_max_[0], pre.scalers[1][1].data_max_[0]]
            self.sample_range = np.where(np.isnan(np.array(self.sample_range)), np.array([mx, My]), np.array(self.sample_range))
        return self

    def transform(self, X):
        """
        Compute the persistence silhouette for each persistence diagram individually and concatenate the results.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
    
        Returns:
            numpy array with shape (number of diagrams) x (**resolution**): output persistence silhouettes.
        """
        num_diag, Xfit = len(X), []
        x_values = np.linspace(self.sample_range[0], self.sample_range[1], self.resolution)
        step_x = x_values[1] - x_values[0]

        for i in range(num_diag):

            diagram, num_pts_in_diag = X[i], X[i].shape[0]

            sh, weights = np.zeros(self.resolution), np.zeros(num_pts_in_diag)
            for j in range(num_pts_in_diag):
                weights[j] = self.weight(diagram[j,:])
            total_weight = np.sum(weights)

            for j in range(num_pts_in_diag):

                [px,py] = diagram[j,:2]
                weight  = weights[j] / total_weight
                min_idx = np.clip(np.ceil((px          - self.sample_range[0]) / step_x).astype(int), 0, self.resolution)
                mid_idx = np.clip(np.ceil((0.5*(py+px) - self.sample_range[0]) / step_x).astype(int), 0, self.resolution)
                max_idx = np.clip(np.ceil((py          - self.sample_range[0]) / step_x).astype(int), 0, self.resolution)

                if min_idx < self.resolution and max_idx > 0:

                    silhouette_value = self.sample_range[0] + min_idx * step_x - px
                    for k in range(min_idx, mid_idx):
                        sh[k] += weight * silhouette_value
                        silhouette_value += step_x

                    silhouette_value = py - self.sample_range[0] - mid_idx * step_x
                    for k in range(mid_idx, max_idx):
                        sh[k] += weight * silhouette_value
                        silhouette_value -= step_x

            Xfit.append(np.reshape(np.sqrt(2) * sh, [1,-1]))

        Xfit = np.concatenate(Xfit, 0)

        return Xfit 

    def __call__(self, diag):
        """
        Apply Silhouette on a single persistence diagram and outputs the result.

        Parameters:
            diag (n x 2 numpy array): input persistence diagram.

        Returns:
            numpy array with shape (**resolution**): output persistence silhouette.
        """
        return self.fit_transform([diag])[0,:]

class BettiCurve(BaseEstimator, TransformerMixin):
    """
    This is a class for computing Betti curves from a list of persistence diagrams. A Betti curve is a 1D piecewise-constant function obtained from the rank function. It is sampled evenly on a given range and the vector of samples is returned. See https://www.researchgate.net/publication/316604237_Time_Series_Classification_via_Topological_Data_Analysis for more details.
    """
    def __init__(self, resolution=100, sample_range=[np.nan, np.nan]):
        """
        Constructor for the BettiCurve class.

        Parameters:
            resolution (int): number of sample for the piecewise-constant function (default 100).
            sample_range ([double, double]): minimum and maximum of the piecewise-constant function domain, of the form [x_min, x_max] (default [numpy.nan, numpy.nan]). It is the interval on which samples will be drawn evenly. If one of the values is numpy.nan, it can be computed from the persistence diagrams with the fit() method.
        """
        self.resolution, self.sample_range = resolution, sample_range

    def fit(self, X, y=None):
        """
        Fit the BettiCurve class on a list of persistence diagrams: if any of the values in **sample_range** is numpy.nan, replace it with the corresponding value computed on the given list of persistence diagrams.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
            y (n x 1 array): persistence diagram labels (unused).
        """
        if np.isnan(np.array(self.sample_range)).any():
            pre = DiagramScaler(use=True, scalers=[([0], MinMaxScaler()), ([1], MinMaxScaler())]).fit(X,y)
            [mx,my],[Mx,My] = [pre.scalers[0][1].data_min_[0], pre.scalers[1][1].data_min_[0]], [pre.scalers[0][1].data_max_[0], pre.scalers[1][1].data_max_[0]]
            self.sample_range = np.where(np.isnan(np.array(self.sample_range)), np.array([mx, My]), np.array(self.sample_range))
        return self

    def transform(self, X):
        """
        Compute the Betti curve for each persistence diagram individually and concatenate the results.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
    
        Returns:
            numpy array with shape (number of diagrams) x (**resolution**): output Betti curves.
        """
        num_diag, Xfit = len(X), []
        x_values = np.linspace(self.sample_range[0], self.sample_range[1], self.resolution)
        step_x = x_values[1] - x_values[0]

        for i in range(num_diag):

            diagram, num_pts_in_diag = X[i], X[i].shape[0]

            bc =  np.zeros(self.resolution)
            for j in range(num_pts_in_diag):
                [px,py] = diagram[j,:2]
                min_idx = np.clip(np.ceil((px - self.sample_range[0]) / step_x).astype(int), 0, self.resolution)
                max_idx = np.clip(np.ceil((py - self.sample_range[0]) / step_x).astype(int), 0, self.resolution)
                for k in range(min_idx, max_idx):
                    bc[k] += 1

            Xfit.append(np.reshape(bc,[1,-1]))

        Xfit = np.concatenate(Xfit, 0)

        return Xfit

    def __call__(self, diag):
        """
        Apply BettiCurve on a single persistence diagram and outputs the result.

        Parameters:
            diag (n x 2 numpy array): input persistence diagram.

        Returns:
            numpy array with shape (**resolution**): output Betti curve.
        """
        return self.fit_transform([diag])[0,:]

class Entropy(BaseEstimator, TransformerMixin):
    """
    This is a class for computing persistence entropy. Persistence entropy is a statistic for persistence diagrams inspired from Shannon entropy. This statistic can also be used to compute a feature vector, called the entropy summary function. See https://arxiv.org/pdf/1803.08304.pdf for more details. Note that a previous implementation was contributed by Manuel Soriano-Trigueros.
    """
    def __init__(self, mode="scalar", normalized=True, resolution=100, sample_range=[np.nan, np.nan]):
        """
        Constructor for the Entropy class.

        Parameters:
            mode (string): what entropy to compute: either "scalar" for computing the entropy statistics, or "vector" for computing the entropy summary functions (default "scalar").
            normalized (bool): whether to normalize the entropy summary function (default True). Used only if **mode** = "vector". 
            resolution (int): number of sample for the entropy summary function (default 100). Used only if **mode** = "vector".
            sample_range ([double, double]): minimum and maximum of the entropy summary function domain, of the form [x_min, x_max] (default [numpy.nan, numpy.nan]). It is the interval on which samples will be drawn evenly. If one of the values is numpy.nan, it can be computed from the persistence diagrams with the fit() method. Used only if **mode** = "vector".
        """
        self.mode, self.normalized, self.resolution, self.sample_range = mode, normalized, resolution, sample_range

    def fit(self, X, y=None):
        """
        Fit the Entropy class on a list of persistence diagrams.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
            y (n x 1 array): persistence diagram labels (unused).
        """
        if np.isnan(np.array(self.sample_range)).any():
            pre = DiagramScaler(use=True, scalers=[([0], MinMaxScaler()), ([1], MinMaxScaler())]).fit(X,y)
            [mx,my],[Mx,My] = [pre.scalers[0][1].data_min_[0], pre.scalers[1][1].data_min_[0]], [pre.scalers[0][1].data_max_[0], pre.scalers[1][1].data_max_[0]]
            self.sample_range = np.where(np.isnan(np.array(self.sample_range)), np.array([mx, My]), np.array(self.sample_range))
        return self

    def transform(self, X):
        """
        Compute the entropy for each persistence diagram individually and concatenate the results.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
    
        Returns:
            numpy array with shape (number of diagrams) x (1 if **mode** = "scalar" else **resolution**): output entropy.
        """
        num_diag, Xfit = len(X), []
        x_values = np.linspace(self.sample_range[0], self.sample_range[1], self.resolution)
        step_x = x_values[1] - x_values[0]
        new_X = BirthPersistenceTransform().fit_transform(X)        

        for i in range(num_diag):

            orig_diagram, diagram, num_pts_in_diag = X[i], new_X[i], X[i].shape[0]
            new_diagram = DiagramScaler(use=True, scalers=[([1], MaxAbsScaler())]).fit_transform([diagram])[0]

            if self.mode == "scalar":
                ent = - np.sum( np.multiply(new_diagram[:,1], np.log(new_diagram[:,1])) )
                Xfit.append(np.array([[ent]]))

            else:
                ent = np.zeros(self.resolution)
                for j in range(num_pts_in_diag):
                    [px,py] = orig_diagram[j,:2]
                    min_idx = np.clip(np.ceil((px - self.sample_range[0]) / step_x).astype(int), 0, self.resolution)
                    max_idx = np.clip(np.ceil((py - self.sample_range[0]) / step_x).astype(int), 0, self.resolution)
                    for k in range(min_idx, max_idx):
                        ent[k] += (-1) * new_diagram[j,1] * np.log(new_diagram[j,1])
                    if self.normalized:
                        ent = ent / np.linalg.norm(ent, ord=1)
                    Xfit.append(np.reshape(ent,[1,-1]))

        Xfit = np.concatenate(Xfit, 0)

        return Xfit

    def __call__(self, diag):
        """
        Apply Entropy on a single persistence diagram and outputs the result.

        Parameters:
            diag (n x 2 numpy array): input persistence diagram.

        Returns:
            numpy array with shape (1 if **mode** = "scalar" else **resolution**): output entropy.
        """
        return self.fit_transform([diag])[0,:]

class TopologicalVector(BaseEstimator, TransformerMixin):
    """
    This is a class for computing topological vectors from a list of persistence diagrams. The topological vector associated to a persistence diagram is the sorted vector of a slight modification of the pairwise distances between the persistence diagram points. See https://diglib.eg.org/handle/10.1111/cgf12692 for more details.
    """
    def __init__(self, threshold=10):
        """
        Constructor for the TopologicalVector class.

        Parameters:
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
            numpy array with shape (number of diagrams) x (**threshold**): output topological vectors.
        """
        if self.threshold == -1:
            thresh = np.array([X[i].shape[0] for i in range(len(X))]).max()
        else:
            thresh = self.threshold

        num_diag = len(X)
        Xfit = np.zeros([num_diag, thresh])

        for i in range(num_diag):

            diagram, num_pts_in_diag = X[i], X[i].shape[0]
            pers = 0.5 * (diagram[:,1]-diagram[:,0])
            min_pers = np.minimum(pers,np.transpose(pers))
            distances = DistanceMetric.get_metric("chebyshev").pairwise(diagram)
            vect = np.flip(np.sort(np.triu(np.minimum(distances, min_pers)), axis=None), 0)
            dim = min(len(vect), thresh)
            Xfit[i, :dim] = vect[:dim]

        return Xfit

    def __call__(self, diag):
        """
        Apply TopologicalVector on a single persistence diagram and outputs the result.

        Parameters:
            diag (n x 2 numpy array): input persistence diagram.

        Returns:
            numpy array with shape (**threshold**): output topological vector.
        """
        return self.fit_transform([diag])[0,:]

class ComplexPolynomial(BaseEstimator, TransformerMixin):
    """
    This is a class for computing complex polynomials from a list of persistence diagrams. The persistence diagram points are seen as the roots of some complex polynomial, whose coefficients are returned in a complex vector. See https://link.springer.com/chapter/10.1007%2F978-3-319-23231-7_27 for more details.
    """
    def __init__(self, polynomial_type="R", threshold=10):
        """
        Constructor for the ComplexPolynomial class.

        Parameters:
           polynomial_type (char): either "R", "S" or "T" (default "R"). Type of complex polynomial that is going to be computed (explained in https://link.springer.com/chapter/10.1007%2F978-3-319-23231-7_27).
           threshold (int): number of coefficients (default 10). This is the dimension of the complex vector of coefficients, i.e. the number of coefficients corresponding to the largest degree terms of the polynomial. If -1, this threshold is computed from the list of persistence diagrams by considering the one with the largest number of points and using the dimension of its corresponding complex vector of coefficients as threshold. 
        """
        self.threshold, self.polynomial_type = threshold, polynomial_type

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
            numpy array with shape (number of diagrams) x (**threshold**): output complex vectors of coefficients.
        """
        if self.threshold == -1:
            thresh = np.array([X[i].shape[0] for i in range(len(X))]).max()
        else:
            thresh = self.threshold

        Xfit = np.zeros([len(X), thresh]) + 1j * np.zeros([len(X), thresh])
        for d in range(len(X)):
            D, N = X[d], X[d].shape[0]
            if self.polynomial_type == "R":
                roots = D[:,0] + 1j * D[:,1]
            elif self.polynomial_type == "S":
                alpha = np.linalg.norm(D, axis=1)
                alpha = np.where(alpha==0, np.ones(N), alpha)
                roots = np.multiply( np.multiply(  (D[:,0]+1j*D[:,1]), (D[:,1]-D[:,0])  ), 1./(np.sqrt(2)*alpha) )
            elif self.polynomial_type == "T":
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

    def __call__(self, diag):
        """
        Apply ComplexPolynomial on a single persistence diagram and outputs the result.

        Parameters:
            diag (n x 2 numpy array): input persistence diagram.

        Returns:
            numpy array with shape (**threshold**): output complex vector of coefficients.
        """
        return self.fit_transform([diag])[0,:]

def _lapl_contrast(measure, centers, inertias):
    """contrast function for vectorising `measure` in ATOL"""
    return np.exp(-pairwise.pairwise_distances(measure, Y=centers) / inertias)

def _gaus_contrast(measure, centers, inertias):
    """contrast function for vectorising `measure` in ATOL"""
    return np.exp(-pairwise.pairwise_distances(measure, Y=centers)**2 / inertias**2)

def _indicator_contrast(diags, centers, inertias):
    """contrast function for vectorising `measure` in ATOL"""
    robe_curve = np.clip(2-pairwise.pairwise_distances(diags, Y=centers)/inertias, 0, 1)
    return robe_curve

def _cloud_weighting(measure):
    """automatic uniform weighting with mass 1 for `measure` in ATOL"""
    return np.ones(shape=measure.shape[0])

def _iidproba_weighting(measure):
    """automatic uniform weighting with mass 1/N for `measure` in ATOL"""
    return np.ones(shape=measure.shape[0]) / measure.shape[0]

class Atol(BaseEstimator, TransformerMixin):
    """
    This class allows to vectorise measures (e.g. point clouds, persistence diagrams, etc) after a quantisation step.

    ATOL paper: :cite:`royer2019atol`

    Example
    --------
    >>> from sklearn.cluster import KMeans
    >>> from gudhi.representations.vector_methods import Atol
    >>> import numpy as np
    >>> a = np.array([[1, 2, 4], [1, 4, 0], [1, 0, 4]])
    >>> b = np.array([[4, 2, 0], [4, 4, 0], [4, 0, 2]])
    >>> c = np.array([[3, 2, -1], [1, 2, -1]])
    >>> atol_vectoriser = Atol(quantiser=KMeans(n_clusters=2, random_state=202006))
    >>> atol_vectoriser.fit(X=[a, b, c]).centers
    array([[ 2.        ,  0.66666667,  3.33333333],
           [ 2.6       ,  2.8       , -0.4       ]])
    >>> atol_vectoriser(a)
    array([1.18168665, 0.42375966])
    >>> atol_vectoriser(c)
    array([0.02062512, 1.25157463])
    >>> atol_vectoriser.transform(X=[a, b, c])
    array([[1.18168665, 0.42375966],
           [0.29861028, 1.06330156],
           [0.02062512, 1.25157463]])
    """
    def __init__(self, quantiser, weighting_method="cloud", contrast="gaussian"):
        """
        Constructor for the Atol measure vectorisation class.

        Parameters:
            quantiser (Object): Object with `fit` (sklearn API consistent) and `cluster_centers` and `n_clusters`
                attributes, e.g. sklearn.cluster.KMeans. It will be fitted when the Atol object function `fit` is called.
            weighting_method (string): constant generic function for weighting the measure points
                choose from {"cloud", "iidproba"}
                (default: constant function, i.e. the measure is seen as a point cloud by default).
                This will have no impact if weights are provided along with measures all the way: `fit` and `transform`.
            contrast (string): constant function for evaluating proximity of a measure with respect to centers
                choose from {"gaussian", "laplacian", "indicator"}
                (default: gaussian contrast function, see page 3 in the ATOL paper).
        """
        self.quantiser = quantiser
        self.contrast = {
            "gaussian": _gaus_contrast,
            "laplacian": _lapl_contrast,
            "indicator": _indicator_contrast,
        }.get(contrast, _gaus_contrast)
        self.weighting_method = {
            "cloud"   : _cloud_weighting,
            "iidproba": _iidproba_weighting,
        }.get(weighting_method, _cloud_weighting)

    def fit(self, X, y=None, sample_weight=None):
        """
        Calibration step: fit centers to the sample measures and derive inertias between centers.

        Parameters:
            X (list N x d numpy arrays): input measures in R^d from which to learn center locations and inertias
                (measures can have different N).
            y: Ignored, present for API consistency by convention.
            sample_weight (list of numpy arrays): weights for each measure point in X, optional.
                If None, the object's weighting_method will be used.

        Returns:
            self
        """
        if not hasattr(self.quantiser, 'fit'):
            raise TypeError("quantiser %s has no `fit` attribute." % (self.quantiser))
        if sample_weight is None:
            sample_weight = np.concatenate([self.weighting_method(measure) for measure in X])

        measures_concat = np.concatenate(X)
        self.quantiser.fit(X=measures_concat, sample_weight=sample_weight)
        self.centers = self.quantiser.cluster_centers_
        dist_centers = pairwise.pairwise_distances(self.centers)
        np.fill_diagonal(dist_centers, np.inf)
        self.inertias = np.min(dist_centers, axis=0)/2
        return self

    def __call__(self, measure, sample_weight=None):
        """
        Apply measure vectorisation on a single measure.

        Parameters:
            measure (n x d numpy array): input measure in R^d.

        Returns:
            numpy array in R^self.quantiser.n_clusters.
        """
        if sample_weight is None:
            sample_weight = self.weighting_method(measure)
        return np.sum(sample_weight * self.contrast(measure, self.centers, self.inertias.T).T, axis=1)

    def transform(self, X, sample_weight=None):
        """
        Apply measure vectorisation on a list of measures.

        Parameters:
            X (list N x d numpy arrays): input measures in R^d from which to learn center locations and inertias
                (measures can have different N).
            sample_weight (list of numpy arrays): weights for each measure point in X, optional.
                If None, the object's weighting_method will be used.

        Returns:
            numpy array with shape (number of measures) x (self.quantiser.n_clusters).
        """
        if sample_weight is None:
            sample_weight = [self.weighting_method(measure) for measure in X]
        return np.stack([self(measure, sample_weight=weight) for measure, weight in zip(X, sample_weight)])
