# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Mathieu Carrière, Martin Royer, Gard Spreemann
#
# Copyright (C) 2018-2020 Inria
#
# Modification(s):
#   - 2020/06 Martin: ATOL integration
#   - 2020/12 Gard: A more flexible Betti curve class capable of computing exact curves.
#   - 2021/11 Vincent Rouvreau: factorize _automatic_sample_range

import numpy as np
from scipy.spatial.distance import cdist
from sklearn.base          import BaseEstimator, TransformerMixin
from sklearn.exceptions    import NotFittedError
from sklearn.preprocessing import MinMaxScaler, MaxAbsScaler
from sklearn.metrics       import pairwise
from sklearn.cluster import KMeans

try:
    # New location since 1.0
    from sklearn.metrics     import DistanceMetric
except ImportError:
    # Will be removed in 1.3
    from sklearn.neighbors     import DistanceMetric

from .preprocessing import DiagramScaler, BirthPersistenceTransform, _maybe_fit_transform

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
            if all(len(d) == 0 for d in X):
                self.im_range_fixed_ = self.im_range
            else:
                new_X = BirthPersistenceTransform().fit_transform(X)
                pre = DiagramScaler(use=True, scalers=[([0], MinMaxScaler()), ([1], MinMaxScaler())]).fit(new_X,y)
                [mx,my],[Mx,My] = [pre.scalers[0][1].data_min_[0], pre.scalers[1][1].data_min_[0]], [pre.scalers[0][1].data_max_[0], pre.scalers[1][1].data_max_[0]]
                self.im_range_fixed_ = np.where(np.isnan(np.array(self.im_range)), np.array([mx, Mx, my, My]), np.array(self.im_range))
        else:
            self.im_range_fixed_ = self.im_range
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

            x_values, y_values = np.linspace(self.im_range_fixed_[0], self.im_range_fixed_[1], self.resolution[0]), np.linspace(self.im_range_fixed_[2], self.im_range_fixed_[3], self.resolution[1])
            Xs, Ys = np.tile((diagram[:,0][:,np.newaxis,np.newaxis]-x_values[np.newaxis,np.newaxis,:]),[1,self.resolution[1],1]), np.tile(diagram[:,1][:,np.newaxis,np.newaxis]-y_values[np.newaxis,:,np.newaxis],[1,1,self.resolution[0]])
            image = np.tensordot(w, np.exp((-np.square(Xs)-np.square(Ys))/(2*np.square(self.bandwidth)))/(np.square(self.bandwidth)*2*np.pi), 1)

            Xfit.append(image.reshape(1,-1))

        Xfit = np.concatenate(Xfit, 0)

        return Xfit

    def __call__(self, diag):
        """
        Apply PersistenceImage on a single persistence diagram and outputs the result.
        If :func:`fit` hasn't been run, this uses `fit_transform` on a clone of the object and thus does not affect later calls.

        Parameters:
            diag (n x 2 numpy array): input persistence diagram.

        Returns:
            numpy array with shape (number of pixels = **resolution[0]** x **resolution[1]**):: output persistence image.
        """
        return _maybe_fit_transform(self, 'im_range_fixed_', diag)

def _automatic_sample_range(sample_range, X):
        """
        Compute and returns sample range from the persistence diagrams if one of the sample_range values is numpy.nan.

        Parameters:
            sample_range (a numpy array of 2 float): minimum and maximum of all piecewise-linear function domains, of
                the form [x_min, x_max].
            X (list of n x 2 numpy arrays): input persistence diagrams.
            y (n x 1 array): persistence diagram labels (unused).
        """
        nan_in_range = np.isnan(sample_range)
        if nan_in_range.any():
            try:
                pre = DiagramScaler(use=True, scalers=[([0], MinMaxScaler()), ([1], MinMaxScaler())]).fit(X)
                [mx,my] = [pre.scalers[0][1].data_min_[0], pre.scalers[1][1].data_min_[0]]
                [Mx,My] = [pre.scalers[0][1].data_max_[0], pre.scalers[1][1].data_max_[0]]
                return np.where(nan_in_range, np.array([mx, My]), sample_range)
            except ValueError:
                b = np.nanmax([sample_range[0], sample_range[1], -np.inf])
                print(f"Empty list or empty diagrams: sample range is [{b}, {b}]")
                return np.array([b, b])
        return sample_range


def _trim_endpoints(x, are_endpoints_nan):
    if are_endpoints_nan[0]:
        x = x[1:]
    if are_endpoints_nan[1]:
        x = x[:-1]
    return x


def _grid_from_sample_range(self, X):
    sample_range = np.array(self.sample_range)
    self.nan_in_range_ = np.isnan(sample_range)
    self.new_resolution_ = self.resolution
    if not self.keep_endpoints:
        self.new_resolution_ += self.nan_in_range_.sum()
    self.sample_range_fixed_ = _automatic_sample_range(sample_range, X)
    if self.sample_range_fixed_[0] != self.sample_range_fixed_[1]:
        self.grid_ = np.linspace(self.sample_range_fixed_[0], self.sample_range_fixed_[1], self.new_resolution_)
    else:
        print('First value and second value in range are the same: grid is made of resolution copies of this value')
        self.grid_ = np.full(shape=[self.new_resolution_], fill_value=self.sample_range_fixed_[0])
    if not self.keep_endpoints:
        self.grid_ = _trim_endpoints(self.grid_, self.nan_in_range_)


class Landscape(BaseEstimator, TransformerMixin):
    """
    This is a class for computing persistence landscapes from a list of persistence diagrams. A persistence landscape is a collection of 1D piecewise-linear functions computed from the rank function associated to the persistence diagram. These piecewise-linear functions are then sampled evenly on a given range and the corresponding vectors of samples are concatenated and returned. See http://jmlr.org/papers/v16/bubenik15a.html for more details.

    Attributes:
        grid_ (1d array): The grid on which the landscapes are computed.
    """
    def __init__(self, num_landscapes=5, resolution=100, sample_range=[np.nan, np.nan], *, keep_endpoints=False):
        """
        Constructor for the Landscape class.

        Parameters:
            num_landscapes (int): number of piecewise-linear functions to output (default 5).
            resolution (int): number of sample for all piecewise-linear functions (default 100).
            sample_range ([double, double]): minimum and maximum of all piecewise-linear function domains, of the form [x_min, x_max] (default [numpy.nan, numpy.nan]). It is the interval on which samples will be drawn evenly. If one of the values is numpy.nan, it can be computed from the persistence diagrams with the fit() method.
            keep_endpoints (bool): when computing `sample_range`, use the exact extremities (where the value is always 0). This is mostly useful for plotting, the default is to use a slightly smaller range.
        """
        self.num_landscapes, self.resolution, self.sample_range = num_landscapes, resolution, sample_range
        self.keep_endpoints = keep_endpoints

    def fit(self, X, y=None):
        """
        Fit the Landscape class on a list of persistence diagrams: if any of the values in **sample_range** is numpy.nan, replace it with the corresponding value computed on the given list of persistence diagrams.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
            y (n x 1 array): persistence diagram labels (unused).
        """
        _grid_from_sample_range(self, X)
        return self

    def transform(self, X):
        """
        Compute the persistence landscape for each persistence diagram individually and concatenate the results.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
    
        Returns:
            numpy array with shape (number of diagrams) x (number of samples = **num_landscapes** x **resolution**): output persistence landscapes.
        """

        Xfit = []
        x_values = self.grid_
        for diag in X:
            midpoints, heights = (diag[:, 0] + diag[:, 1]) / 2., (diag[:, 1] - diag[:, 0]) / 2.
            tent_functions = np.maximum(heights[None, :] - np.abs(x_values[:, None] - midpoints[None, :]), 0)
            n_points = diag.shape[0]
            # Complete the array with zeros to get the right number of landscapes
            if self.num_landscapes > n_points:
                tent_functions = np.concatenate(
                    [tent_functions, np.zeros((tent_functions.shape[0], self.num_landscapes-n_points))],
                    axis=1
                )
            tent_functions.partition(tent_functions.shape[1]-self.num_landscapes, axis=1)
            landscapes = np.sort(tent_functions[:, -self.num_landscapes:], axis=1)[:, ::-1].T

            landscapes = np.sqrt(2) * np.ravel(landscapes)
            Xfit.append(landscapes)

        return np.stack(Xfit, axis=0)

    def __call__(self, diag):
        """
        Apply Landscape on a single persistence diagram and outputs the result.
        If :func:`fit` hasn't been run, this uses `fit_transform` on a clone of the object and thus does not affect later calls.

        Parameters:
            diag (n x 2 numpy array): input persistence diagram.

        Returns:
            numpy array with shape (number of samples = **num_landscapes** x **resolution**): output persistence landscape.
        """
        return _maybe_fit_transform(self, 'grid_', diag)

class Silhouette(BaseEstimator, TransformerMixin):
    """
    This is a class for computing persistence silhouettes from a list of persistence diagrams. A persistence silhouette is computed by taking a weighted average of the collection of 1D piecewise-linear functions given by the persistence landscapes, and then by evenly sampling this average on a given range. Finally, the corresponding vector of samples is returned. See https://arxiv.org/abs/1312.0308 for more details.

    Attributes:
        grid_ (1d array): The grid on which the silhouette is computed.
    """
    def __init__(self, weight=lambda x: 1, resolution=100, sample_range=[np.nan, np.nan], *, keep_endpoints=False):
        """
        Constructor for the Silhouette class.

        Parameters:
            weight (function): weight function for the persistence diagram points (default constant function, ie lambda x: 1). This function must be defined on 2D points, ie on lists or numpy arrays of the form [p_x,p_y].
            resolution (int): number of samples for the weighted average (default 100).
            sample_range ([double, double]): minimum and maximum for the weighted average domain, of the form [x_min, x_max] (default [numpy.nan, numpy.nan]). It is the interval on which samples will be drawn evenly. If one of the values is numpy.nan, it can be computed from the persistence diagrams with the fit() method.
            keep_endpoints (bool): when computing `sample_range`, use the exact extremities (where the value is always 0). This is mostly useful for plotting, the default is to use a slightly smaller range.
        """
        self.weight, self.resolution, self.sample_range = weight, resolution, sample_range
        self.keep_endpoints = keep_endpoints

    def fit(self, X, y=None):
        """
        Fit the Silhouette class on a list of persistence diagrams: if any of the values in **sample_range** is numpy.nan, replace it with the corresponding value computed on the given list of persistence diagrams.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
            y (n x 1 array): persistence diagram labels (unused).
        """
        _grid_from_sample_range(self, X)
        return self

    def transform(self, X):
        """
        Compute the persistence silhouette for each persistence diagram individually and concatenate the results.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
    
        Returns:
            numpy array with shape (number of diagrams) x (**resolution**): output persistence silhouettes.
        """
        Xfit = []
        x_values = self.grid_

        for diag in X:
            midpoints, heights = (diag[:, 0] + diag[:, 1]) / 2., (diag[:, 1] - diag[:, 0]) / 2.
            weights = np.array([self.weight(pt) for pt in diag])
            total_weight = np.sum(weights)

            tent_functions = np.maximum(heights[None, :] - np.abs(x_values[:, None] - midpoints[None, :]), 0)
            silhouette = np.sum(weights[None, :] / total_weight * tent_functions, axis=1)
            Xfit.append(silhouette * np.sqrt(2))

        return np.stack(Xfit, axis=0)

    def __call__(self, diag):
        """
        Apply Silhouette on a single persistence diagram and outputs the result.
        If :func:`fit` hasn't been run, this uses `fit_transform` on a clone of the object and thus does not affect later calls.

        Parameters:
            diag (n x 2 numpy array): input persistence diagram.

        Returns:
            numpy array with shape (**resolution**): output persistence silhouette.
        """
        return _maybe_fit_transform(self, 'grid_', diag)


class BettiCurve(BaseEstimator, TransformerMixin):
    """
    Compute Betti curves from persistence diagrams. There are several modes of operation: with a given resolution (with or without a sample_range), with a predefined grid, and with none of the previous. With a predefined grid, the class computes the Betti numbers at those grid points. Without a predefined grid, if the resolution is set to None, it can be fit to a list of persistence diagrams and produce a grid that consists of (at least) the filtration values at which at least one of those persistence diagrams changes Betti numbers, and then compute the Betti numbers at those grid points. In the latter mode, the exact Betti curve is computed for the entire real line. Otherwise, if the resolution is given, the Betti curve is obtained by sampling evenly using either the given sample_range or based on the persistence diagrams.

    Examples
    --------
    If pd is a persistence diagram and xs is a nonempty grid of finite values such that xs[0] >= pd.min(), then the results of:

    >>> bc = BettiCurve(predefined_grid=xs) # doctest: +SKIP
    >>> result = bc(pd) # doctest: +SKIP

    and

    >>> from scipy.interpolate import interp1d # doctest: +SKIP
    >>> bc = BettiCurve(resolution=None, predefined_grid=None) # doctest: +SKIP
    >>> bettis = bc.fit_transform([pd]) # doctest: +SKIP
    >>> interp = interp1d(bc.grid_, bettis[0, :], kind="previous", fill_value="extrapolate") # doctest: +SKIP
    >>> result = np.array(interp(xs), dtype=int) # doctest: +SKIP

    are the same.

    Attributes
    ----------
    grid_ : 1d array
        The grid on which the Betti numbers are computed. If predefined_grid was specified, `grid_` will always be that grid, independently of data. If not and resolution is None, the grid is fitted to capture all filtration values at which the Betti numbers change.
    """

    def __init__(self, resolution=100, sample_range=[np.nan, np.nan], predefined_grid=None, *, keep_endpoints=False):
        """
        Constructor for the BettiCurve class.

        Parameters:
            resolution (int): number of samples for the piecewise-constant function (default 100), or None for the exact curve.
            sample_range ([double, double]): minimum and maximum of the piecewise-constant function domain, of the form [x_min, x_max] (default [numpy.nan, numpy.nan]). It is the interval on which samples will be drawn evenly. If one of the values is numpy.nan, it can be computed from the persistence diagrams with the fit() method.
            predefined_grid (1d array or None, default=None): Predefined filtration grid points at which to compute the Betti curves. Must be strictly ordered. Infinities are ok. If None (default), and resolution is given, the grid will be uniform from x_min to x_max in 'resolution' steps, otherwise a grid will be computed that captures all changes in Betti numbers in the provided data.
            keep_endpoints (bool): when computing `sample_range` (fixed `resolution`, no `predefined_grid`), use the exact extremities. This is mostly useful for plotting, the default is to use a slightly smaller range.
        """

        if (predefined_grid is not None) and (not isinstance(predefined_grid, np.ndarray)):
            raise ValueError("Expected predefined_grid as array or None.")

        self.predefined_grid = predefined_grid
        self.resolution = resolution
        self.sample_range = sample_range
        self.keep_endpoints = keep_endpoints

    def is_fitted(self):
        return hasattr(self, "grid_")

    def fit(self, X, y = None):
        """
        Fit the BettiCurve class on a list of persistence diagrams: if any of the values in **sample_range** is numpy.nan, replace it with the corresponding value computed on the given list of persistence diagrams. When no predefined grid is provided and resolution set to None, compute a filtration grid that captures all changes in Betti numbers for all the given persistence diagrams.

        Parameters:
            X (list of 2d arrays): Persistence diagrams.
            y (None): Ignored.
        """

        if self.predefined_grid is None:
            if self.resolution is None: # Flexible/exact version
                self.grid_ = np.unique(np.concatenate([pd.ravel() for pd in X] + [[-np.inf]], axis=0)) 
            else:
                _grid_from_sample_range(self, X)
        else:
            self.grid_ = self.predefined_grid # Get the predefined grid from user

        return self

    def transform(self, X):
        """
        Compute Betti curves.

        Parameters:
            X (list of 2d arrays): Persistence diagrams.

        Returns:
            `len(X).len(self.grid_)` array of ints: Betti numbers of the given persistence diagrams at the grid points given in `self.grid_`
        """

        if not self.is_fitted():
            raise NotFittedError("Not fitted.")

        N = len(X)

        if N == 0:

            print("Empty list: output has shape [0, len(grid)]")
            return np.zeros((N, len(self.grid_)))
            
        else:

            events = np.concatenate([pd.ravel(order="F") for pd in X], axis=0)
            sorting = np.argsort(events)
            offsets = np.zeros(1 + N, dtype=int)
            for i in range(0, N):
                offsets[i+1] = offsets[i] + 2*X[i].shape[0]
            starts = offsets[0:N]
            ends = offsets[1:N + 1] - 1

            bettis = [[0] for i in range(0, N)]

            i = 0
            for x in self.grid_:
                while i < len(sorting) and events[sorting[i]] <= x:
                    j = np.searchsorted(ends, sorting[i])
                    delta = 1 if sorting[i] - starts[j] < len(X[j]) else -1
                    bettis[j][-1] += delta
                    i += 1
                for k in range(0, N):
                    bettis[k].append(bettis[k][-1])
    
            return np.array(bettis, dtype=int)[:, 0:-1]

    def fit_transform(self, X, y = None):
        """
        The result is the same as fit(X) followed by transform(X), but potentially faster.
        """

        if self.predefined_grid is None and self.resolution is None:

            N = len(X)

            if sum([len(x) for x in X]) == 0:
                print("Empty list or empty diagrams: evaluation grid only contains -infinity and output contains only zeros")
                self.grid_ = np.array([-np.inf])
                return np.zeros((N, 1))

            else:

                events = np.concatenate([pd.ravel(order="F") for pd in X], axis=0)
                sorting = np.argsort(events)
                offsets = np.zeros(1 + N, dtype=int)
                for i in range(0, N):
                    offsets[i+1] = offsets[i] + 2*X[i].shape[0]
                starts = offsets[0:N]
                ends = offsets[1:N + 1] - 1

                xs = [-np.inf]
                bettis = [[0] for i in range(0, N)]

                for i in sorting:
                    j = np.searchsorted(ends, i)
                    delta = 1 if i - starts[j] < len(X[j]) else -1
                    if events[i] == xs[-1]:
                        bettis[j][-1] += delta
                    else:
                        xs.append(events[i])
                        for k in range(0, j):
                            bettis[k].append(bettis[k][-1])
                        bettis[j].append(bettis[j][-1] + delta)
                        for k in range(j+1, N):
                            bettis[k].append(bettis[k][-1])

                self.grid_ = np.array(xs)
                return np.array(bettis, dtype=int)

        else:

            return self.fit(X).transform(X)

    def __call__(self, diag):
        """
        Shorthand for transform on a single persistence diagram.
        If :func:`fit` hasn't been run, this uses `fit_transform` on a clone of the object and thus does not affect later calls.
        """
        return _maybe_fit_transform(self, 'grid_', diag)



class Entropy(BaseEstimator, TransformerMixin):
    """
    This is a class for computing persistence entropy. Persistence entropy is a statistic for persistence diagrams inspired from Shannon entropy. This statistic can also be used to compute a feature vector, called the entropy summary function. See https://arxiv.org/pdf/1803.08304.pdf for more details. Note that a previous implementation was contributed by Manuel Soriano-Trigueros.

    Attributes:
        grid_ (1d array): In vector mode, the grid on which the entropy summary function is computed.
    """
    def __init__(self, mode="scalar", normalized=True, resolution=100, sample_range=[np.nan, np.nan], *, keep_endpoints=False):
        """
        Constructor for the Entropy class.

        Parameters:
            mode (string): what entropy to compute: either "scalar" for computing the entropy statistics, or "vector" for computing the entropy summary functions (default "scalar").
            normalized (bool): whether to normalize the entropy summary function (default True). Used only if **mode** = "vector". 
            resolution (int): number of sample for the entropy summary function (default 100). Used only if **mode** = "vector".
            sample_range ([double, double]): minimum and maximum of the entropy summary function domain, of the form [x_min, x_max] (default [numpy.nan, numpy.nan]). It is the interval on which samples will be drawn evenly. If one of the values is numpy.nan, it can be computed from the persistence diagrams with the fit() method. Used only if **mode** = "vector".
            keep_endpoints (bool): when computing `sample_range`, use the exact extremities. This is mostly useful for plotting, the default is to use a slightly smaller range.
        """
        self.mode, self.normalized, self.resolution, self.sample_range = mode, normalized, resolution, sample_range
        self.keep_endpoints = keep_endpoints

    def fit(self, X, y=None):
        """
        Fit the Entropy class on a list of persistence diagrams.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
            y (n x 1 array): persistence diagram labels (unused).
        """
        if self.mode == "vector":
            _grid_from_sample_range(self, X)
            if self.sample_range_fixed_[0] != self.sample_range_fixed_[1]:
                self.step_ = self.grid_[1] - self.grid_[0]
            else:
                self.step_ = 0.
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
        new_X = BirthPersistenceTransform().fit_transform(X)        

        for i in range(num_diag):
            orig_diagram, new_diagram, num_pts_in_diag = X[i], new_X[i], X[i].shape[0]
                
            p = new_diagram[:,1]
            p = p/np.sum(p)
            if self.mode == "scalar":
                ent = -np.dot(p, np.log(p))
                Xfit.append(np.array([[ent]]))
            else:
                ent = np.zeros(self.resolution)
                for j in range(num_pts_in_diag):
                    [px,py] = orig_diagram[j,:2]
                    min_idx = np.clip(np.ceil((px - self.sample_range_fixed_[0]) / self.step_).astype(int), 0, self.resolution)
                    max_idx = np.clip(np.ceil((py - self.sample_range_fixed_[0]) / self.step_).astype(int), 0, self.resolution)
                    ent[min_idx:max_idx]-=p[j]*np.log(p[j])
                if self.normalized:
                    ent = ent / np.linalg.norm(ent, ord=1)
                Xfit.append(np.reshape(ent,[1,-1]))

        Xfit = np.concatenate(Xfit, axis=0)
        return Xfit

    def __call__(self, diag):
        """
        Apply Entropy on a single persistence diagram and outputs the result.
        If :func:`fit` hasn't been run, this uses `fit_transform` on a clone of the object and thus does not affect later calls.

        Parameters:
            diag (n x 2 numpy array): input persistence diagram.

        Returns:
            numpy array with shape (1 if **mode** = "scalar" else **resolution**): output entropy.
        """
        return _maybe_fit_transform(self, 'grid_', diag)

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
            # Works fine with sklearn 1.0, but an ValueError exception is thrown on past versions
            try:
                distances = DistanceMetric.get_metric("chebyshev").pairwise(diagram)
            except ValueError:
                # Empty persistence diagram case - https://github.com/GUDHI/gudhi-devel/issues/507
                assert len(diagram) == 0
                distances = np.empty(shape = [0, 0])
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
        return self.transform([diag])[0,:]

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
        return self.transform([diag])[0,:]

def _lapl_contrast(measure, centers, inertias):
    """contrast function for vectorising `measure` in ATOL
    we use cdist so as to accept 'inf' values instead of raising ValueError with sklearn.pairwise"""
    return np.exp(-cdist(XA=measure, XB=centers) / inertias)

def _gaus_contrast(measure, centers, inertias):
    """contrast function for vectorising `measure` in ATOL
    we use cdist so as to accept 'inf' values instead of raising ValueError with sklearn.pairwise"""
    return np.exp(-cdist(XA=measure, XB=centers, metric="sqeuclidean") / inertias**2)

def _indicator_contrast(diags, centers, inertias):
    """contrast function for vectorising `measure` in ATOL
    we use cdist so as to accept 'inf' values instead of raising ValueError with sklearn.pairwise"""
    robe_curve = np.clip(2-cdist(XA=diags, XB=centers)/inertias, 0, 1)
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
    >>> atol_vectoriser = Atol(quantiser=KMeans(n_clusters=2, random_state=202006, n_init=10))
    >>> atol_vectoriser.fit(X=[a, b, c]).centers
    array([[ 2.6       ,  2.8       , -0.4       ],
           [ 2.        ,  0.66666667,  3.33333333]])
    >>> atol_vectoriser(a)
    array([0.42375966, 1.18168665])
    >>> atol_vectoriser(c)
    array([1.25157463, 0.02062512])
    >>> atol_vectoriser.transform(X=[a, b, c])
    array([[0.42375966, 1.18168665],
           [1.06330156, 0.29861028],
           [1.25157463, 0.02062512]])
    """
    # Note the example above must be up to date with the one in tests called test_atol_doc
    def __init__(
            self,
            quantiser=KMeans(n_clusters=2, n_init="auto"),
            weighting_method="cloud",
            contrast="gaussian"
    ):
        """
        Constructor for the Atol measure vectorisation class.

        Parameters:
            quantiser (Object): Object with `fit` (sklearn API consistent) and `cluster_centers` and `n_clusters`
                attributes, e.g. `sklearn.cluster.KMeans`. It will be fitted when the Atol object function `fit` is
                called. Users are encouraged to provide their own quantiser, and in particular increase the number
                of clusters.
            weighting_method (string): constant generic function for weighting the measure points
                choose from {"cloud", "iidproba"}
                (default: constant function, i.e. the measure is seen as a point cloud by default).
                This will have no impact if weights are provided along with measures all the way: `fit` and `transform`.
            contrast (string): constant function for evaluating proximity of a measure with respect to centers
                choose from {"gaussian", "laplacian", "indicator"}
                (default: gaussian contrast function, see page 3 in the ATOL paper).
        """
        self.quantiser = quantiser
        self.contrast = contrast
        self.weighting_method = weighting_method
        self._running_transform_names = ""

    def get_contrast(self):
        return {
            "gaussian": _gaus_contrast,
            "laplacian": _lapl_contrast,
            "indicator": _indicator_contrast,
        }.get(self.contrast, _gaus_contrast)

    def get_weighting_method(self):
        return {
            "cloud"   : _cloud_weighting,
            "iidproba": _iidproba_weighting,
        }.get(self.weighting_method, _cloud_weighting)

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

        # In fitting we remove infinite death time points so that every center is finite
        X = [dgm[~np.isinf(dgm).any(axis=1), :] for dgm in X]

        if sample_weight is None:
            sample_weight = [self.get_weighting_method()(measure) for measure in X]

        measures_concat = np.concatenate(X)
        weights_concat = np.concatenate(sample_weight)

        self.quantiser.fit(X=measures_concat, sample_weight=weights_concat)

        self.centers = self.quantiser.cluster_centers_
        # Hack, but some people are unhappy if the order depends on the version of sklearn
        self.centers = self.centers[np.lexsort(self.centers.T)]
        if self.quantiser.n_clusters == 1:
            dist_centers = pairwise.pairwise_distances(measures_concat)
            np.fill_diagonal(dist_centers, 0)
            best_inertia = np.max(dist_centers)/2 if np.max(dist_centers)/2 > 0 else 1
            self.inertias = np.array([best_inertia])
        else:
            dist_centers = pairwise.pairwise_distances(self.centers)
            dist_centers[dist_centers == 0] = np.inf
            self.inertias = np.min(dist_centers, axis=0)/2
        return self

    def __call__(self, measure, sample_weight=None):
        """
        Apply measure vectorisation on a single measure. Only available after `fit` has been called.

        Parameters:
            measure (n x d numpy array): input measure in R^d.

        Returns:
            numpy array in R^self.quantiser.n_clusters.
        """
        if sample_weight is None:
            sample_weight = self.get_weighting_method()(measure)
        return np.sum(sample_weight * self.get_contrast()(measure, self.centers, self.inertias.T).T, axis=1)

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
            sample_weight = [self.get_weighting_method()(measure) for measure in X]
        self._running_transform_names = [f"Atol Center {i + 1}" for i in range(self.quantiser.n_clusters)]
        return np.stack([self(measure, sample_weight=weight) for measure, weight in zip(X, sample_weight)])

    def get_feature_names_out(self):
        return self._running_transform_names
