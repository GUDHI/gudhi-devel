# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Mathieu Carrière, Vincent Rouvreau
#
# Copyright (C) 2018-2019 Inria
#
# Modification(s):
#   - 2021/10 Vincent Rouvreau: Add DimensionSelector
#   - YYYY/MM Author: Description of the modification

import numpy as np
from sklearn.base          import BaseEstimator, TransformerMixin
from sklearn.preprocessing import StandardScaler


#############################################
# Utils #####################################
#############################################

def _maybe_fit_transform(obj, attr, diag):
    """
    In __call__, use transform on the object itself if it has been fitted,
    otherwise fit_transform on a clone of the object so it doesn't affect future calls.
    """
    if hasattr(obj, attr):
        result = obj.transform([diag])
    else:
        result = obj.__class__(**obj.get_params()).fit_transform([diag])
    return result[0]

class Clamping(BaseEstimator, TransformerMixin):
    """
    This is a class for clamping a list of values. It is not meant to be called directly on (a list of) persistence diagrams, but it is rather meant to be used as a parameter for the DiagramScaler class. As such it has the same methods and purpose as common scalers from sklearn.preprocessing such as MinMaxScaler, RobustScaler, StandardScaler, etc. A typical use would be for instance if you want to clamp abscissae or ordinates (or both) of persistence diagrams within a pre-defined interval.
    """
    def __init__(self, minimum=-np.inf, maximum=np.inf):
        """
        Constructor for the Clamping class.

        Parameters:
            limit (float): clamping value (default np.inf).
        """
        self.minimum = minimum
        self.maximum = maximum

    def fit(self, X, y=None):
        """
        Fit the Clamping class on a list of values (this function actually does nothing but is useful when Clamping is included in a scikit-learn Pipeline).

        Parameters:
            X (numpy array of size n): input values.
            y (n x 1 array): value labels (unused).
        """
        return self

    def transform(self, X):
        """
        Clamp list of values.

        Parameters:
            X (numpy array of size n): input list of values.

        Returns:
            numpy array of size n: output list of values.
        """
        Xfit = np.clip(X, self.minimum, self.maximum)
        return Xfit


#############################################
# Preprocessing #############################
#############################################

class BirthPersistenceTransform(BaseEstimator, TransformerMixin):
    """
    This is a class for the affine transformation (x,y) -> (x,y-x) to be applied on persistence diagrams.
    """
    def __init__(self):
        """
        Constructor for BirthPersistenceTransform class.
        """
        return None

    def fit(self, X, y=None):
        """
        Fit the BirthPersistenceTransform class on a list of persistence diagrams (this function actually does nothing but is useful when BirthPersistenceTransform is included in a scikit-learn Pipeline).

        Parameters:
            X (list of n x 2 numpy array): input persistence diagrams.
            y (n x 1 array): persistence diagram labels (unused).
        """
        return self

    def transform(self, X):
        """
        Apply the BirthPersistenceTransform function on the persistence diagrams.

        Parameters:
            X (list of n x 2 numpy array): input persistence diagrams.

        Returns:
            list of n x 2 numpy array: transformed persistence diagrams.
        """
        Xfit = []
        for diag in X:
            #new_diag = np.empty(diag.shape)
            #np.copyto(new_diag, diag)
            new_diag = np.copy(diag)
            new_diag[:,1] = new_diag[:,1] - new_diag[:,0]
            Xfit.append(new_diag)
        return Xfit

    def __call__(self, diag):
        """
        Apply BirthPersistenceTransform on a single persistence diagram and outputs the result.

        Parameters:
            diag (n x 2 numpy array): input persistence diagram.

        Returns:
            n x 2 numpy array: transformed persistence diagram.
        """
        return self.transform([diag])[0]

class DiagramScaler(BaseEstimator, TransformerMixin):
    """
    This is a class for preprocessing persistence diagrams with a given list of scalers, such as those included in scikit-learn.
    """
    def __init__(self, use=False, scalers=[]):
        """
        Constructor for the DiagramScaler class.

        Parameters:
            use (bool): whether to use the class or not (default False).
            scalers (list of classes): list of scalers to be fit on the persistence diagrams (default []). Each element of the list is a tuple with two elements: the first one is a list of coordinates, and the second one is a scaler (i.e. a class with fit() and transform() methods) that is going to be applied to these coordinates. Common scalers can be found in the scikit-learn library (such as MinMaxScaler for instance).
        """
        self.scalers  = scalers
        self.use      = use

    def fit(self, X, y=None):
        """
        Fit the DiagramScaler class on a list of persistence diagrams: persistence diagrams are concatenated in a big numpy array, and scalers are fit (by calling their fit() method) on their corresponding coordinates in this big array.

        Parameters:
            X (list of n x 2 or n x 1 numpy arrays): input persistence diagrams.
            y (n x 1 array): persistence diagram labels (unused).
        """
        self.is_fitted_ = True
        if self.use:
            if len(X) == 1:
                P = X[0]
            else:
                P = np.concatenate(X,0)
            for (indices, scaler) in self.scalers:
                scaler.fit(np.reshape(P[:,indices], [-1, 1]))
        return self

    def transform(self, X):
        """
        Apply the DiagramScaler function on the persistence diagrams. The fitted scalers are applied (by calling their transform() method) to their corresponding coordinates in each persistence diagram individually.  

        Parameters:
            X (list of n x 2 or n x 1 numpy arrays): input persistence diagrams.

        Returns:
            list of n x 2 or n x 1 numpy arrays: transformed persistence diagrams.
        """
        Xfit = [np.copy(d) for d in X]
        if self.use:
            for i in range(len(Xfit)):
                if Xfit[i].shape[0] > 0:
                    for (indices, scaler) in self.scalers:
                        for I in indices:
                            Xfit[i][:,I] = np.squeeze(scaler.transform(np.reshape(Xfit[i][:,I], [-1,1])))
        return Xfit

    def __call__(self, diag):
        """
        Apply DiagramScaler on a single persistence diagram and outputs the result.
        If :func:`fit` hasn't been run, this uses `fit_transform` on a clone of the object and thus does not affect later calls.

        Parameters:
            diag (n x 2 numpy array): input persistence diagram.

        Returns:
            n x 2 numpy array: transformed persistence diagram.
        """
        return _maybe_fit_transform(self, 'is_fitted_', diag)

class Padding(BaseEstimator, TransformerMixin):
    """
    This is a class for padding a list of persistence diagrams with dummy points, so that all persistence diagrams end up with the same number of points.
    """
    def __init__(self, use=False):
        """
        Constructor for the Padding class.

        Parameters:
            use (bool): whether to use the class or not (default False).
        """
        self.use = use

    def fit(self, X, y=None):
        """
        Fit the Padding class on a list of persistence diagrams.

        Parameters:
            X (list of n x 2 or n x 1 numpy arrays): input persistence diagrams.
            y (n x 1 array): persistence diagram labels (unused).
        """
        self.max_pts_ = max(len(diag) for diag in X)
        return self

    def transform(self, X):
        """
        Add dummy points to each persistence diagram so that they all have the same cardinality. All points are given an additional coordinate indicating if the point was added after padding (0) or already present before (1).  

        Parameters:
            X (list of n x 2 or n x 1 numpy arrays): input persistence diagrams.

        Returns:
            list of n x 3 or n x 2 numpy arrays: padded persistence diagrams.
        """
        if self.use:
            Xfit, num_diag = [], len(X)
            for diag in X:
                diag_pad = np.pad(diag, ((0,max(0, self.max_pts_ - diag.shape[0])), (0,1)), "constant", constant_values=((0,0),(0,0)))
                diag_pad[:diag.shape[0],2] = np.ones(diag.shape[0])
                Xfit.append(diag_pad)                    
        else:
            Xfit = X
        return Xfit

    def __call__(self, diag):
        """
        Apply Padding on a single persistence diagram and outputs the result.
        If :func:`fit` hasn't been run, this uses `fit_transform` on a clone of the object and thus does not affect later calls.

        Parameters:
            diag (n x 2 numpy array): input persistence diagram.

        Returns:
            n x 2 numpy array: padded persistence diagram.
        """
        return _maybe_fit_transform(self, 'max_pts_', diag)

class ProminentPoints(BaseEstimator, TransformerMixin):
    """
    This is a class for removing points that are close or far from the diagonal in persistence diagrams.  If persistence diagrams are n x 2 numpy arrays (i.e. persistence diagrams with ordinary features), points are ordered and thresholded by distance-to-diagonal. If persistence diagrams are n x 1 numpy arrays (i.e. persistence diagrams with essential features), points are not ordered and thresholded by first coordinate.
    """
    def __init__(self, use=False, num_pts=10, threshold=-1, location="upper"):
        """
        Constructor for the ProminentPoints class.
     
        Parameters:
            use (bool): whether to use the class or not (default False).
            location (string): either "upper" or "lower" (default "upper"). Whether to keep the points that are far away ("upper") or close ("lower") to the diagonal.
            num_pts (int): cardinality threshold (default 10). If location == "upper", keep the top **num_pts** points that are the farthest away from the diagonal. If location == "lower", keep the top **num_pts** points that are the closest to the diagonal. 
            threshold (float): distance-to-diagonal threshold (default -1). If location == "upper", keep the points that are at least at a distance **threshold** from the diagonal. If location == "lower", keep the points that are at most at a distance **threshold** from the diagonal. 
        """
        self.num_pts    = num_pts
        self.threshold  = threshold
        self.use        = use
        self.location   = location

    def fit(self, X, y=None):
        """
        Fit the ProminentPoints class on a list of persistence diagrams (this function actually does nothing but is useful when ProminentPoints is included in a scikit-learn Pipeline).

        Parameters:
            X (list of n x 2 or n x 1 numpy arrays): input persistence diagrams.
            y (n x 1 array): persistence diagram labels (unused).
        """
        return self

    def transform(self, X):
        """
        If location == "upper", first select the top **num_pts** points that are the farthest away from the diagonal, then select and return from these points the ones that are at least at distance **threshold** from the diagonal for each persistence diagram individually. If location == "lower", first select the top **num_pts** points that are the closest to the diagonal, then select and return from these points the ones that are at most at distance **threshold** from the diagonal for each persistence diagram individually.

        Parameters:
            X (list of n x 2 or n x 1 numpy arrays): input persistence diagrams.

        Returns:
            list of n x 2 or n x 1 numpy arrays: thresholded persistence diagrams.
        """
        if self.use:
            Xfit, num_diag = [], len(X)
            for i in range(num_diag):
                diag = X[i]
                if diag.shape[1] >= 2:
                    if diag.shape[0] > 0:
                        pers       = np.abs(diag[:,1] - diag[:,0])
                        idx_thresh = pers >= self.threshold
                        thresh_diag, thresh_pers  = diag[idx_thresh], pers[idx_thresh]
                        sort_index  = np.flip(np.argsort(thresh_pers, axis=None), 0)
                        if self.location == "upper":
                            new_diag = thresh_diag[sort_index[:min(self.num_pts, thresh_diag.shape[0])],:]
                        if self.location == "lower":
                            new_diag = np.concatenate( [ thresh_diag[sort_index[min(self.num_pts, thresh_diag.shape[0]):],:], diag[~idx_thresh] ], axis=0)
                    else:
                        new_diag = diag

                else:
                    if diag.shape[0] > 0:
                        birth      = diag[:,:1]
                        idx_thresh = birth >= self.threshold
                        thresh_diag, thresh_birth  = diag[idx_thresh], birth[idx_thresh]
                        if self.location == "upper":
                            new_diag = thresh_diag[:min(self.num_pts, thresh_diag.shape[0]),:]
                        if self.location == "lower":
                            new_diag = np.concatenate( [ thresh_diag[min(self.num_pts, thresh_diag.shape[0]):,:], diag[~idx_thresh] ], axis=0)
                    else:
                        new_diag = diag

                Xfit.append(new_diag)                    
        else:
            Xfit = X
        return Xfit

    def __call__(self, diag):
        """
        Apply ProminentPoints on a single persistence diagram and outputs the result.

        Parameters:
            diag (n x 2 numpy array): input persistence diagram.

        Returns:
            n x 2 numpy array: thresholded persistence diagram.
        """
        return self.transform([diag])[0]

class DiagramSelector(BaseEstimator, TransformerMixin):
    """
    This is a class for extracting finite or essential points in persistence diagrams.
    """
    def __init__(self, use=False, limit=np.inf, point_type="finite"):
        """
        Constructor for the DiagramSelector class.

        Parameters:
            use (bool): whether to use the class or not (default False).
            limit (float): second coordinate value that is the criterion for being an essential point (default numpy.inf).
            point_type (string): either "finite" or "essential". The type of the points that are going to be extracted.
        """
        self.use, self.limit, self.point_type = use, limit, point_type

    def fit(self, X, y=None):
        """
        Fit the DiagramSelector class on a list of persistence diagrams (this function actually does nothing but is useful when DiagramSelector is included in a scikit-learn Pipeline).

        Parameters:
            X (list of n x 2 or n x 1 numpy arrays): input persistence diagrams.
            y (n x 1 array): persistence diagram labels (unused).
        """
        return self

    def transform(self, X):
        """
        Extract and return the finite or essential points of each persistence diagram individually.

        Parameters:
            X (list of n x 2 or n x 1 numpy arrays): input persistence diagrams.

        Returns:
            list of n x 2 or n x 1 numpy arrays: extracted persistence diagrams.
        """
        if self.use:
            Xfit, num_diag = [], len(X)
            if self.point_type == "finite":
                Xfit = [ diag[diag[:,1] < self.limit] if diag.shape[0] != 0 else diag for diag in X]
            else:
                Xfit = [ diag[diag[:,1] >= self.limit, 0:1] if diag.shape[0] != 0 else diag for diag in X]
        else:
            Xfit = X
        return Xfit

    def __call__(self, diag):
        """
        Apply DiagramSelector on a single persistence diagram and outputs the result.

        Parameters:
            diag (n x 2 numpy array): input persistence diagram.

        Returns:
            n x 2 numpy array: extracted persistence diagram.
        """
        return self.transform([diag])[0]


# Mermaid sequence diagram - https://mermaid-js.github.io/mermaid-live-editor/
# sequenceDiagram
#     USER->>DimensionSelector: fit_transform(<br/>[[array( Hi(X0) ), array( Hj(X0) ), ...],<br/> [array( Hi(X1) ), array( Hj(X1) ), ...],<br/> ...])
#     DimensionSelector->>thread1: _transform([array( Hi(X0) ), array( Hj(X0) )], ...)
#     DimensionSelector->>thread2: _transform([array( Hi(X1) ), array( Hj(X1) )], ...)
#     Note right of DimensionSelector: ...
#     thread1->>DimensionSelector: array( Hn(X0) )
#     thread2->>DimensionSelector: array( Hn(X1) )
#     Note right of DimensionSelector: ...
#     DimensionSelector->>USER: [array( Hn(X0) ), <br/> array( Hn(X1) ), <br/> ...]

class DimensionSelector(BaseEstimator, TransformerMixin):
    """
    This is a class to select persistence diagrams in a specific dimension from its index.
    """

    def __init__(self, index=0):
        """
        Constructor for the DimensionSelector class.

        Parameters:
            index (int): The returned persistence diagrams dimension index. Default value is `0`.
        """
        self.index = index

    def fit(self, X, Y=None):
        """
        Nothing to be done, but useful when included in a scikit-learn Pipeline.
        """
        return self

    def transform(self, X, Y=None):
        """
        Select persistence diagrams from its dimension.

        Parameters:
            X (list of list of tuple): List of list of persistence pairs, i.e.
                `[[array( Hi(X0) ), array( Hj(X0) ), ...], [array( Hi(X1) ), array( Hj(X1) ), ...], ...]` 

        Returns:
            list of tuple:
            Persistence diagrams in a specific dimension. i.e. if `index` was set to `m` and `Hn` is at index `m` of
            the input, it returns `[array( Hn(X0) ), array( Hn(X1), ...]`
        """

        return [persistence[self.index] for persistence in X]
