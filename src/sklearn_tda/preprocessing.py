"""
@author: Mathieu Carriere
All rights reserved
"""

import numpy as np
from sklearn.base          import BaseEstimator, TransformerMixin
from sklearn.preprocessing import StandardScaler

#############################################
# Preprocessing #############################
#############################################

class BirthPersistenceTransform(BaseEstimator, TransformerMixin):
    """
    This is a class for the affine transformation (x,y) -> (x,y-x) to be applied on persistence diagrams. It is a particular scaler for persistence diagram that can be given as input for the DiagramPreprocessor class.
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
            X (list of n x 2 or n x 1 numpy arrays): input persistence diagrams.
            y (n x 1 array): persistence diagram labels (unused).
        """
        return self

    def transform(self, X):
        """
        Apply the BirthPersistenceTransform function on the persistence diagrams.

        Parameters:
            X (list of n x 2 or n x 1 numpy arrays): input persistence diagrams.

        Returns:
            Xfit (list of n x 2 numpy arrays): transformed persistence diagrams.
        """
        Xfit = np.matmul(X, np.array([[1., -1.],[0., 1.]]))
        return Xfit


class DiagramPreprocessor(BaseEstimator, TransformerMixin):
    """
    This is a class for preprocessing persistence diagrams with a given list of scalers, such as those included in scikit-learn.
    """
    def __init__(self, use=False, scalers=[]):
        """
        Constructor for the DiagramPreprocessor class.

        Attributes:
            use (bool): whether to use the class or not (default False).
            scalers (list of classes): list of scalers to be fit on the persistence diagrams (default []). Each element of the list is a tuple with two elements: the first  one is a list of coordinates, and the second one is a scaler (i.e. a class with fit() and transform() methods) that is going to be applied to these coordinates.
        """
        self.scalers = scalers
        self.use     = use

    def fit(self, X, y=None):
        """
        Fit the DiagramPreprocessor class on a list of persistence diagrams. Persistence diagrams are concatenated in a big numpy array, and scalers are fit (by calling their fit() method) on their corresponding coordinates in this big array.

        Parameters:
            X (list of n x 2 or n x 1 numpy arrays): input persistence diagrams.
            y (n x 1 array): persistence diagram labels (unused).
        """
        if self.use:
            if len(X) == 1:
                P = X[0]
            else:
                P = np.concatenate(X,0)
            for (indices, scaler) in self.scalers:
                scaler.fit(P[:,indices])
        return self

    def transform(self, X):
        """
        Apply the DiagramPreprocessor function on the persistence diagrams. The fitted scalers are applied (by calling their transform() method) to their corresponding coordinates in each persistence diagram individually.  

        Parameters:
            X (list of n x 2 or n x 1 numpy arrays): input persistence diagrams.

        Returns:
            Xfit (list of n x 2 or n x 1 numpy arrays): transformed persistence diagrams.
        """
        Xfit = [np.copy(d) for d in X]
        if self.use:
            for i in range(len(Xfit)):
                if Xfit[i].shape[0] > 0:
                    for (indices, scaler) in self.scalers:
                        Xfit[i][:,indices] = scaler.transform(Xfit[i][:,indices])
        return Xfit

class Padding(BaseEstimator, TransformerMixin):
    """
    This is a class for padding a list of persistence diagrams with dummy points, so that all persistence diagrams end up with the same number of points.
    """
    def __init__(self, use=False):
        """
        Constructor for the Padding class.

        Attributes:
            use (bool): whether to use the class or not (default False).
        """
        self.use = use

    def fit(self, X, y=None):
        """
        Fit the Padding class on a list of persistence diagrams (this function actually does nothing but is useful when Padding is included in a scikit-learn Pipeline).

        Parameters:
            X (list of n x 2 or n x 1 numpy arrays): input persistence diagrams.
            y (n x 1 array): persistence diagram labels (unused).
        """
        return self

    def transform(self, X):
        """
        Add dummy points to each persistence diagram so that they all have the same cardinality. All points are given an additional coordinate indicating if the point was added after padding (0) or already present before (1).  

        Parameters:
            X (list of n x 2 or n x 1 numpy arrays): input persistence diagrams.

        Returns:
            Xfit (list of n x 3 or n x 2 numpy arrays): padded persistence diagrams.
        """
        if self.use:
            Xfit, num_diag = [], len(X)
            max_card = max([len(diag) for diag in X])
            for diag in X:
                [num_pts, dim] = diag.shape
                diag_pad = np.zeros([max_card, dim+1])
                diag_pad[:num_pts,:dim] = diag
                diag_pad[:num_pts, dim] = np.ones(num_pts)
                Xfit.append(diag_pad)                    
        else:
            Xfit = X
        return Xfit

class ProminentPoints(BaseEstimator, TransformerMixin):
    """
    This is a class for removing points that are close or far from the diagonal in persistence diagrams.
    """
    def __init__(self, use=False, num_pts=10, threshold=-1, location="upper", point_type="finite"):
        """
        Constructor for the ProminentPoints class.
     
        Attributes:
            use (bool): whether to use the class or not (default False).
            location (string): either "upper" or "lower" (default "upper"). Whether to keep the points that are far away ("upper") or close ("lower") to the diagonal.
            num_pts (int): cardinality threshold (default 10). If location == "upper", keep the top **num_pts** points that are the farthest away from the diagonal. If location == "lower", keep the top **num_pts** points that are the closest to the diagonal. 
            threshold (double): distance-to-diagonal threshold (default -1). If location == "upper", keep the points that are at least at a distance **threshold** from the diagonal. If location == "lower", keep the points that are at most at a distance **threshold** from the diagonal. 
            point_type (string): either "finite" if persistence diagrams are n x 2 numpy arrays, or "essential" if persistence diagrams are n x 1 numpy arrays (default "finite"). If "finite", points are ordered and thresholded by distance-to-diagonal. If "essential", points are ordered and thresholded by first coordinate.
        """
        self.num_pts    = num_pts
        self.threshold  = threshold
        self.use        = use
        self.location   = location
        self.point_type = point_type

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
        If location == "upper", first select the top **num_pts** points that are the farthest away from the diagonal, then select and return from these points the ones that are at least at distance **threshold** from the diagonal for each persistence diagram individually. If location == "lower", first select the top **num_pts** points that are the closest to the diagonal, then select and return from these points the ones that are at most at distance **threshold** from the diagonal for each persistence diagram individually. If point_type == "essential", do the same with first coordinate instead of distance-to-diagonal.

        Parameters:
            X (list of n x 2 or n x 1 numpy arrays): input persistence diagrams.

        Returns:
            Xfit (list of n x 2 or n x 1 numpy arrays): thresholded persistence diagrams.
        """
        if self.use:
            Xfit, num_diag = [], len(X)
            for i in range(num_diag):
                diag = X[i]
                if self.point_type == "finite":
                    if diag.shape[0] > 0:
                        pers       = np.abs(np.matmul(diag[:,:2], [-1., 1.]))
                        idx_thresh = pers >= self.threshold
                        thresh_diag, thresh_pers  = diag[idx_thresh.flatten()], pers[idx_thresh.flatten()]
                        sort_index  = np.flip(np.argsort(thresh_pers, axis=None), 0)
                        if self.location == "upper":
                            new_diag = thresh_diag[sort_index[:min(self.num_pts, thresh_diag.shape[0])],:]
                        if self.location == "lower":
                            new_diag = np.concatenate( [ thresh_diag[sort_index[min(self.num_pts, thresh_diag.shape[0]):],:], diag[~idx_thresh.flatten()] ], axis=0)
                    else:
                        new_diag = diag

                else:
                    if diag.shape[0] > 0:
                        birth      = diag[:,:1]
                        idx_thresh = birth >= self.threshold
                        thresh_diag, thresh_birth  = diag[idx_thresh.flatten()], birth[idx_thresh.flatten()]
                        if self.location == "upper":
                            new_diag = thresh_diag[:min(self.num_pts, thresh_diag.shape[0]),:]
                        if self.location == "lower":
                            new_diag = np.concatenate( [ thresh_diag[min(self.num_pts, thresh_diag.shape[0]):,:], diag[~idx_thresh.flatten()] ], axis=0)
                    else:
                        new_diag = diag

                Xfit.append(new_diag)                    
        else:
            Xfit = X
        return Xfit

class DiagramSelector(BaseEstimator, TransformerMixin):
    """
    This is a class for extracting finite or essential points in persistence diagrams.
    """
    def __init__(self, use=False, limit=np.inf, point_type="finite"):
        """
        Constructor for the DiagramSelector class.

        Attributes:
            use (bool): whether to use the class or not (default False).
            limit (double): second coordinate value that is the criterion for being an essential point (default numpy.inf).
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
            Xfit (list of n x 2 or n x 1 numpy arrays): extracted persistence diagrams.
        """
        if self.use:
            Xfit, num_diag = [], len(X)
            if self.point_type == "finite":
                for i in range(num_diag):
                    diag = X[i]
                    if diag.shape[0] != 0:
                        idx_fin = diag[:,1] != self.limit
                        Xfit.append(diag[idx_fin,:])
                    else:
                        Xfit.append(diag)
            if self.point_type == "essential":
                for i in range(num_diag):
                    diag = X[i]
                    if diag.shape[0] != 0:
                        idx_ess = diag[:,1] == self.limit
                        Xfit.append(np.delete(diag,1,1)[idx_ess,:])
                    else:
                        Xfit.append(np.delete(diag,1,1))
        else:
            Xfit = X
        return Xfit
