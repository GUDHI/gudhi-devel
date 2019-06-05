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

    def __init__(self):
        return None

    def fit(self, X, y=None):
        return self

    def transform(self, X):
        return np.matmul(X, np.array([[1., -1.],[0., 1.]]))


class DiagramPreprocessor(BaseEstimator, TransformerMixin):

    def __init__(self, use=False, scalers=[]):
        self.scalers = scalers
        self.use     = use

    def fit(self, X, y=None):
        if self.use:
            if len(X) == 1:
                P = X[0]
            else:
                P = np.concatenate(X,0)
            for (indices, scaler) in self.scalers:
                scaler.fit(P[:,indices])
        return self

    def transform(self, X):
        Xfit = [np.copy(d) for d in X]
        if self.use:
            for i in range(len(Xfit)):
                if Xfit[i].shape[0] > 0:
                    for (indices, scaler) in self.scalers:
                        Xfit[i][:,indices] = scaler.transform(Xfit[i][:,indices])
        return Xfit

class Padding(BaseEstimator, TransformerMixin):

    def __init__(self, use=False):
        self.use = use

    def fit(self, X, y=None):
        return self

    def transform(self, X):
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

    def __init__(self, use=False, num_pts=10, threshold=-1, location="upper", point_type="finite"):
        self.num_pts    = num_pts
        self.threshold  = threshold
        self.use        = use
        self.location   = location
        self.point_type = point_type

    def fit(self, X, y=None):
        return self

    def transform(self, X):
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

    def __init__(self, use=False, limit=np.inf, point_type="finite"):
        self.use, self.limit, self.point_type = use, limit, point_type

    def fit(self, X, y=None):
        return self

    def transform(self, X):
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
