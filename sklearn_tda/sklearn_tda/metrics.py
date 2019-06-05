"""
@author: Mathieu Carriere
All rights reserved
"""

import sys
import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin

try:
    from gudhi import bottleneck_distance
    USE_GUDHI = True
except ImportError:
    USE_GUDHI = False
    print("Gudhi not found: BottleneckDistance will return null matrix")

#############################################
# Metrics ###################################
#############################################

def compute_bott_matrix(diags1, diags2, epsilon=1e-3):

    num_diag1 = len(diags1)

    if np.array_equal(np.concatenate(diags1,0), np.concatenate(diags2,0)):
        matrix = np.zeros((num_diag1, num_diag1))

        if USE_GUDHI:
            for i in range(num_diag1):
                #sys.stdout.write( str(i*1.0 / num_diag1) + "\r")
                for j in range(i+1, num_diag1):
                    matrix[i,j] = bottleneck_distance(diags1[i], diags1[j], delta)
                    matrix[j,i] = matrix[i,j]
        else:
            print("Gudhi required---returning null matrix")

    else:
        num_diag2 = len(diags2)
        matrix = np.zeros((num_diag1, num_diag2))

        if USE_GUDHI:
            for i in range(num_diag1):
                #sys.stdout.write( str(i*1.0 / num_diag1) + "\r")
                for j in range(num_diag2):
                    matrix[i,j] = bottleneck_distance(diags1[i], diags2[j], epsilon)
        else:
            print("Gudhi required---returning null matrix")

    return matrix

class BottleneckDistance(BaseEstimator, TransformerMixin):

    def __init__(self, epsilon=1e-3):
        self.epsilon = epsilon

    def fit(self, X, y=None):
        self.diagrams_ = X
        return self

    def transform(self, X):
        return compute_bott_matrix(X, self.diagrams_, self.epsilon)
