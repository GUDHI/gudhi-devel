# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Marc Glisse
#
# Copyright (C) 2023 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

from sklearn.base import BaseEstimator, TransformerMixin

import numpy as np


class Reshape(BaseEstimator, TransformerMixin):
    """
    This class allows reshaping arrays in a Scikit-Learn pipeline.
    """

    def __init__(self, newshape=(-1,), order="C"):
        """
        The arguments are the same as `numpy.reshape`.
        """
        self.newshape = newshape
        self.order = order

    def fit(self, X, Y=None):
        """
        Nothing to be done, but useful when included in a scikit-learn Pipeline.
        """
        return self

    def transform(self, X, Y=None):
        """Replace each x in X with `numpy.reshape(x, newshape, order)`.

        :param X: List of arrays to be reshaped
        :type X: list of array_like

        :return: List of reshaped arrays
        :rtype: list of numpy.ndarray
        """
        return [np.reshape(x, self.newshape, self.order) for x in X]
