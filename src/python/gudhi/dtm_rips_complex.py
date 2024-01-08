# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Yuichi Ike, Raphaël Tinarrage
#
# Copyright (C) 2020 Inria, Copyright (C) 2020 FUjitsu Laboratories Ltd.
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification


from gudhi.weighted_rips_complex import WeightedRipsComplex
from gudhi.point_cloud.dtm import DistanceToMeasure
from scipy.spatial.distance import cdist
import numpy as np


class DTMRipsComplex(WeightedRipsComplex):
    """
    Class to generate a DTM Rips complex from a distance matrix or a point set,
    in the way described in :cite:`dtmfiltrations`.
    Remark that all the filtration values are doubled compared to the definition in the paper
    for the consistency with RipsComplex.
    
    :Requires: `SciPy <installation.html#scipy>`_
    """

    def __init__(self, points=None, distance_matrix=None, k=1, q=2, max_filtration=float("inf")):
        """
        Args:
            points (numpy.ndarray): array of points.
            distance_matrix (numpy.ndarray): full distance matrix.
            k (int): number of neighbors for the computation of DTM. Defaults to 1, which is equivalent to the usual
                Rips complex.
            q (float): order used to compute the distance to measure. Defaults to 2.
            max_filtration (float): specifies the maximal filtration value to be considered.
        """
        if distance_matrix is None:
            if points is not None:
                distance_matrix = cdist(points, points)
            else:
                # Empty Rips construction
                distance_matrix = np.ndarray((0,0))
        self.distance_matrix = distance_matrix

        # TODO: address the error when k is too large
        if k <= 1:
            self.weights = [0] * len(distance_matrix)
        else:
            dtm = DistanceToMeasure(k, q=q, metric="precomputed")
            self.weights = dtm.fit_transform(distance_matrix)
        self.max_filtration = max_filtration
