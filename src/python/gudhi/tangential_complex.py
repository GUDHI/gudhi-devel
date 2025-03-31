# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2025 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

import os

from gudhi import _tangential_complex_ext as t

from gudhi.simplex_tree import SimplexTree

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2025 Inria"
__license__ = "GPL v3"

# TangentialComplex python interface
class TangentialComplex(t._Tangential_complex_interface):
    """The class Tangential_complex represents a tangential complex. After the
    computation of the complex, an optional post-processing called perturbation
    can be run to attempt to remove inconsistencies.
    """

    def __init__(self, intrisic_dim, points=None, off_file=''):
        """TangentialComplex constructor.

        :param intrisic_dim: Intrinsic dimension of the manifold.
        :type intrisic_dim: integer

        :param points: A list of points in d-Dimension.
        :type points: list of list of double

        Or

        :param off_file: An OFF file style name.
        :type off_file: string
        """
        if off_file:
            if os.path.isfile(off_file):
                super().__init__(intrisic_dim, off_file, True)
            else:
                print("file " + off_file + " not found.")
        else:
            if points is None:
                # Empty tangential construction
                points=[]
            super().__init__(intrisic_dim, points)


    def create_simplex_tree(self):
        """
        Exports the complex into a simplex tree.
        Returns:
            SimplexTree: A simplex tree created from the complex
        """
        stree = SimplexTree()
        super().create_simplex_tree(stree)
        return stree
