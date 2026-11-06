# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2016 Inria
#
# Modification(s):
#   - 2025/03 Thibaud Kloczko: Use nanobind instead of Cython for python bindings.
#   - 2025/04 Hannah Schreiber: Re-add possibility of tensors (numpy, torch etc.) as input.
#   - YYYY/MM Author: Description of the modification

__license__ = "GPL v3"

import os
import warnings

from gudhi import _tangential_complex_ext as t
from gudhi.simplex_tree import SimplexTree


# TangentialComplex python interface
class TangentialComplex(t._Tangential_complex_interface):
    """The class Tangential_complex represents a tangential complex. After the
    computation of the complex, an optional post-processing called perturbation
    can be run to attempt to remove inconsistencies.
    """

    def __init__(self, intrinsic_dim, intrisic_dim=None, points=None, off_file=""):
        """TangentialComplex constructor.

        :param intrinsic_dim: Intrinsic dimension of the manifold.
        :type intrinsic_dim: integer

        :param intrisic_dim: **[deprecated]** consider using `intrinsic_dim` instead.

        :param points: A list of points in d-Dimension.
        :type points (Sequence[Sequence[float]]): list of list of double

        Or

        :param off_file: An OFF file style name.
        :type off_file: string
        """
        if intrisic_dim is not None:
            warnings.warn(
                "The 'intrisic_dim' argument is deprecated and will be removed in a future version. Use 'intrinsic_dim' instead.",
                DeprecationWarning
            )
            intrinsic_dim = intrisic_dim
        if off_file:
            if os.path.isfile(off_file):
                super().__init__(intrinsic_dim, off_file)
            else:
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), off_file)
        else:
            if points is None:
                super().__init__(intrinsic_dim)
            else:
                super().__init__(intrinsic_dim, points)

    def create_simplex_tree(self):
        """Exports the complex into a simplex tree.

        :returns: A simplex tree created from the complex.
        :rtype: SimplexTree
        """
        stree = SimplexTree()
        super().create_simplex_tree(stree)
        return stree
