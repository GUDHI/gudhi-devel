# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2026 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

__license__ = "MIT"

from ._default_random_seed_ext import _set_seed as _set_gudhi_cxx_seed


def set_seed(seed: int):
    """
    Sets the seed for internal GUDHI functionnalities that are using randomness.

    :param seed: The new seed value.
    
    .. note::

        This function also sets the default random CGAL seed (cf.
        `CGAL::Random <https://doc.cgal.org/latest/Generator/classCGAL_1_1Random.html>`_).

    """
    _set_gudhi_cxx_seed(seed)
    # Nothing to be done if no CGAL
    try:
        from ._default_cgal_random_seed_ext import _set_seed as _set_cgal_seed
        _set_cgal_seed(seed)
    except ImportError:
        pass
