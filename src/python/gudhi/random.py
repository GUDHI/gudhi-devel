# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):  Vincent Rouvreau
#
# Copyright (C) 2025 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

# Based on https://doc.cgal.org/latest/Generator/index.html
__license__ = "LGPL v3"

from typing import Optional

from gudhi._random_ext import GudhiRandomGenerator

import numpy as np


class GudhiBitGenerator(np.random.BitGenerator):
    """GUDHI random generator that fits numpy random generator requirements."""
    def __init__(self, seed: Optional[int] = None):
        """GudhiBitGenerator constructor.

        :param seed: A seed to initialize the BitGenerator. If None, then it will be set with a fresh and random seed.
            Default is None.
        :type seed: Optional[int]
        """
        super().__init__(0)
        if seed is None:
            self.rng = GudhiRandomGenerator()
        else:
            self.rng = GudhiRandomGenerator(seed)
        self.rng.setup_bitgen(self.capsule)
