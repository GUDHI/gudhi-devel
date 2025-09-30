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

from gudhi import _random_ext as t

import numpy as np

class GudhiBitGenerator(np.random.BitGenerator):
    def __init__(self, seed: Optional[int] = None):
        super().__init__(0)
        if seed is None:
            self.rng = t.RandomGenerator()
        else:
            self.rng = t.RandomGenerator(seed)
        self.rng.setup_bitgen(self.capsule)
