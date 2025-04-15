#!/usr/bin/env python

""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Francois Godi, Vincent Rouvreau

    Copyright (C) 2016 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

# Because of bottleneck_distance
# Would be "BSD-3-Clause" if using "from gudhi.hera import bottleneck_distance"
__license__ = "GPL v3"


import numpy as np
from gudhi import bottleneck_distance


diag1 = np.array([[2.7, 3.7], [9.6, 14.0], [34.2, 34.974], [3.0, float("Inf")]])

diag2 = np.array([[2.8, 4.45], [9.5, 14.1], [3.2, float("Inf")]])

print(f"diag1={diag1}")

print(f"diag2={diag2}")

print(f"Bottleneck distance approximation={bottleneck_distance(diag1, diag2, 0.1)}")

print(f"Bottleneck distance exact value={bottleneck_distance(diag1, diag2)}")
