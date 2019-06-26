#!/usr/bin/env python

import gudhi

""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2016 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

__author__ = "Francois Godi, Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "MIT"

diag1 = [[2.7, 3.7],[9.6, 14.],[34.2, 34.974], [3.,float('Inf')]]

diag2 = [[2.8, 4.45],[9.5, 14.1],[3.2,float('Inf')]]

message = "diag1=" + repr(diag1)
print(message)

message = "diag2=" + repr(diag2)
print(message)

message = "Bottleneck distance approximation=" + repr(gudhi.bottleneck_distance(diag1, diag2, 0.1))
print(message)

message = "Bottleneck distance exact value=" + repr(gudhi.bottleneck_distance(diag1, diag2))
print(message)

