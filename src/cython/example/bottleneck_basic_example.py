#!/usr/bin/env python

import gudhi

"""This file is part of the Gudhi Library. The Gudhi library
   (Geometric Understanding in Higher Dimensions) is a generic C++
   library for computational topology.

   Author(s):       Francois Godi, Vincent Rouvreau

   Copyright (C) 2016 INRIA

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

__author__ = "Francois Godi, Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 INRIA"
__license__ = "GPL v3"

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

