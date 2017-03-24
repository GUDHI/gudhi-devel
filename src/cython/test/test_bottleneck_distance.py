import gudhi

"""This file is part of the Gudhi Library. The Gudhi library
   (Geometric Understanding in Higher Dimensions) is a generic C++
   library for computational topology.

   Author(s):       Vincent Rouvreau

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

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 INRIA"
__license__ = "GPL v3"


def test_basic_bottleneck():
    diag1 = [[2.7, 3.7],[9.6, 14.],[34.2, 34.974], [3.,float('Inf')]]
    diag2 = [[2.8, 4.45],[9.5, 14.1],[3.2,float('Inf')]]

    assert(gudhi.bottleneck_distance(diag1, diag2, 0.1) == 0.8081763781405569)
    assert(gudhi.bottleneck_distance(diag1, diag2) == 0.75)
