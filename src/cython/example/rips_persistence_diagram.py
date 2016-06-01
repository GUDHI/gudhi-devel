#!/usr/bin/env python

import gudhi
import matplotlib.pyplot as plt

"""This file is part of the Gudhi Library. The Gudhi library
   (Geometric Understanding in Higher Dimensions) is a generic C++
   library for computational topology.

   Author(s):       Marc Glisse

   Copyright (C) 2016  INRIA Saclay (France)

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

__author__ = "Marc Glisse"
__copyright__ = "Copyright (C) 2016  INRIA Saclay (France)"
__license__ = "GPL v3"

print("#####################################################################")
print("RipsComplex creation from points")
rips = gudhi.RipsComplex(points=[[0, 0], [1, 0], [0, 1], [1, 1]],
                         max_dimension=1, max_edge_length=42)

diag = rips.persistence(homology_coeff_field=2, min_persistence=0)
print("diag=", diag)

diag0=[i[1] for i in diag if i[0]==0]
diag1=[i[1] for i in diag if i[0]==1]

plt.plot([i[0] for i in diag0],[min(i[1],1.5) for i in diag0],'ro',
         [i[0] for i in diag1],[i[1] for i in diag1],'bs',[-.2,1.5],[-.2,1.5],
         '-')
plt.axis([-.2,1.5,-.2,1.5])
plt.show()
