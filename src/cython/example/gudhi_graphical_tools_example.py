#!/usr/bin/env python

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

print("#####################################################################")
print("Show palette colors values for dimension")

gudhi.show_palette_values()

print("#####################################################################")
print("Show barcode persistence example")

persistence = [(2, (1.0, float('inf'))), (1, (1.4142135623730951, float('inf'))),
               (1, (1.4142135623730951, float('inf'))), (0, (0.0, float('inf'))),
               (0, (0.0, 1.0)), (0, (0.0, 1.0)), (0, (0.0, 1.0))]
gudhi.barcode_persistence(persistence)

print("#####################################################################")
print("Show diagram persistence example")

gudhi.diagram_persistence(persistence)
