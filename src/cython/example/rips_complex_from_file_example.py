#!/usr/bin/env python

import gudhi
import pandas
import argparse

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
print("RipsComplex creation from points read in a file")

parser = argparse.ArgumentParser(description='RipsComplex creation from '
                                 'points read in a file.',
                                 epilog='Example: '
                                 'example/rips_complex_from_file_example.py '
                                 'data/2000_random_points_on_3D_Torus.csv '
                                 '- Constructs a rips complex with the '
                                 'points from the given file. File format '
                                 'is X1, X2, ..., Xn')
parser.add_argument('file', type=argparse.FileType('r'))
args = parser.parse_args()

points = pandas.read_csv(args.file, header=None)

print("RipsComplex with max_edge_length=0.7")

rips_complex = gudhi.RipsComplex(points=points.values,
                                 max_dimension=len(points.values[0]), max_edge_length=0.7)

rips_complex.initialize_filtration()
diag = rips_complex.persistence(homology_coeff_field=2, min_persistence=0.3)

print("betti_numbers()=")
print(rips_complex.betti_numbers())

gudhi.diagram_persistence(diag)

gudhi.barcode_persistence(diag)
