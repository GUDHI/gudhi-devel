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
print("AlphaComplex creation from points read in a file")

parser = argparse.ArgumentParser(description='AlphaComplex creation from '
                                 'points read in a file.',
                                 epilog='Example: '
                                 'example/alpha_complex_from_file_example.py '
                                 'data/500_random_points_on_3D_Torus.csv '
                                 '- Constructs a alpha complex with the '
                                 'points from the given file. File format '
                                 'is X1, X2, ..., Xn')
parser.add_argument('file', type=argparse.FileType('r'))
args = parser.parse_args()

points = pandas.read_csv(args.file, header=None)

alpha_complex = gudhi.AlphaComplex(points=points.values,
                                   max_alpha_square=0.5)

print("dimension=", alpha_complex.dimension())
print("point[0]=", alpha_complex.get_point(0))
print("point[5]=", alpha_complex.get_point(5))

alpha_complex.initialize_filtration()
diag = alpha_complex.persistence(homology_coeff_field=2, min_persistence=0.1)

print("betti_numbers()=")
print(alpha_complex.betti_numbers())

print("star([0])=", alpha_complex.get_star_tree([0]))
print("coface([0], 1)=", alpha_complex.get_coface_tree([0], 1))
