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
print("WitnessComplex creation from points read in a file")

parser = argparse.ArgumentParser(description='WitnessComplex creation from '
                                 'points read in a file.',
                                 epilog='Example: '
                                 'example/witness_complex_from_file_example.py'
                                 ' data/2000_random_points_on_3D_Torus.csv '
                                 '- Constructs a witness complex with the '
                                 'points from the given file. File format '
                                 'is X1, X2, ..., Xn')
parser.add_argument('file', type=argparse.FileType('r'))
args = parser.parse_args()

points = pandas.read_csv(args.file, header=None)

print("WitnessComplex with number_of_landmarks=100")

witness_complex = gudhi.WitnessComplex(points=points.values,
                                       number_of_landmarks=100)

witness_complex.initialize_filtration()

print("filtered_tree=", witness_complex.get_filtered_tree())
