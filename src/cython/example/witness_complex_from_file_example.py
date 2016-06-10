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

print("WitnessComplex with number_of_landmarks=100 alpha=0.7 epsilon_mu=0.001 max_dim=10")

witness_complex = gudhi.WitnessComplex(points=points.values,
                                       number_of_landmarks=100,
                                       max_alpha_square=0.7,
                                       mu_epsilon=0.001,
                                       dimension_limit=10)

witness_complex.initialize_filtration()
diag = witness_complex.persistence(homology_coeff_field=2, min_persistence=0.1)

print("betti_numbers()=")
print(witness_complex.betti_numbers())

gudhi.diagram_persistence(diag)

gudhi.barcode_persistence(diag)
