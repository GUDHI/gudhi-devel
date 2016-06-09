#!/usr/bin/env python

import gudhi
import numpy
import argparse
import operator


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

parser = argparse.ArgumentParser(description='Cubical complex from a perseus'
                                 'file style name.',
                                 epilog='Example: '
                                 './cubical_complex_from_perseus_file_example.py'
                                 '../data/bitmap/CubicalTwoSphere.txt')
parser.add_argument('perseus_file', type=str, required=True,
                    help='Perseus file style name')

args = parser.parse_args()

print("#####################################################################")
print("CubicalComplex creation")
cubical_complex = gudhi.CubicalComplex(perseus_file=args.perseus_file)

print("persistence(homology_coeff_field=2, min_persistence=0)=")
print(cubical_complex.persistence(homology_coeff_field=2, min_persistence=0))

print("betti_numbers()=")
print(cubical_complex.betti_numbers())