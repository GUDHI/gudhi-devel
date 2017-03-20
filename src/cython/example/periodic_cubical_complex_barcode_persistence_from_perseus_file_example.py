#!/usr/bin/env python

import gudhi
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

def is_file_perseus(file):
    num_lines = open(file).read().count('\n')
    try:
        f = open(file)
        num_dim = int(f.readline())
        coeff = 1
        for dim in range(0, num_dim):
            try:
                line = int(f.readline())
                coeff *= abs(line)
            except ValueError:
                return False
        if num_lines == (1 + num_dim + coeff):
            return True
        else:
            return False
    except ValueError:
        return False

parser = argparse.ArgumentParser(description='Periodic cubical complex from a '
                                 'perseus file style name.',
                                 epilog='Example: '
                                 './periodic_cubical_complex_barcode_persistence_from_perseus_file_example.py'
                                 ' -f ../data/bitmap/CubicalTwoSphere.txt')

parser.add_argument("-f", "--file", type=str, required=True)
parser.add_argument('--no-barcode', default=False, action='store_true' , help='Flag for not to display the barcodes')

args = parser.parse_args()

if is_file_perseus(args.file):
    print("#####################################################################")
    print("PeriodicCubicalComplex creation")
    periodic_cubical_complex = gudhi.PeriodicCubicalComplex(perseus_file=args.file)

    print("persistence(homology_coeff_field=3, min_persistence=0)=")
    diag = periodic_cubical_complex.persistence(homology_coeff_field=3, min_persistence=0)
    print(diag)

    print("betti_numbers()=")
    print(periodic_cubical_complex.betti_numbers())
    if args.no_barcode == False:
        gudhi.plot_barcode_persistence(diag)
else:
    print(args.file, "is not a valid perseus style file")
