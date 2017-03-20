#!/usr/bin/env python

import gudhi
import argparse
import math

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

parser = argparse.ArgumentParser(description='AlphaComplex and RipsComplex '
                                 'persistence creation from points read in '
                                 'a OFF file. Bottleneck distance computation'
                                 ' on each dimension',
                                 epilog='Example: '
                                 'example/alpha_rips_persistence_bottleneck_distance.py '
                                 '-f ../data/points/tore3D_1307.off -t 0.15 -d 3')
parser.add_argument("-f", "--file", type=str, required=True)
parser.add_argument("-t", "--threshold", type=float, default=0.5)
parser.add_argument("-d", "--max_dimension", type=int, default=1)

args = parser.parse_args()
with open(args.file, 'r') as f:
    first_line = f.readline()
    if (first_line == 'OFF\n') or (first_line == 'nOFF\n'):
        print("#####################################################################")
        print("RipsComplex creation from points read in a OFF file")

        message = "RipsComplex with max_edge_length=" + repr(args.threshold)
        print(message)

        rips_complex = gudhi.RipsComplex(off_file=args.file,
                                         max_edge_length=args.threshold)

        rips_stree = rips_complex.create_simplex_tree(max_dimension=args.max_dimension)

        message = "Number of simplices=" + repr(rips_stree.num_simplices())
        print(message)

        rips_diag = rips_stree.persistence()

        print("#####################################################################")
        print("AlphaComplex creation from points read in a OFF file")

        message = "AlphaComplex with max_edge_length=" + repr(args.threshold)
        print(message)

        alpha_complex = gudhi.AlphaComplex(off_file=args.file)
        alpha_stree = alpha_complex.create_simplex_tree(max_alpha_square=(args.threshold * args.threshold))

        message = "Number of simplices=" + repr(alpha_stree.num_simplices())
        print(message)

        alpha_diag = alpha_stree.persistence()

        max_b_distance = 0.0
        for dim in range(args.max_dimension):
            # Alpha persistence values needs to be transform because filtration
            # values are alpha square values
            funcs = [math.sqrt, math.sqrt]
            alpha_intervals = []
            for interval in alpha_stree.persistence_intervals_in_dimension(dim):
                alpha_intervals.append(map(lambda func,value: func(value), funcs, interval))

            rips_intervals = rips_stree.persistence_intervals_in_dimension(dim)
            bottleneck_distance = gudhi.bottleneck_distance(rips_intervals, alpha_intervals)
            message = "In dimension " + repr(dim) + ", bottleneck distance = " + repr(bottleneck_distance)
            print(message)
            max_b_distance = max(bottleneck_distance, max_b_distance)

        print("================================================================================")
        message = "Bottleneck distance is " + repr(max_b_distance)
        print(message)

    else:
        print(args.file, "is not a valid OFF file")

    f.close()
