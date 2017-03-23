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

parser = argparse.ArgumentParser(description='AlphaComplex creation from '
                                 'points read in a OFF file.',
                                 epilog='Example: '
                                 'example/alpha_complex_diagram_persistence_from_off_file_example.py '
                                 '-f ../data/points/tore3D_300.off -a 0.6'
                                 '- Constructs a alpha complex with the '
                                 'points from the given OFF file.')
parser.add_argument("-f", "--file", type=str, required=True)
parser.add_argument("-a", "--max_alpha_square", type=float, default=0.5)
parser.add_argument('--no-diagram', default=False, action='store_true' , help='Flag for not to display the diagrams')

args = parser.parse_args()

with open(args.file, 'r') as f:
    first_line = f.readline()
    if (first_line == 'OFF\n') or (first_line == 'nOFF\n'):
        print("#####################################################################")
        print("AlphaComplex creation from points read in a OFF file")
        
        message = "AlphaComplex with max_edge_length=" + repr(args.max_alpha_square)
        print(message)
        
        alpha_complex = gudhi.AlphaComplex(off_file=args.file)
        simplex_tree = alpha_complex.create_simplex_tree(max_alpha_square=args.max_alpha_square)
    
        message = "Number of simplices=" + repr(simplex_tree.num_simplices())
        print(message)
        
        diag = simplex_tree.persistence()
    
        print("betti_numbers()=")
        print(simplex_tree.betti_numbers())
    
        if args.no_diagram == False:
            gudhi.plot_persistence_diagram(diag)
    else:
        print(args.file, "is not a valid OFF file")

    f.close()
