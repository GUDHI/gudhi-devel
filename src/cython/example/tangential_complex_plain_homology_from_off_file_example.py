#!/usr/bin/env python

import gudhi
import argparse

""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2016 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "MIT"

parser = argparse.ArgumentParser(description='TangentialComplex creation from '
                                 'points read in a OFF file.',
                                 epilog='Example: '
                                 'example/tangential_complex_plain_homology_from_off_file_example.py '
                                 '-f ../data/points/tore3D_300.off -i 3'
                                 '- Constructs a tangential complex with the '
                                 'points from the given OFF file')
parser.add_argument("-f", "--file", type=str, required=True)
parser.add_argument("-i", "--intrisic_dim", type=int, required=True)
parser.add_argument("-b", "--band", type=float, default=0.)
parser.add_argument('--no-diagram', default=False, action='store_true' , help='Flag for not to display the diagrams')

args = parser.parse_args()

with open(args.file, 'r') as f:
    first_line = f.readline()
    if (first_line == 'OFF\n') or (first_line == 'nOFF\n'):
        print("#####################################################################")
        print("TangentialComplex creation from points read in a OFF file")
        
        tc = gudhi.TangentialComplex(intrisic_dim = args.intrisic_dim, off_file=args.file)
        tc.compute_tangential_complex()
        st = tc.create_simplex_tree()
    
        message = "Number of simplices=" + repr(st.num_simplices())
        print(message)
        
        diag = st.persistence(persistence_dim_max = True)
    
        print("betti_numbers()=")
        print(st.betti_numbers())
    
        if args.no_diagram == False:
            pplot = gudhi.plot_persistence_diagram(diag, band=args.band)
            pplot.show()
    else:
        print(args.file, "is not a valid OFF file")

    f.close()
