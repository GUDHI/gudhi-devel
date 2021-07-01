#!/usr/bin/env python

import argparse
import gudhi as gd

""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ -
    which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full
    license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2016 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "MIT"

parser = argparse.ArgumentParser(
    description="AlphaComplex creation from " "points read in a OFF file.",
    epilog="Example: "
    "example/alpha_complex_diagram_persistence_from_off_file_example.py "
    "-f ../data/points/tore3D_300.off"
    "- Constructs a alpha complex with the "
    "points from the given OFF file.",
)
parser.add_argument("-f", "--file", type=str, required=True)
parser.add_argument("-a", "--max_alpha_square", type=float, required=False)
parser.add_argument("-b", "--band", type=float, default=0.0)
parser.add_argument(
    "--no-diagram",
    default=False,
    action="store_true",
    help="Flag for not to display the diagrams",
)

args = parser.parse_args()

print("##############################################################")
print("AlphaComplex creation from points read in a OFF file")

points = gd.read_points_from_off_file(off_file = args.file)
alpha_complex = gd.AlphaComplex(points = points)
if args.max_alpha_square is not None:
    print("with max_edge_length=", args.max_alpha_square)
    simplex_tree = alpha_complex.create_simplex_tree(
        max_alpha_square=args.max_alpha_square
    )
else:
    simplex_tree = alpha_complex.create_simplex_tree()

print("Number of simplices=", simplex_tree.num_simplices())

diag = simplex_tree.persistence()
print("betti_numbers()=", simplex_tree.betti_numbers())
if args.no_diagram == False:
    import matplotlib.pyplot as plot
    gd.plot_persistence_diagram(diag, band=args.band)
    plot.show()
