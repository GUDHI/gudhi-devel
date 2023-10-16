#!/usr/bin/env python

import argparse
import errno
import os
import gudhi

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
    description="RipsComplex creation from " "points read in a OFF file.",
    epilog="Example: "
    "example/rips_complex_diagram_persistence_from_off_file_example.py "
    "-f ../data/points/tore3D_300.off -a 0.6"
    "- Constructs a Rips complex with the "
    "points from the given OFF file.",
)
parser.add_argument("-f", "--file", type=str, required=True)
parser.add_argument("-e", "--max_edge_length", type=float, default=0.5)
parser.add_argument("-d", "--max_dimension", type=int, default=1)
parser.add_argument("-b", "--band", type=float, default=0.0)
parser.add_argument(
    "--no-diagram",
    default=False,
    action="store_true",
    help="Flag for not to display the diagrams",
)

args = parser.parse_args()

with open(args.file) as f:
    first_line = f.readline()
    if (first_line == "OFF\n") or (first_line == "nOFF\n"):
        print("##############################################################")
        print("RipsComplex creation from points read in a OFF file")

        message = "RipsComplex with max_edge_length=" + \
            repr(args.max_edge_length)
        print(message)

        point_cloud = gudhi.read_points_from_off_file(off_file=args.file)
        rips_complex = gudhi.RipsComplex(
            points=point_cloud, max_edge_length=args.max_edge_length
        )
        simplex_tree = rips_complex.create_simplex_tree(
            max_dimension=args.max_dimension
        )

        message = "Number of simplices=" + repr(simplex_tree.num_simplices())
        print(message)

        diag = simplex_tree.persistence()

        print("betti_numbers()=")
        print(simplex_tree.betti_numbers())

        if args.no_diagram == False:
            import matplotlib.pyplot as plot
            gudhi.plot_persistence_diagram(diag, band=args.band)
            plot.show()
    else:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                                args.file)

    f.close()
