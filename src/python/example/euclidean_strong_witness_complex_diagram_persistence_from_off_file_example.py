#!/usr/bin/env python

import argparse
import errno
import os
import matplotlib.pyplot as plot
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
    description="EuclideanStrongWitnessComplex creation from "
    "points read in a OFF file.",
    epilog="Example: "
    "example/euclidean_strong_witness_complex_diagram_persistence_from_off_file_example.py "
    "-f ../data/points/tore3D_300.off -a 1.0 -n 20 -d 2"
    "- Constructs a strong witness complex with the "
    "points from the given OFF file.",
)
parser.add_argument("-f", "--file", type=str, required=True)
parser.add_argument("-a", "--max_alpha_square", type=float, required=True)
parser.add_argument("-n", "--number_of_landmarks", type=int, required=True)
parser.add_argument("-d", "--limit_dimension", type=int, required=True)
parser.add_argument("-b", "--band", type=float, default=0.0)
parser.add_argument(
    "--no-diagram",
    default=False,
    action="store_true",
    help="Flag for not to display the diagrams",
)

args = parser.parse_args()

with open(args.file, "r") as f:
    first_line = f.readline()
    if (first_line == "OFF\n") or (first_line == "nOFF\n"):
        print("##############################################################")
        print("EuclideanStrongWitnessComplex creation from points read "\
            "in a OFF file")

        witnesses = gudhi.read_points_from_off_file(off_file=args.file)
        landmarks = gudhi.pick_n_random_points(
            points=witnesses, nb_points=args.number_of_landmarks
        )

        message = (
            "EuclideanStrongWitnessComplex with max_edge_length="
            + repr(args.max_alpha_square)
            + " - Number of landmarks="
            + repr(args.number_of_landmarks)
        )
        print(message)

        witness_complex = gudhi.EuclideanStrongWitnessComplex(
            witnesses=witnesses, landmarks=landmarks
        )
        simplex_tree = witness_complex.create_simplex_tree(
            max_alpha_square=args.max_alpha_square,
            limit_dimension=args.limit_dimension
        )

        message = "Number of simplices=" + repr(simplex_tree.num_simplices())
        print(message)

        diag = simplex_tree.persistence()

        print("betti_numbers()=")
        print(simplex_tree.betti_numbers())

        if args.no_diagram == False:
            gudhi.plot_persistence_diagram(diag, band=args.band)
            plot.show()
    else:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                                args.file)

    f.close()
