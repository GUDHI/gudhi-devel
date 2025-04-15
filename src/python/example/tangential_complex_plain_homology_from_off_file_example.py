#!/usr/bin/env python

""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2016 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""
__license__ = "GPL v3"  # Because of TangentialComplex


import argparse
import errno
import os
import gudhi as gd


parser = argparse.ArgumentParser(
    description="TangentialComplex creation from points read in a OFF file.",
    epilog="Example: "
    "example/tangential_complex_plain_homology_from_off_file_example.py "
    "-f ../data/points/tore3D_300.off -i 3"
    "- Constructs a tangential complex with the "
    "points from the given OFF file",
)
parser.add_argument("-f", "--file", type=str, required=True)
parser.add_argument("-i", "--intrisic_dim", type=int, required=True)
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
        print("TangentialComplex creation from points read in a OFF file")

        tc = gd.TangentialComplex(intrisic_dim=args.intrisic_dim, off_file=args.file)
        tc.compute_tangential_complex()
        st = tc.create_simplex_tree()

        print(f"Number of simplices={st.num_simplices()}")

        diag = st.persistence(persistence_dim_max=True)

        print(f"betti_numbers()={st.betti_numbers()}")

        if args.no_diagram == False:
            import matplotlib.pyplot as plot

            gd.plot_persistence_diagram(diag, band=args.band)
            plot.show()
    else:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.file)

    f.close()
