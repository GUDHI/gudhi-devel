#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plot
import gudhi

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


def is_file_perseus(file):
    num_lines = open(file).read().count("\n")
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


parser = argparse.ArgumentParser(
    description="Periodic cubical complex from a " "Perseus-style file name.",
    epilog="Example: "
    "./periodic_cubical_complex_barcode_persistence_from_perseus_file_example.py"
    " -f ../data/bitmap/CubicalTwoSphere.txt",
)

parser.add_argument("-f", "--file", type=str, required=True)
parser.add_argument(
    "--no-barcode",
    default=False,
    action="store_true",
    help="Flag for not to display the barcodes",
)

args = parser.parse_args()

if is_file_perseus(args.file):
    print("#####################################################################")
    print("PeriodicCubicalComplex creation")
    periodic_cubical_complex = gudhi.PeriodicCubicalComplex(perseus_file=args.file)

    print("persistence(homology_coeff_field=3, min_persistence=0)=")
    diag = periodic_cubical_complex.persistence(
        homology_coeff_field=3, min_persistence=0
    )
    print(diag)

    print("betti_numbers()=")
    print(periodic_cubical_complex.betti_numbers())
    if args.no_barcode == False:
        gudhi.plot_persistence_barcode(diag)
        plot.show()
else:
    print(args.file, "is not a valid perseus style file")
