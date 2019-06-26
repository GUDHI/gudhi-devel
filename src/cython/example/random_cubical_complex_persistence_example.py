#!/usr/bin/env python

import gudhi
import numpy
from functools import reduce
import argparse
import operator


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

parser = argparse.ArgumentParser(description='Random cubical complex.',
                                 epilog='Example: '
                                 './random_cubical_complex_persistence_example.py'
                                 ' 10 10 10 - Constructs a random cubical '
                                 'complex in a dimension [10, 10, 10] (aka. '
                                 '1000 random top dimensional cells).')
parser.add_argument('dimension', type=int, nargs="*",
                    help='Cubical complex dimensions')

args = parser.parse_args()
dimension_multiplication = reduce(operator.mul, args.dimension, 1)

if dimension_multiplication > 1: 
    print("#####################################################################")
    print("CubicalComplex creation")
    cubical_complex = gudhi.CubicalComplex(dimensions=args.dimension,
                                           top_dimensional_cells = numpy.random.rand(dimension_multiplication))

    print("persistence(homology_coeff_field=2, min_persistence=0)=")
    print(cubical_complex.persistence(homology_coeff_field=2, min_persistence=0))

    print("betti_numbers()=")
    print(cubical_complex.betti_numbers())
