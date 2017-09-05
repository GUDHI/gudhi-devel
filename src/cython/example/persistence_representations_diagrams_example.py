#!/usr/bin/env python

import gudhi
import argparse

"""This file is part of the Gudhi Library. The Gudhi library
   (Geometric Understanding in Higher Dimensions) is a generic C++
   library for computational topology.

   Author(s):       Pawel Dlotko

   Copyright (C) 2017 Swansea University

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

__author__ = "Pawel Dlotko"
__copyright__ = "Copyright (C) 2017 Swansea University"
__license__ = "GPL v3"

print("#####################################################################")
print("Persistence representations diagrams example")


parser = argparse.ArgumentParser(description='Statistics of persistence diagrams from file ',
                                 epilog='Example: '
                                 'example/persistence_representations_diagrams_example.py '
                                 '-f file_with_diagram -d 1')
parser.add_argument("-f", "--file", type=str, required=True)
parser.add_argument("-d", "--dimension", type=int, default=0)

args = parser.parse_args()

print "Here are the parameters of the program: ",args.file," , " ,args.dimension

p = gudhi.PersistenceIntervals(None,args.dimension,args.file);
min_max_ = p.get_x_range();
print "Birth-death range : ", min_max_ 

dominant_ten_intervals_length = p.length_of_dominant_intervals(10)
print "Length of ten dominant intervals : ",  dominant_ten_intervals_length

ten_dominant_intervals = p.dominant_intervals(10);
print "Here are the dominant intervals : " , ten_dominant_intervals

histogram = p.histogram_of_lengths(10);
print "Here is the histogram of barcode's length : ", histogram

cumulative_histogram = p.cumulative_histogram_of_lengths(10)
print "Cumulative histogram : " ,cumulative_histogram

char_funct_diag = p.characteristic_function_of_diagram(min_max_[0], min_max_[1],None)
print "Characteristic function of diagram : ",char_funct_diag

cumul_char_funct_diag = p.cumulative_characteristic_function_of_diagram(min_max_[0], min_max_[1],None)
print "Cumulative characteristic function of diagram : ",cumul_char_funct_diag

pbns = p.compute_persistent_betti_numbers()
print "Persistence Betti numbers ", pbns
