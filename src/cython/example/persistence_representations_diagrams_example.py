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











parser = argparse.ArgumentParser(description='Statistics od persistence diagrams from file ',
                                 epilog='Example: '
                                 'example/persistence_representations_diagrams_example.py '
                                 '-f file_with_diagram -d 1')
parser.add_argument("-f", "--file", type=str, required=True)
parser.add_argument("-d", "--dimension", type=int, default=0)

args = parser.parse_args()

print "Here are the parameters of the program: ",args.file," , " ,args.dimension


p = gudhi.PersistenceIntervals(args.file,args.dimension);
min_max_ = p.get_x_range();
print( "Birth-death range : ", min_max_)
"""
std::vector<double> dominant_ten_intervals_length = p.length_of_dominant_intervals(10);
std::cout << "Length of ten dominant intervals : " << std::endl;
for (size_t i = 0; i != dominant_ten_intervals_length.size(); ++i) {
std::cout << dominant_ten_intervals_length[i] << std::endl;
}

std::vector<std::pair<double, double> > ten_dominant_intervals = p.dominant_intervals(10);
std::cout << "Here are the dominant intervals : " << std::endl;
for (size_t i = 0; i != ten_dominant_intervals.size(); ++i) {
std::cout << "( " << ten_dominant_intervals[i].first << "," << ten_dominant_intervals[i].second << std::endl;
}

std::vector<size_t> histogram = p.histogram_of_lengths(10);
std::cout << "Here is the histogram of barcode's length : " << std::endl;
for (size_t i = 0; i != histogram.size(); ++i) {
std::cout << histogram[i] << " ";
}
std::cout << std::endl;

std::vector<size_t> cumulative_histogram = p.cumulative_histogram_of_lengths(10);
std::cout << "Cumulative histogram : " << std::endl;
for (size_t i = 0; i != cumulative_histogram.size(); ++i) {
std::cout << cumulative_histogram[i] << " ";
}
std::cout << std::endl;

std::vector<double> char_funct_diag = p.characteristic_function_of_diagram(min_max_.first, min_max_.second);
std::cout << "Characteristic function of diagram : " << std::endl;
for (size_t i = 0; i != char_funct_diag.size(); ++i) {
std::cout << char_funct_diag[i] << " ";
}
std::cout << std::endl;

std::vector<double> cumul_char_funct_diag =
  p.cumulative_characteristic_function_of_diagram(min_max_.first, min_max_.second);
std::cout << "Cumulative characteristic function of diagram : " << std::endl;
for (size_t i = 0; i != cumul_char_funct_diag.size(); ++i) {
std::cout << cumul_char_funct_diag[i] << " ";
}
std::cout << std::endl;

std::cout << "Persistence Betti numbers \n";
std::vector<std::pair<double, size_t> > pbns = p.compute_persistent_betti_numbers();
for (size_t i = 0; i != pbns.size(); ++i) {
std::cout << pbns[i].first << " " << pbns[i].second << std::endl;
}
"""
