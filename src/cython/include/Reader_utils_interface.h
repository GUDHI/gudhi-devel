/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2017 INRIA
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef INCLUDE_READER_UTILS_INTERFACE_H_
#define INCLUDE_READER_UTILS_INTERFACE_H_

#include <gudhi/reader_utils.h>

#include <iostream>
#include <vector>
#include <string>

namespace Gudhi {

// Redefine functions with a different name in order the original name can be used in the Python version.
std::vector<std::vector<double>> read_matrix_from_csv_file(const std::string& filename,
                                                                const char separator = ';') {
  return read_lower_triangular_matrix_from_csv_file<double>(filename, separator);
}

inline std::map<int, std::vector<std::pair<double, double>>>
    read_pers_intervals_grouped_by_dimension(std::string const& filename) {
  return read_persistence_intervals_grouped_by_dimension(filename);
}

inline std::vector<std::pair<double, double>>
    read_pers_intervals_in_dimension(std::string const& filename, int only_this_dim = -1) {
  return read_persistence_intervals_in_dimension(filename, only_this_dim);
}


}  // namespace Gudhi


#endif  // INCLUDE_READER_UTILS_INTERFACE_H_
