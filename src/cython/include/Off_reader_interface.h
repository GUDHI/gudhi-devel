/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016 INRIA
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

#ifndef OFF_READER_INTERFACE_H
#define	OFF_READER_INTERFACE_H

#include <gudhi/Points_off_io.h>

#include <iostream>
#include <vector>
#include <string>

namespace Gudhi {

std::vector<std::vector<double>> read_points_from_OFF_file(const std::string& off_file) {
  Gudhi::Points_off_reader<std::vector<double>> off_reader(off_file);
  return off_reader.get_point_cloud();
}

}  // namespace Gudhi

#endif  // OFF_READER_INTERFACE_H

