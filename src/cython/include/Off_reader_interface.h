/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef INCLUDE_OFF_READER_INTERFACE_H_
#define INCLUDE_OFF_READER_INTERFACE_H_

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

#endif  // INCLUDE_OFF_READER_INTERFACE_H_

