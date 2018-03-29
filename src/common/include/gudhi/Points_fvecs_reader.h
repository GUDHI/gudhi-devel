/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2015  INRIA Saclay (France)
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
#ifndef POINTS_FVECS_READER_H_
#define POINTS_FVECS_READER_H_

namespace Gudhi {

template<typename Kernel, typename OutputIteratorPoints>
bool load_points_from_fvecs_file(const std::string &filename, OutputIteratorPoints points, int only_the_first_n_points = -1)
{
  typedef typename Kernel::Point_d Point;

  std::ifstream in(filename, std::ios::binary);
  if (!in.is_open()) {
    std::cerr << "Could not open '" << filename << "'" << std::endl;
    return false;
  }

  Kernel k;
  unsigned long pt_dim = 0;

  in.read(reinterpret_cast<char*>(&pt_dim), 4);
  std::vector<float> current_pt;
  current_pt.reserve(pt_dim);
  for (int c = 0; !in.fail() && c != only_the_first_n_points; ++c) {

    for (int j = 0; j < pt_dim; ++j)
    {
      float coord = 0.f;
      in.read(reinterpret_cast<char*>(&coord), 4);
      current_pt.push_back(coord);
    }

    *points++ = Point(current_pt.begin(), current_pt.end());
    current_pt.clear();
    in.read(reinterpret_cast<char*>(&pt_dim), 4);
  }

#ifdef DEBUG_TRACES
  std::cerr << "'" << filename << "' loaded." << std::endl;
#endif

  return true;
}


}  // namespace Gudhi

#endif  // POINTS_FVECS_READER_H_
