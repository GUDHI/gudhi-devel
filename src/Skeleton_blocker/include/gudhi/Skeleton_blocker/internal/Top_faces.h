/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014  INRIA Sophia Antipolis-Mediterranee (France)
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
#ifndef SKELETON_BLOCKER_INTERNAL_TOP_FACES_H_
#define SKELETON_BLOCKER_INTERNAL_TOP_FACES_H_

#include <list>
#include <vector>
#include <set>

namespace Gudhi {

namespace skbl {

template<typename SimplexHandle>
std::list<SimplexHandle> subfaces(SimplexHandle top_face) {
  std::list<SimplexHandle> res;
  if (top_face.dimension() == -1) return res;
  if (top_face.dimension() == 0) {
    res.push_back(top_face);
    return res;
  } else {
    auto first_vertex = top_face.first_vertex();
    top_face.remove_vertex(first_vertex);
    res = subfaces(top_face);
    std::list<SimplexHandle> copy = res;
    for (auto& simplex : copy) {
      simplex.add_vertex(first_vertex);
    }
    res.push_back(SimplexHandle(first_vertex));
    res.splice(res.end(), copy);
    return res;
  }
}

/**
 * add all faces of top_face in simplices_per_dimension
 */
template<typename SimplexHandle>
void register_faces(std::vector< std::set<SimplexHandle> >& simplices_per_dimension,
                    const SimplexHandle& top_face) {
  std::list<SimplexHandle> subfaces_list = subfaces(top_face);
  for (auto& simplex : subfaces_list) {
    simplices_per_dimension[simplex.dimension()].insert(simplex);
  }
}

}  // namespace skbl

}  // namespace Gudhi

#endif  // SKELETON_BLOCKER_INTERNAL_TOP_FACES_H_
