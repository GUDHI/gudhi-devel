/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef SKELETON_BLOCKER_INTERNAL_TOP_FACES_H_
#define SKELETON_BLOCKER_INTERNAL_TOP_FACES_H_

#include <list>
#include <vector>
#include <set>

namespace Gudhi {

namespace skeleton_blocker {

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

}  // namespace skeleton_blocker

namespace skbl = skeleton_blocker;

}  // namespace Gudhi

#endif   // SKELETON_BLOCKER_INTERNAL_TOP_FACES_H_
