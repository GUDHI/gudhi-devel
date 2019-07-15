/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2019 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PERMUTAHEDRAL_REPRESENTATION_JOIN_H_
#define PERMUTAHEDRAL_REPRESENTATION_JOIN_H_

#include <set>

namespace Gudhi {

namespace coxeter_triangulation {

 /** \brief Computes the permutahedral representation of the join 
  *  (the minimal common coface by inclusion) if exists, returns an
  *  empty permutahedral representation otherwise. 
  *
  * \tparam Permutahedral_representation_range A range of elements
  *  of type Permutahedral_representation.
  *
  * @param[in] simplex_range Range of input simplices.
  */
template <class Permutahedral_representation_range>
typename Permutahedral_representation_range::value_type
join(const Permutahedral_representation_range& simplex_range) {
  using Permutahedral_representation = typename Permutahedral_representation_range::value_type;
  using Vertex = typename Permutahedral_representation::Vertex;
  using Ordered_partition = typename Permutahedral_representation::OrderedSetPartition;
  using Part = typename Ordered_partition::value_type;
  if (simplex_range.empty())
    return Permutahedral_representation();
  std::set<Vertex> vertices;
  for (auto s: simplex_range)
    for (auto v: s.vertex_range())
      vertices.insert(v);
  Permutahedral_representation result;
  result.vertex() = *vertices.begin();
  std::size_t d = vertices.begin()->size();
  auto curr_it = vertices.begin();
  auto prev_it = curr_it++;
  std::set<std::size_t> indices;
  for (std::size_t k = 0; k < d+1; ++k)
    indices.insert(k);
  while (curr_it != vertices.end()) {
    Part part;
    result.partition().emplace_back(Part());
    for (std::size_t i = 0; i < d; ++i)
      if (curr_it->at(i) < prev_it->at(i))
	return Permutahedral_representation();
      else if (curr_it->at(i) == prev_it->at(i) + 1) {
	result.partition().back().push_back(i);
	indices.erase(i);
      }
      else if (curr_it->at(i) > prev_it->at(i) + 1)
	return Permutahedral_representation();
    curr_it++;
    prev_it++;
  }
  result.partition().emplace_back(Part());
  for (std::size_t k: indices)
    result.partition().back().push_back(k);
  return result;
}

} // namespace coxeter_triangulation

} // namespace Gudhi


#endif
