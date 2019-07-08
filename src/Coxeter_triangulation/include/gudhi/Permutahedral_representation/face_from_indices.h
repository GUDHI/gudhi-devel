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

#ifndef PERMUTAHEDRAL_REPRESENTATION_FACE_FROM_INDICES_H_
#define PERMUTAHEDRAL_REPRESENTATION_FACE_FROM_INDICES_H_

namespace Gudhi {

namespace coxeter_triangulation {

 /** \brief Computes the permutahedral representation of a face of a given simplex
  *  and a range of the vertex indices that compose the face.
  *
  * \tparam Permutahedral_representation has to be Permutahedral_representation
  * \tparam Index_range is a range of unsigned integers taking values in 0,...,k,
  * where k is the dimension of the simplex simplex.
  *
  * @param[in] simplex Input simplex.
  * @param[in] indices Input range of indices.
  */
template <class Permutahedral_representation,
	  class Index_range>
Permutahedral_representation face_from_indices(const Permutahedral_representation& simplex,
					       const Index_range& indices) {
  using range_index = typename Index_range::value_type;
  using Ordered_set_partition = typename Permutahedral_representation::OrderedSetPartition;
  using Part = typename Ordered_set_partition::value_type;
  using part_index = typename Part::value_type;
  Permutahedral_representation value;
  std::size_t d = simplex.vertex().size();
  value.vertex() = simplex.vertex();
  std::size_t k = indices.size()-1;
  value.partition().resize(k+1);
  std::size_t l = simplex.partition().size()-1;
  for (std::size_t h = 1; h < k+1; h++)
    for (range_index i = indices[h-1]; i < indices[h]; i++)
      for (part_index j: simplex.partition()[i])
	value.partition()[h-1].push_back(j);
  for (range_index i = indices[k]; i < l+1; i++)
    for (part_index j: simplex.partition()[i])
      value.partition()[k].push_back(j);
  for (range_index i = 0; i < indices[0]; i++)
    for (part_index j: simplex.partition()[i]) {
      if (j != d)
	value.vertex()[j]++;
      else
	for (std::size_t l = 0; l < d; l++)
	  value.vertex()[l]--;
      value.partition()[k].push_back(j);
    }
  // sort the values in each part (probably not needed)
  for (auto& part: value.partition())
    std::sort(part.begin(), part.end());
  return value;
}

} // namespace coxeter_triangulation

} // namespace Gudhi


#endif
