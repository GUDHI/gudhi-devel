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

#ifndef PERMUTAHEDRAL_REPRESENTATION_SIMPLEX_COMPARATOR_H_
#define PERMUTAHEDRAL_REPRESENTATION_SIMPLEX_COMPARATOR_H_

namespace Gudhi {

namespace coxeter_triangulation {

/** \class Simplex_comparator 
 *  \brief A comparator class for Permutahedral_representation.
 *   The comparison is in lexicographic order first on
 *   vertices and then on ordered partitions with sorted parts.
 *   The lexicographic order forces that any face is larger than
 *   a coface.
 *
 *  \tparam Permutahdral_representation_ Needs to be 
 *   Permutahedral_representation<Vertex_, Ordered_set_partition_>
 *
 *  \ingroup coxeter_triangulation
 */
template <class Permutahedral_representation_>
struct Simplex_comparator {

  /** \brief Comparison between two permutahedral representations.
   *  Both permutahedral representations need to be valid and
   *  the vertices of both permutahedral representations need to be of the same size.
   */
  bool operator() (const Permutahedral_representation_& lhs,
		   const Permutahedral_representation_& rhs) const {
    if (lhs.vertex() < rhs.vertex())
      return true;
    if (lhs.vertex() > rhs.vertex())
      return false;
    if (lhs.partition() < lhs.partition())
      return true;
    if (lhs.partition() > lhs.partition())
      return false;
    return false;
  }
};

} // namespace coxeter_triangulation 

} // namespace Gudhi

#endif
