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

#ifndef PERMUTAHEDRAL_REPRESENTATION_H_
#define PERMUTAHEDRAL_REPRESENTATION_H_

#include <gudhi/Permutahedral_representation/Permutahedral_representation_iterators.h>

#include <unordered_set>

namespace Gudhi {

namespace coxeter_triangulation {

/** 
 * \class Permutahedral_representation
 * \brief A class that stores the permutahedral representation of a simplex
 * in a Coxeter triangulation or a Freudenthal-Kuhn triangulation.
 *
 * \ingroup coxeter_triangulation
 *
 * \details The data structure is a record consisting of a range that
 * represents the vertex and a range that represents the ordered set
 * partition, both of which identify the simplex in the triangulation.
 *
 * \tparam Vertex_ needs to be a random-access range.
 * \tparam Ordered_set_partition_ needs to be a a random-access range that consists of
 * random-access ranges. 
 */
template <class Vertex_,
	  class Ordered_set_partition_>
class Permutahedral_representation {

public:
  /**
   * \brief Type of the vertex.
   */
  typedef Vertex_ Vertex;
  
  /**
   * \brief Type of the ordered partition.
   */
  typedef Ordered_set_partition_ OrderedSetPartition;

  /** \brief Permutahedral_representation constructor from a vertex and an ordered set partition.
   *
   * @param[in] vertex Vertex.
   * @param[in] partition Ordered set partition.
   * 
   * \details If the size of vertex is d, the ranges in partition must consist
   * of the integers 0,...,d without repetition or collision between the ranges.
   */
  Permutahedral_representation(const Vertex& vertex, const OrderedSetPartition& partition)
    : vertex_(vertex), partition_(partition) {}

  /** \brief Dimension of the simplex.
   */
  unsigned dimension() const {
    return partition_.size() - 1;
  }

  /** \brief Identifying vertex.
   */
  Vertex& vertex() {
    return vertex_;
  }

  /** \brief Identifying vertex.
   */
  const Vertex& vertex() const {
    return vertex_;
  }
  
private:
  Vertex vertex_;
  Ordered_set_partition partition_;

};

} // namespace coxeter_triangulation

} // namespace Gudhi

#endif // PERMUTAHEDRAL_REPRESENTATION_H_
