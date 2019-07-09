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
#include <iostream>

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

  typedef Permutahedral_representation<Vertex_, Ordered_set_partition_> Self;

public:
  /** \brief Type of the vertex. */
  typedef Vertex_ Vertex;
  
  /** \brief Type of the ordered partition. */
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

  /** \brief Constructor for an empty permutahedral representation that does not correspond
   *  to any simplex.
   */
  Permutahedral_representation() {}
  
  /** \brief Dimension of the simplex. */
  unsigned dimension() const {
    return partition_.size() - 1;
  }

  /** \brief Lexicographically-minimal vertex. */
  Vertex& vertex() {
    return vertex_;
  }

  /** \brief Lexicographically-minimal vertex. */
  const Vertex& vertex() const {
    return vertex_;
  }

  /** \brief Ordered set partition. */
  OrderedSetPartition& partition() {
    return partition_;
  }

  /** \brief Identifying vertex. */
  const OrderedSetPartition& partition() const {
    return partition_;
  }

  /** \brief Equality operator.
   * Returns true if an only if both vertex and the ordered set partition coincide.
   */
  bool operator==(const Permutahedral_representation& other) const {
    if (dimension() != other.dimension())
      return false;
    if (vertex_ != other.vertex_)
      return false;
    for (std::size_t k = 0; k < partition_.size(); ++k)
      if (partition_[k] != other.partition_[k])
	return false;
    return true;
  }

  /** \brief Inequality operator.
   * Returns true if an only if either vertex or the ordered set partition are different.
   */
  bool operator!=(const Permutahedral_representation& other) const {
    return !(*this == other);
  }

  typedef Gudhi::coxeter_triangulation::Vertex_iterator<Self> Vertex_iterator;
  typedef boost::iterator_range<Vertex_iterator> Vertex_range;
  /** \brief Returns a range of vertices of the simplex.
   * The type of vertices is Vertex.
   */
  Vertex_range vertex_range() const {
    return Vertex_range(Vertex_iterator(*this),
  			Vertex_iterator());
  }

  typedef Gudhi::coxeter_triangulation::Face_iterator<Self> Face_iterator;
  typedef boost::iterator_range<Face_iterator> Face_range;
  /** \brief Returns a range of permutahedral representations of faces of the simplex.
   * @param[in] value_dim The dimension of the faces. Must be between 0 and the dimension of the simplex.
   */
  Face_range face_range(std::size_t value_dim) const {
    return Face_range(Face_iterator(*this, value_dim),
		      Face_iterator());
  }

  /** \brief Returns a range of permutahedral representations of facets of the simplex.
   * The dimension of the simplex must be strictly positive.
   */
  Face_range facet_range() const {
    return Face_range(Face_iterator(*this, dimension()-1),
		      Face_iterator());
  }
  
  typedef Gudhi::coxeter_triangulation::Coface_iterator<Self> Coface_iterator;
  typedef boost::iterator_range<Coface_iterator> Coface_range;
  /** \brief Returns a range of permutahedral representations of cofaces of the simplex.
   * @param[in] value_dim The dimension of the cofaces. Must be between the dimension of the simplex and the ambient dimension (the size of the vertex).
   */
  Coface_range coface_range(std::size_t value_dim) const {
    return Coface_range(Coface_iterator(*this, value_dim),
			Coface_iterator());
  }

  /** \brief Returns a range of permutahedral representations of cofacets of the simplex.
   * The dimension of the simplex must be strictly different from the ambient dimension (the size of the vertex).
   */
  Coface_range cofacet_range() const {
    return Coface_range(Coface_iterator(*this, dimension()+1),
			Coface_iterator());
  }

  
    bool is_face_of(const Permutahedral_representation& other) const {
    if (other.dimension() < dimension())
      return false;
    if (other.vertex_.size() != vertex_.size())
      std::cerr << "Error: Permutahedral_representation::is_face_of: incompatible ambient dimensions.\n";
    std::size_t d = vertex_.size();    

    /* check the vertex */    
    Vertex v_copy = vertex_;
    std::size_t other_i = 0, self_i = 0;
    bool other_chosen = false;
    bool other_start = 0;
    while (other_i < other.partition_.size()) {
      auto part_it = std::find(partition_[partition_.size()-1].begin(),
			       partition_[partition_.size()-1].end(),
			       *other.partition_[other_i].begin());      
      if (!other_chosen && part_it != partition_[partition_.size()-1].end()) {
	bool cyclic_shifted = false;
	for (const auto& i: other.partition_[other_i])
	  if (i == d) {
	    if (cyclic_shifted)
	      return false;
	    if (!other_chosen) {
	      other_chosen = true;
	      other_start = other_i;
	    }
	  }
	  else if (vertex_[i] == other.vertex_[i]) {
	    if (cyclic_shifted)
	      return false;
	    other_chosen = true;
	    other_start = other_i;
	  }
	  else if (vertex_[i] == other.vertex_[i] + 1) {
	    if (other_chosen)
	      return false;
	    cyclic_shifted = true;
	  }
	  else
	    return false;	
      }
      else {
	if (!other_chosen) {
	  other_chosen = true;
	  other_start = other_i;
	}
	for (const auto& i: other.partition_[other_i])
	  if (i != d && vertex_[i] != other.vertex_[i])
	    return false;		
      }
      other_i++;
    }
    other_i = other_start;
    if (other_i == other.partition_.size())
      std::cerr << "Error: Permutahedral_representation::is_face_of: wrong partition.\n";

    std::unordered_set<std::size_t> part;
    for (const auto& i: partition_[self_i])
      part.insert(i);
    
    /* check the partitions */    
    while (self_i < partition_.size()) {
      for (const auto& i: other.partition_[other_i]) {
	auto part_it = part.find(i);
	if (part_it == part.end()) //mismatch
	  return false;
	else
	  part.erase(part_it);
      }
      other_i = (other_i + 1) % other.partition_.size();
      if (part.empty()) {
	self_i++;
	if (self_i < partition_.size())
	  for (const auto& i: partition_[self_i])
	    part.insert(i);
      }
    }
    return true;
  }
  
private:
  Vertex vertex_;
  OrderedSetPartition partition_;

};

/** \brief Print a permutahedral representation to a stream.
 * \ingroup coxeter_triangulation
 *
 * @param[in] os The output stream.
 * @param[in] simplex A simplex represented by its permutahedral representation.
 */
template <class Vertex,
	  class OrderedSetPartition>
std::ostream& operator<<(std::ostream& os,
			 const Permutahedral_representation<Vertex, OrderedSetPartition>& simplex) {
  // vertex part
  os << "(";
  if (simplex.vertex().empty()) {
    std::cout << ")";
    return os;
  }
  auto v_it = simplex.vertex().begin();
  os << *v_it++;
  for (; v_it != simplex.vertex().end(); ++v_it)
    os << ", " << *v_it;
  os << ")";
  
  // ordered partition part
  using Part = typename OrderedSetPartition::value_type;
  auto print_part =
    [&os](const Part& p) {
      os << "{";
      if (p.empty()) {
	std::cout << "}";
      }
      auto p_it = p.begin();
      os << *p_it++;
      for (; p_it != p.end(); ++p_it)
	os << ", " << *p_it;  
      os << "}";
    };
  os << " [";
  if (simplex.partition().empty()) {
    std::cout << "]";
    return os;
  }
  auto o_it = simplex.partition().begin();
  print_part(*o_it++);
  for (; o_it != simplex.partition().end(); ++o_it) {
    os << ", ";
    print_part(*o_it);
  }
  os << "]";
  return os;
}

} // namespace coxeter_triangulation

} // namespace Gudhi

#endif // PERMUTAHEDRAL_REPRESENTATION_H_
