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

#ifndef PERMUTAHEDRAL_REPRESENTATION_PERMUTAHEDRAL_REPRESENTATION_ITERATORS_H_
#define PERMUTAHEDRAL_REPRESENTATION_PERMUTAHEDRAL_REPRESENTATION_ITERATORS_H_

#include <gudhi/Permutahedral_representation/Size_range.h>
#include <gudhi/Permutahedral_representation/Ordered_set_partition_iterator.h>
#include <gudhi/Permutahedral_representation/Integer_combination_iterator.h>
#include <gudhi/Permutahedral_representation/Combination_iterator.h>
#include <gudhi/Permutahedral_representation/face_from_indices.h>
#include <boost/iterator/iterator_facade.hpp>

#include <vector>
#include <iostream>
#include <algorithm>

namespace Gudhi {

namespace coxeter_triangulation {

/* \addtogroup coxeter_triangulation
 * Iterator types for Permutahedral_representation
 * @{
 */

/* \brief Iterator over the vertices of a simplex
 * represented by its permutahedral representation.
 *
 * Forward iterator, 'value_type' is Permutahedral_representation::Vertex.*/
template <class Permutahedral_representation>
class Vertex_iterator : public boost::iterator_facade< Vertex_iterator<Permutahedral_representation>,
						       typename Permutahedral_representation::Vertex const,
						       boost::forward_traversal_tag> {
protected:
  friend class boost::iterator_core_access;

  using Vertex = typename Permutahedral_representation::Vertex;
  using Ordered_partition = typename Permutahedral_representation::OrderedSetPartition;

  typedef Vertex value_t;

  
  bool equal(Vertex_iterator const& other) const {
    return (is_end_ && other.is_end_);
  }

  value_t const& dereference() const {
    return value_;
  }

  void update_value() {
    std::size_t d = value_.size();
    for (auto i: *o_it_)
      if (i != d)
	value_[i]++;
      else
	for (std::size_t j = 0; j < d; j++)
	  value_[j]--;
  }
		       
  void increment() {
    if (is_end_)
      return;
    update_value();
    if (++o_it_ == o_end_)
      is_end_ = true;
  }

public:  
  Vertex_iterator(const Permutahedral_representation& simplex)
    :
    o_it_(simplex.partition().begin()),
    o_end_(simplex.partition().end()),
    value_(simplex.vertex()),
    is_end_(o_it_ == o_end_)
  {}

  Vertex_iterator() : is_end_(true) {}
  
  
protected:
  typename Ordered_partition::const_iterator o_it_, o_end_; 
  value_t value_;
  bool is_end_;

}; // Vertex_iterator

/*---------------------------------------------------------------------------*/
/* \brief Iterator over the k-faces of a simplex
 *  given by its permutahedral representation.
 *
 * Forward iterator, value_type is Permutahedral_representation. */
template <class Permutahedral_representation>
class Face_iterator : public boost::iterator_facade< Face_iterator<Permutahedral_representation>,
						     Permutahedral_representation const,
						     boost::forward_traversal_tag> {
  typedef Permutahedral_representation value_t;
  
protected:
  friend class boost::iterator_core_access;

  using Vertex = typename Permutahedral_representation::Vertex;
  using Ordered_partition = typename Permutahedral_representation::OrderedSetPartition;

  bool equal(Face_iterator const& other) const {
    return (is_end_ && other.is_end_);
  }

  value_t const& dereference() const {
    return value_;
  }

  void increment() {
    if (++c_it_ == c_end_) {
      is_end_ = true;
      return;
    }
    update_value();
  }

  void update_value() {
    // Combination *c_it_ is supposed to be sorted in increasing order
    value_ = face_from_indices<Permutahedral_representation>(simplex_, *c_it_);
  }
  
public:
  
  Face_iterator(const Permutahedral_representation& simplex, const uint& k)
    : simplex_(simplex),
      k_(k),
      l_(simplex.dimension()),
      c_it_(l_+1, k_+1),
      is_end_(k_ > l_),
      value_({Vertex(simplex.vertex().size()), Ordered_partition(k+1)})
  {
    update_value();
  }

  // Used for the creating an end iterator
  Face_iterator() : is_end_(true) {}
  
protected:
  Permutahedral_representation simplex_; // Input simplex
  uint k_;
  uint l_; // Dimension of the input simplex
  Combination_iterator c_it_, c_end_; // indicates the vertices in the current face

  bool is_end_;   // is true when the current permutation is the final one
  value_t value_; // the dereference value

}; // Face_iterator

/*---------------------------------------------------------------------------*/
/* \brief Iterator over the k-cofaces of a simplex
 *  given by its permutahedral representation.
 *
 * Forward iterator, value_type is Permutahedral_representation. */
template <class Permutahedral_representation>
class Coface_iterator : public boost::iterator_facade< Coface_iterator<Permutahedral_representation>,
						       Permutahedral_representation const,
						       boost::forward_traversal_tag> {
  typedef Permutahedral_representation value_t;
  
protected:
  friend class boost::iterator_core_access;

  using Vertex = typename Permutahedral_representation::Vertex;
  using Ordered_partition = typename Permutahedral_representation::OrderedSetPartition;

  bool equal(Coface_iterator const& other) const {
    return (is_end_ && other.is_end_);
  }

  value_t const& dereference() const {
    return value_;
  }

  void increment() {
    uint i = 0;
    for (; i < k_+1; i++) {
      if (++(o_its_[i]) != o_end_)
	break;
    }
    if (i == k_+1) {
      if (++i_it_ == i_end_) {
	is_end_ = true;
	return;
      }
      o_its_.clear();
      for (uint j = 0; j < k_ + 1; j++)
	o_its_.emplace_back(Ordered_set_partition_iterator(simplex_.partition()[j].size(), (*i_it_)[j]+1));
    }
    else
      for (uint j = 0; j < i; j++)
	o_its_[j].reinitialize();
    update_value();
  }

  void update_value() {
    value_.vertex() = simplex_.vertex();
    for (auto& p: value_.partition())
      p.clear();
    uint u_ = 0; // the part in o_its_[k_] that contains t_
    for (; u_ <= (*i_it_)[k_]; u_++) {
      auto range = (*o_its_[k_])[u_];
      if (std::find(range.begin(), range.end(), t_) != range.end())
	break;
    }
    uint i = 0;
    for (uint j = u_+1; j <= (*i_it_)[k_]; j++, i++)
      for (uint b: (*o_its_[k_])[j]) {
	uint c = simplex_.partition()[k_][b];
        value_.partition()[i].push_back(c);
	value_.vertex()[c]--;
      }
    for (uint h = 0; h < k_; h++)
      for (uint j = 0; j <= (*i_it_)[h]; j++, i++) {
	for (uint b: (*o_its_[h])[j])
	  value_.partition()[i].push_back(simplex_.partition()[h][b]);
      }
    for (uint j = 0; j <= u_; j++, i++)
      for (uint b: (*o_its_[k_])[j])
	value_.partition()[i].push_back(simplex_.partition()[k_][b]);
    // sort the values in each part (probably not needed)
    for (auto& part: value_.partition())
      std::sort(part.begin(), part.end());
  }
  
public:
  
  Coface_iterator(const Permutahedral_representation& simplex, const uint& l)
    : simplex_(simplex),
      d_(simplex.vertex().size()),
      l_(l),
      k_(simplex.dimension()),
      i_it_(l_-k_ , k_+1, Size_range<Ordered_partition>(simplex.partition())),
      is_end_(k_ > l_),
      value_({Vertex(d_), Ordered_partition(l_+1)})
  {
    uint j = 0;
    for (; j < simplex_.partition()[k_].size(); j++)
      if (simplex_.partition()[k_][j] == d_) {
	t_ = j;
	break;
      }
    if (j == simplex_.partition()[k_].size()) {
      std::cerr << "Coface iterator: the argument simplex is not a permutahedral representation\n";
      is_end_ = true;
      return;
    }
    for (uint i = 0; i < k_+1; i++)
      o_its_.emplace_back(Ordered_set_partition_iterator(simplex_.partition()[i].size(), (*i_it_)[i]+1));
    update_value();
  }

  // Used for the creating an end iterator
  Coface_iterator() : is_end_(true) {}
  
protected:
  Permutahedral_representation simplex_; // Input simplex
  uint d_; // Ambient dimension
  uint l_; // Dimension of the coface
  uint k_; // Dimension of the input simplex
  uint t_; // The position of d in simplex_.partition()[k_]
  Integer_combination_iterator i_it_, i_end_; // indicates in how many parts each simplex_[i] is subdivided
  std::vector<Ordered_set_partition_iterator> o_its_; // indicates subdivision for each simplex_[i]
  Ordered_set_partition_iterator o_end_;      // one end for all o_its_

  bool is_end_;   // is true when the current permutation is the final one
  value_t value_; // the dereference value

}; // Coface_iterator

} // namespace coxeter_triangulation

} // namespace Gudhi
  
#endif
