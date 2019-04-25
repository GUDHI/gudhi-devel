#ifndef COFACE_ITERATOR_H_
#define COFACE_ITERATOR_H_

#include <vector>
#include <limits>
#include <iostream>

#include <gudhi/Coxeter_triangulation/Size_range.h>
#include <gudhi/Coxeter_triangulation/Integer_combination_iterator.h>
#include <gudhi/Coxeter_triangulation/Ordered_set_partition_iterator.h>
#include <boost/range/iterator_range.hpp>

typedef unsigned uint;

/** \brief Coface iterator class.
 *   
*/
namespace Gudhi {

template <class Freudenthal_representation>
class Coface_iterator : public boost::iterator_facade< Coface_iterator<Freudenthal_representation>,
						       Freudenthal_representation const,
						       boost::forward_traversal_tag> {
  typedef Freudenthal_representation value_t;
  
protected:
  friend class boost::iterator_core_access;

  using Vertex = typename Freudenthal_representation::vertex_type;
  using Ordered_partition = typename Freudenthal_representation::partition_type;

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
	o_its_.emplace_back(Ordered_set_partition_iterator(fr_.partition[j].size(), (*i_it_)[j]+1));
    }
    else
      for (uint j = 0; j < i; j++)
	o_its_[j].reinitialize();
    update_value();
  }

  void update_value() {
    value_.vertex = fr_.vertex;
    for (auto& p: value_.partition)
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
	uint c = fr_.partition[k_][b];
        value_.partition[i].push_back(c);
	if (c != d_)
	  value_.vertex[j]++;
	else
	  for (uint l = 0; l < d_; l++)
	    value_.vertex[l]--;
      }
    for (uint h = 0; h < k_; h++)
      for (uint j = 0; j <= (*i_it_)[h]; j++, i++) {
	for (uint b: (*o_its_[h])[j])
	  value_.partition[i].push_back(fr_.partition[h][b]);
      }
    for (uint j = 0; j <= u_; j++, i++)
      for (uint b: (*o_its_[k_])[j])
	value_.partition[i].push_back(fr_.partition[k_][b]);
  }
  
public:
  
  Coface_iterator(const Freudenthal_representation& fr, const uint& l)
    : fr_(fr),
      d_(fr.vertex.size()),
      l_(l),
      k_(fr.dimension()),
      i_it_(l_-k_ , k_+1, Size_range<Ordered_partition>(fr.partition)),
      is_end_(k_ > l_),
      value_({Vertex(d_), Ordered_partition(l_+1)})
  {
    uint j = 0;
    for (; j < fr_.partition[k_].size(); j++)
      if (fr_.partition[k_][j] == d_) {
	t_ = j;
	break;
      }
    if (j == fr_.partition[k_].size()) {
      std::cerr << "Coface iterator: the argument fr is not a Freudenthal's representation\n";
      is_end_ = true;
      return;
    }
    for (uint i = 0; i < k_+1; i++)
      o_its_.emplace_back(Ordered_set_partition_iterator(fr_.partition[i].size(), (*i_it_)[i]+1));
    update_value();
  }

  // Used for the creating an end iterator
  Coface_iterator() : is_end_(true) {}
  
protected:
  Freudenthal_representation fr_; // Input simplex
  uint d_; // Ambient dimension
  uint l_; // Dimension of the coface
  uint k_; // Dimension of the input simplex
  uint t_; // The position of d in fr_.partition[k_]
  Integer_combination_iterator i_it_, i_end_; // indicates in how many parts each fr_[i] is subdivided
  std::vector<Ordered_set_partition_iterator> o_its_; // indicates subdivision for each fr_[i]
  Ordered_set_partition_iterator o_end_;      // one end for all o_its_

  bool is_end_;   // is true when the current permutation is the final one
  value_t value_; // the dereference value
};

}

#endif
