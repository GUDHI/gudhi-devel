#ifndef FACE_ITERATOR_H_
#define FACE_ITERATOR_H_

#include <vector>
#include <limits>
#include <gudhi/Coxeter_triangulation/Combination_iterator.h>
#include <gudhi/Coxeter_triangulation/Permutation_iterator.h>
#include <gudhi/Coxeter_triangulation/Set_partition_iterator.h>
#include <boost/range/iterator_range.hpp>

typedef unsigned uint;

/** \brief Face iterator class.
 *   
*/
namespace Gudhi {

template <class Freudenthal_representation>
class Face_iterator : public boost::iterator_facade< Face_iterator<Freudenthal_representation>,
						     Freudenthal_representation const,
						     boost::forward_traversal_tag> {
  typedef Freudenthal_representation value_t;
  
protected:
  friend class boost::iterator_core_access;

  using Vertex = typename Freudenthal_representation::vertex_type;
  using Ordered_partition = typename Freudenthal_representation::partition_type;

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
    uint d = fr_.vertex.size();
    value_.vertex = fr_.vertex;
    for (auto& p: value_.partition)
      p.clear();
    for (uint h = 1; h < k_+1; h++)
      for (uint i = (*c_it_)[h-1]; i < (*c_it_)[h]; i++)
	for (uint j: fr_.partition[i])
	  value_.partition[h-1].push_back(j);
    for (uint i = (*c_it_)[k_]; i < l_+1; i++)
      for (uint j: fr_.partition[i])
    	value_.partition[k_].push_back(j);
    for (uint i = 0; i < (*c_it_)[0]; i++)
      for (uint j: fr_.partition[i]) {
	if (j != d)
	  value_.vertex[j]++;
	else
	  for (uint l = 0; l < d; l++)
	    value_.vertex[l]--;
	value_.partition[k_].push_back(j);
      }
  }
  
public:
  
  Face_iterator(const Freudenthal_representation& fr, const uint& k)
    : fr_(fr),
      k_(k),
      l_(fr.dimension()),
      c_it_(l_+1, k_+1),
      is_end_(k_ > l_),
      value_({Vertex(fr.vertex.size()), Ordered_partition(k+1)})
  {
    update_value();
  }

  // Used for the creating an end iterator
  Face_iterator() : is_end_(true) {}
  
protected:
  Freudenthal_representation fr_; // Input simplex
  uint k_;
  uint l_; // Dimension of the input simplex
  Combination_iterator c_it_, c_end_; // indicates the vertices in the current face

  bool is_end_;   // is true when the current permutation is the final one
  value_t value_; // the dereference value
};

}

#endif
