#ifndef VERTEX_ITERATOR_H_
#define VERTEX_ITERATOR_H_

#include <boost/iterator/iterator_facade.hpp>

namespace Gudhi {

template <class Freudenthal_representation>
class Vertex_iterator : public boost::iterator_facade< Vertex_iterator<Freudenthal_representation>,
						       typename Freudenthal_representation::vertex_type const,
						       boost::forward_traversal_tag> {
protected:
  friend class boost::iterator_core_access;

  using Vertex = typename Freudenthal_representation::vertex_type;
  using Ordered_partition = typename Freudenthal_representation::partition_type;

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
    if (++o_it_ == o_end_) {
      is_end_ = true;
      return;
    }
    update_value();
  }

public:  
  Vertex_iterator(const Freudenthal_representation& fr)
    :
    o_it_(fr.partition.begin()),
    o_end_(fr.partition.end()),
    value_(fr.vertex),
    is_end_(o_it_ == o_end_)
  {
    update_value();
  }

  Vertex_iterator() : is_end_(true) {}
  
  
protected:
  typename Ordered_partition::const_iterator o_it_, o_end_; 
  value_t value_;
  bool is_end_;
};

}
  
#endif
