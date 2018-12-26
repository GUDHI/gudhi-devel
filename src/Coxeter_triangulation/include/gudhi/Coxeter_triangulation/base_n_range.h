#ifndef BASE_N_RANGE_H_
#define BASE_N_RANGE_H_

#include <vector>

#include <boost/range/iterator_range.hpp>

class Base_n_iterator : public boost::iterator_facade< Base_n_iterator,
						       std::vector<std::size_t> const,
						       boost::forward_traversal_tag> {
protected:
  friend class boost::iterator_core_access;

  bool equal(Base_n_iterator const& other) const {
    return (is_end_ && other.is_end_);
  }

  std::vector<std::size_t> const& dereference() const {
    return value_;
  }

  void increment() {
    if (is_end_)
      return;
    while (!value_.empty() && ++value_.back() == n_)
      value_.pop_back();
    if (value_.empty()) {
      is_end_ = true;
      return;
    }
    value_.resize(k_);
  }

public:
  Base_n_iterator(std::size_t n, std::size_t k) 
    : n_(n), k_(k), value_(k, 0), is_end_(false) {}

  Base_n_iterator() : is_end_(true) {}
  
  
protected:
  std::size_t n_=0, k_=0;
  std::vector<std::size_t> value_;
  bool is_end_;
};
  
typedef boost::iterator_range<Base_n_iterator> Base_n_range;
Base_n_range base_n_range(std::size_t n, std::size_t k) {
  return Base_n_range(Base_n_iterator(n, k),
		      Base_n_iterator());
}  

#endif
