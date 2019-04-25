#ifndef SET_PARTITION_ITERATOR_H_
#define SET_PARTITION_ITERATOR_H_

#include <vector>
#include <limits>
#include <boost/range/iterator_range.hpp>

typedef unsigned uint;

/** \brief Class that allows the user to generate set partitions of a set {0,...,n-1} in k parts.
 *   
*/
class Set_partition_iterator : public boost::iterator_facade< Set_partition_iterator,
							      std::vector<std::vector<uint> > const,
							      boost::forward_traversal_tag> {
  typedef std::vector<std::vector<uint> > value_t;
  
protected:
  friend class boost::iterator_core_access;

  bool equal(Set_partition_iterator const& other) const {
    return (is_end_ && other.is_end_);
  }

  value_t const& dereference() const {
    return value_;
  }

  void update_value() {
    for (uint i = 0; i < k_; i++)
      value_[i].clear();
    for (uint i = 0; i < n_; i++)
      value_[rgs_[i]].push_back(i);
  }
  
  void increment() {
    if (k_ <= 1) {
      is_end_ = true;
      return;
    }
    uint i = n_ - 1;
    while (rgs_[i] + 1 > max_[i] ||
	   rgs_[i] + 1 >= k_)
      i--;
    if (i == 0) {
      is_end_ = true;
      return;
    }
    rgs_[i]++;
    uint mm = max_[i];
    mm += (rgs_[i] >= mm);
    max_[i+1] = mm;
    while (++i < n_) {
      rgs_[i] = 0;
      max_[i+1] = mm;
    }
    uint p = k_;
    if (mm < p)
      do {
	max_[i] = p;
	--i;
	--p;
	rgs_[i] = p;
      } while (max_[i] < p);
    update_value();
  }
  
public:
  
  Set_partition_iterator(const uint& n, const uint& k)
    :
    value_(k),
    rgs_(n, 0),
    max_(k+1),
    is_end_(n == 0),
    n_(n),
    k_(k)
  {
    max_[0] = std::numeric_limits<uint>::max();
    for (uint i = 0; i <= n-k; ++i)
      value_[0].push_back(i);
    for (uint i = n-k+1, j = 1; i < n; ++i, ++j) {
      rgs_[i] = j;
      value_[j].push_back(i);
    }
    for (uint i = 1; i <= n; i++)
      max_[i] = rgs_[i-1] + 1;
    update_value();
  }

  // Used for the creating an end iterator
  Set_partition_iterator() : is_end_(true), n_(0), k_(0) {}

  void reinitialize() {
    if (n_ > 0)
      is_end_ = false;
    for (uint i = 0; i <= n_-k_; ++i)
      rgs_[i] = 0;
    for (uint i = n_-k_+1, j = 1; i <= n_; ++i, ++j)
      rgs_[i] = j;
    for (uint i = 1; i <= n_; i++)
      max_[i] = rgs_[i-1] + 1;
    update_value();
  }
  
protected:
  value_t value_;         // the dereference value
  std::vector<uint> rgs_; // restricted growth string
  std::vector<uint> max_; // max_[i] = max(rgs_[0],...,rgs[i-1]) + 1
  bool is_end_;   // is true when the current permutation is the final one 

  uint n_;
  uint k_;
};


#endif
