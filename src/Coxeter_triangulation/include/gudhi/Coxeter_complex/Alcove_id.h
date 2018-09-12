#ifndef ALCOVE_ID_H_
#define ALCOVE_ID_H_

#include <vector>

// The structure for the coordinates of a given relatively open simplex (alcove)
namespace Gudhi {
struct Alcove_id {

  typedef long int value_type;

  // typedef typename std::vector<int>::iterator iterator;
  typedef typename std::vector<value_type>::const_iterator const_iterator;
  typedef typename std::vector<bool>::const_iterator fixed_const_iterator;

  Alcove_id()
    : level_(0), dimension_(0) {}

  Alcove_id(double level)
    : level_(level), dimension_(0) {}

  Alcove_id(double level, unsigned dimension)
    : level_(level), dimension_(dimension) {}
 
  // Alcove_id& operator=(const Alcove_id& other) {
  //   return *this;
  // }
  
  value_type operator[] (std::size_t i) const {
    return coords_[i];
  }
  
  void push_back(value_type value, bool fixed = false) {
    coords_.push_back(value);
    fixed_.push_back(fixed);
  }

  void pop_back() {
    coords_.pop_back();
    fixed_.pop_back();
  }
    
  const_iterator begin() const {
    return coords_.begin();
  }

  const_iterator end() const {
    return coords_.end();
  }

  fixed_const_iterator fixed_begin() const {
    return fixed_.begin();
  }

  fixed_const_iterator fixed_end() const {
    return fixed_.end();
  }
  
  std::size_t size() const {
    return coords_.size();
  }

  std::size_t empty() const {
    return coords_.size() == 0;
  }
  
  void reserve(std::size_t new_cap) {
    coords_.reserve(new_cap);
    fixed_.reserve(new_cap);
  }
  
  void resize(std::size_t new_size) {
    coords_.resize(new_size);
    fixed_.resize(new_size);
  }

  void clear() {
    resize(0);
  }
  
  double level() const {
    return level_;
  }
  
  double dimension() const {
    return dimension_;
  }
  
  bool is_fixed(std::size_t i) const {
    return fixed_[i];
  }

  void set_dimension(std::size_t dim) {
    dimension_ = dim;
  }

  bool operator==(const Alcove_id& other) const {
    if (this->size() != other.size())
      return false;
    for (std::size_t k = 0; k < this->size(); ++k)
      if ((*this)[k] != other[k] || this->is_fixed(k) != other.is_fixed(k))
        return false;
    return true;
  }
  
  double level_;
  unsigned dimension_;
  std::vector<value_type> coords_;
  std::vector<bool> fixed_;
};

  bool operator< (const Alcove_id& lhs, const Alcove_id& rhs) {
    return lhs.coords_ < rhs.coords_;
  }

  std::ostream& operator<<(std::ostream & os, const Alcove_id& a_id) {
    os << "[";
    if (a_id.empty()) {
      std::cout << "]_" << a_id.dimension();
      return os;
    }
    if (a_id.is_fixed(0))
        os << "\033[1;31m" << a_id[0] << "\033[0m";
      else
        os << a_id[0];
    for (std::size_t i = 1; i < a_id.size(); ++i)
      if (a_id.is_fixed(i))
        os << ", \033[1;31m" << a_id[i] << "\033[0m";
      else
        os << ", " << a_id[i];
    std::cout << "]_" << a_id.dimension();
    return os;
  }

  
}

#endif
