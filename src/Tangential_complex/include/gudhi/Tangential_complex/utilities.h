/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Clement Jamin
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef TANGENTIAL_COMPLEX_UTILITIES_H_
#define TANGENTIAL_COMPLEX_UTILITIES_H_

#include <CGAL/Dimension.h>
#include <CGAL/Combination_enumerator.h>
#include <CGAL/IO/Triangulation_off_ostream.h>

#include <boost/container/flat_set.hpp>

#include <Eigen/Core>
#include <Eigen/Eigen>

#include <set>
#include <vector>
#include <array>
#include <fstream>
#include <atomic>
#include <cmath>  // for std::sqrt

namespace Gudhi {
namespace tangential_complex {
namespace internal {

// Provides copy constructors to std::atomic so that
// it can be used in a vector
template <typename T>
struct Atomic_wrapper
: public std::atomic<T> {
  typedef std::atomic<T> Base;

  Atomic_wrapper() { }

  Atomic_wrapper(const T &t) : Base(t) { }

  Atomic_wrapper(const std::atomic<T> &a) : Base(a.load()) { }

  Atomic_wrapper(const Atomic_wrapper &other) : Base(other.load()) { }

  Atomic_wrapper &operator=(const T &other) {
    Base::store(other);
    return *this;
  }

  Atomic_wrapper &operator=(const std::atomic<T> &other) {
    Base::store(other.load());
    return *this;
  }

  Atomic_wrapper &operator=(const Atomic_wrapper &other) {
    Base::store(other.load());
    return *this;
  }
};

// Modifies v in-place
template <typename K>
typename K::Vector_d& normalize_vector(typename K::Vector_d& v,
                                       K const& k) {
  v = k.scaled_vector_d_object()(
                                 v, typename K::FT(1) / std::sqrt(k.squared_length_d_object()(v)));
  return v;
}

template<typename Kernel>
struct Basis {
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_d Point;
  typedef typename Kernel::Vector_d Vector;
  typedef typename std::vector<Vector>::const_iterator const_iterator;

  std::size_t m_origin;
  std::vector<Vector> m_vectors;

  std::size_t origin() const {
    return m_origin;
  }

  void set_origin(std::size_t o) {
    m_origin = o;
  }

  const_iterator begin() const {
    return m_vectors.begin();
  }

  const_iterator end() const {
    return m_vectors.end();
  }

  std::size_t size() const {
    return m_vectors.size();
  }

  Vector& operator[](const std::size_t i) {
    return m_vectors[i];
  }

  const Vector& operator[](const std::size_t i) const {
    return m_vectors[i];
  }

  void push_back(const Vector& v) {
    m_vectors.push_back(v);
  }

  void reserve(const std::size_t s) {
    m_vectors.reserve(s);
  }

  Basis() { }

  Basis(std::size_t origin) : m_origin(origin) { }

  Basis(std::size_t origin, const std::vector<Vector>& vectors)
      : m_origin(origin), m_vectors(vectors) { }

  int dimension() const {
    return static_cast<int> (m_vectors.size());
  }
};

// 1st line: number of points
// Then one point per line
template <typename Kernel, typename Point_range>
std::ostream &export_point_set(
                               Kernel const& k,
                               Point_range const& points,
                               std::ostream & os,
                               const char *coord_separator = " ") {
  // Kernel functors
  typename Kernel::Construct_cartesian_const_iterator_d ccci =
      k.construct_cartesian_const_iterator_d_object();

  os << points.size() << "\n";

  typename Point_range::const_iterator it_p = points.begin();
  typename Point_range::const_iterator it_p_end = points.end();
  // For each point p
  for (; it_p != it_p_end; ++it_p) {
    for (auto it = ccci(*it_p); it != ccci(*it_p, 0); ++it)
      os << CGAL::to_double(*it) << coord_separator;

    os << "\n";
  }

  return os;
}

// Compute all the k-combinations of elements
// Output_iterator::value_type must be
// boost::container::flat_set<std::size_t>
template <typename Elements_container, typename Output_iterator>
void combinations(const Elements_container elements, int k,
                  Output_iterator combinations) {
  std::size_t n = elements.size();
  std::vector<bool> booleans(n, false);
  std::fill(booleans.begin() + n - k, booleans.end(), true);
  do {
    boost::container::flat_set<std::size_t> combination;
    typename Elements_container::const_iterator it_elt = elements.begin();
    for (std::size_t i = 0; i < n; ++i, ++it_elt) {
      if (booleans[i])
        combination.insert(*it_elt);
    }
    *combinations++ = combination;
  } while (std::next_permutation(booleans.begin(), booleans.end()));
}

}  // namespace internal
}  // namespace tangential_complex
}  // namespace Gudhi

#endif  // TANGENTIAL_COMPLEX_UTILITIES_H_
