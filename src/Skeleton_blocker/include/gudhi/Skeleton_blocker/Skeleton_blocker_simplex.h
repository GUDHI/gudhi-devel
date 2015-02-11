/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014  INRIA Sophia Antipolis-Mediterranee (France)
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SRC_SKELETON_BLOCKER_INCLUDE_GUDHI_SKELETON_BLOCKER_SKELETON_BLOCKER_SIMPLEX_H_
#define SRC_SKELETON_BLOCKER_INCLUDE_GUDHI_SKELETON_BLOCKER_SKELETON_BLOCKER_SIMPLEX_H_

#include <cassert>
#include <iostream>
#include <set>
#include <vector>
#include <initializer_list>
#include <string>
#include <algorithm>

namespace Gudhi {

namespace skbl {

/**
 *@brief Abstract simplex used in Skeleton blockers data-structure.
 *
 * An abstract simplex is represented as an ordered set of T elements,
 * each element representing a vertex.
 * 
 * The element representing a vertex can be SkeletonBlockerDS::Vertex_handle or SkeletonBlockerDS::Root_vertex_handle.
 *
 *
 */
template<typename T>

class Skeleton_blocker_simplex {
 private:
  std::set<T> simplex_set;

 public:
  typedef typename T::boost_vertex_handle boost_vertex_handle;

  typedef T Vertex_handle;

  typedef typename std::set<T>::const_iterator Simplex_vertex_const_iterator;
  typedef typename std::set<T>::iterator Simplex_vertex_iterator;

  /** @name Constructors / Destructors / Initialization
   */
  //@{

  // Skeleton_blocker_simplex():simplex_set() {}
  void clear() {
    simplex_set.clear();
  }

  Skeleton_blocker_simplex(std::initializer_list<T>& list) {
    for_each(list.begin(), list.end(), add_vertex);
  }

  template<typename ... Args>
  explicit Skeleton_blocker_simplex(Args ... args) {
    add_vertices(args...);
  }

  template<typename ... Args>
  void add_vertices(T v, Args ... args) {
    add_vertex(v);
    add_vertices(args...);
  }

  void add_vertices(T v) {
    add_vertex(v);
  }

  void add_vertices() {
  }

  /**
   * Initialize a simplex with a string such as {0,1,2}
   */
  explicit Skeleton_blocker_simplex(std::string token) {
    clear();
    if ((token[0] == '{') && (token[token.size() - 1] == '}')) {
      token.erase(0, 1);
      token.erase(token.size() - 1, 1);
      while (token.size() != 0) {
        int coma_position = token.find_first_of(',');
        // cout << "coma_position:"<<coma_position<<endl;
        std::string n = token.substr(0, coma_position);
        // cout << "token:"<<token<<endl;
        token.erase(0, n.size() + 1);
        add_vertex((T) (atoi(n.c_str())));
      }
    }
  }

  //@}

  /** @name Simplex manipulation
   */
  //@{
  /**
   * Add the vertex v to the simplex:
   * Note that adding two times the same vertex is
   * the same that adding it once.
   */
  void add_vertex(T v) {
    simplex_set.insert(v);
  }

  /**
   * Remove the vertex v from the simplex:
   */
  void remove_vertex(T v) {
    simplex_set.erase(v);
  }

  /**
   * Intersects the simplex with the simplex a.
   */
  void intersection(const Skeleton_blocker_simplex & a) {
    std::vector<T> v;
    v.reserve((std::min)(simplex_set.size(), a.simplex_set.size()));

    set_intersection(simplex_set.begin(), simplex_set.end(),
                     a.simplex_set.begin(), a.simplex_set.end(),
                     std::back_inserter(v));
    clear();
    for (auto i : v)
      simplex_set.insert(i);
  }

  /**
   * Substracts a from the simplex.
   */
  void difference(const Skeleton_blocker_simplex & a) {
    std::vector<T> v;
    v.reserve(simplex_set.size());

    set_difference(simplex_set.begin(), simplex_set.end(),
                   a.simplex_set.begin(), a.simplex_set.end(),
                   std::back_inserter(v));

    clear();
    for (auto i : v)
      simplex_set.insert(i);
  }

  /**
   * Add vertices of a to the simplex.
   */
  void union_vertices(const Skeleton_blocker_simplex & a) {
    std::vector<T> v;
    v.reserve(simplex_set.size() + a.simplex_set.size());

    set_union(simplex_set.begin(), simplex_set.end(), a.simplex_set.begin(),
              a.simplex_set.end(), std::back_inserter(v));
    clear();
    simplex_set.insert(v.begin(), v.end());
  }

  typename std::set<T>::const_iterator begin() const {
    return simplex_set.cbegin();
  }

  typename std::set<T>::const_iterator end() const {
    return simplex_set.cend();
  }

  typename std::set<T>::const_reverse_iterator rbegin() const {
    return simplex_set.crbegin();
  }

  typename std::set<T>::const_reverse_iterator rend() const {
    return simplex_set.crend();
  }


  typename std::set<T>::iterator begin() {
    return simplex_set.begin();
  }

  typename std::set<T>::iterator end() {
    return simplex_set.end();
  }

  //@}

  /** @name Queries
   */
  //@{
  /**
   * Returns the dimension of the simplex.
   */
  int dimension() const {
    return (simplex_set.size() - 1);
  }

  bool empty() const {
    return simplex_set.empty();
  }

  /**
   * Returns the first vertex of the (oriented) simplex.
   *
   * Be careful : assumes the simplex is non-empty.
   */
  T first_vertex() const {
    assert(!empty());
    return *(begin());
  }

  /**
   * Returns the last vertex of the (oriented) simplex.
   *
   * Be careful : assumes the simplex is non-empty.
   */
  T last_vertex() const {
    assert(!empty());
    return *(simplex_set.rbegin());
  }
  /**
   * @return true iff the simplex contains the simplex a.
   */
  bool contains(const Skeleton_blocker_simplex & a) const {
    return includes(simplex_set.cbegin(), simplex_set.cend(),
                    a.simplex_set.cbegin(), a.simplex_set.cend());
  }

  /**
   * @return true iff the simplex contains the difference \f$ a \setminus b \f$.
   */
  bool contains_difference(const Skeleton_blocker_simplex& a,
                           const Skeleton_blocker_simplex& b) const {
    auto first1 = begin();
    auto last1 = end();

    auto first2 = a.begin();
    auto last2 = a.end();

    while (first2 != last2) {
      // we ignore vertices of b
      if (b.contains(*first2)) {
        ++first2;
      } else {
        if ((first1 == last1) || (*first2 < *first1))
          return false;
        if (!(*first1 < *first2))
          ++first2;
        ++first1;
      }
    }
    return true;
  }

  /**
   * @return true iff the simplex contains the difference \f$ a \setminus \{ x \} \f$.
   */
  bool contains_difference(const Skeleton_blocker_simplex& a, T x) const {
    auto first1 = begin();
    auto last1 = end();

    auto first2 = a.begin();
    auto last2 = a.end();

    while (first2 != last2) {
      // we ignore vertices x
      if (x == *first2) {
        ++first2;
      } else {
        if ((first1 == last1) || (*first2 < *first1))
          return false;
        if (!(*first1 < *first2))
          ++first2;
        ++first1;
      }
    }
    return true;
  }

  /**
   * @return true iff the simplex contains the difference \f$ a \setminus \{ x,y \} \f$.
   */
  bool contains_difference(const Skeleton_blocker_simplex& a, T x, T y) const {
    auto first1 = begin();
    auto last1 = end();

    auto first2 = a.begin();
    auto last2 = a.end();

    while (first2 != last2) {
      // we ignore vertices of x,y
      if (x == *first2 || y == *first2) {
        ++first2;
      } else {
        if ((first1 == last1) || (*first2 < *first1))
          return false;
        if (!(*first1 < *first2))
          ++first2;
        ++first1;
      }
    }
    return true;
  }

  bool contains(T v) const {
    return (simplex_set.find(v) != simplex_set.end());
  }

  bool disjoint(const Skeleton_blocker_simplex& a) const {
    std::vector<T> v;
    v.reserve(std::min(simplex_set.size(), a.simplex_set.size()));

    set_intersection(simplex_set.cbegin(), simplex_set.cend(),
                     a.simplex_set.cbegin(), a.simplex_set.cend(),
                     std::back_inserter(v));

    return (v.size() == 0);
  }

  bool operator==(const Skeleton_blocker_simplex& other) const {
    return (this->simplex_set == other.simplex_set);
  }

  bool operator!=(const Skeleton_blocker_simplex& other) const {
    return (this->simplex_set != other.simplex_set);
  }

  bool operator<(const Skeleton_blocker_simplex& other) const {
    return (std::lexicographical_compare(this->simplex_set.begin(),
                                         this->simplex_set.end(), other.begin(),
                                         other.end()));
  }

  //@}

  friend std::ostream& operator <<(std::ostream& o,
                                   const Skeleton_blocker_simplex & sigma) {
    bool first = true;
    o << "{";
    for (auto i : sigma) {
      if (first)
        first = false;
      else
        o << ",";
      o << i;
    }
    o << "}";
    return o;
  }
};

}  // namespace skbl

}  // namespace Gudhi

#endif  // SRC_SKELETON_BLOCKER_INCLUDE_GUDHI_SKELETON_BLOCKER_SKELETON_BLOCKER_SIMPLEX_H_

