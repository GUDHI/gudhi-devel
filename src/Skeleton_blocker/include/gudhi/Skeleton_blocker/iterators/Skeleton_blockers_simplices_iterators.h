/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
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

#ifndef SKELETON_BLOCKER_ITERATORS_SKELETON_BLOCKERS_SIMPLICES_ITERATORS_H_
#define SKELETON_BLOCKER_ITERATORS_SKELETON_BLOCKERS_SIMPLICES_ITERATORS_H_

#include <gudhi/Skeleton_blocker_link_complex.h>
#include <gudhi/Skeleton_blocker/Skeleton_blocker_link_superior.h>
#include <gudhi/Skeleton_blocker/internal/Trie.h>
#include <gudhi/Debug_utils.h>

#include <boost/iterator/iterator_facade.hpp>

#include <memory>
#include <list>
#include <iostream>

namespace Gudhi {

namespace skeleton_blocker {

/**
 * Link may be Skeleton_blocker_link_complex<SkeletonBlockerComplex> to iterate over all
 * simplices around a vertex OR
 * Skeleton_blocker_superior_link_complex<SkeletonBlockerComplex> to iterate over all
 * superior vertices around a vertex.
 * The iteration is done by computing a trie with the link and doing a breadth-first traversal
 * of the trie.
 */
template<typename SkeletonBlockerComplex, typename Link>
class Simplex_around_vertex_iterator :
public boost::iterator_facade < Simplex_around_vertex_iterator<SkeletonBlockerComplex, Link>
, typename SkeletonBlockerComplex::Simplex
, boost::forward_traversal_tag
, typename SkeletonBlockerComplex::Simplex
> {
  friend class boost::iterator_core_access;
  typedef SkeletonBlockerComplex Complex;
  typedef typename Complex::Vertex_handle Vertex_handle;
  typedef typename Complex::Edge_handle Edge_handle;
  typedef typename Complex::Simplex Simplex;

  // Link_vertex_handle == Complex_Vertex_handle but this renaming helps avoiding confusion
  typedef typename Link::Vertex_handle Link_vertex_handle;

  typedef typename Gudhi::skeleton_blocker::Trie<Simplex> Trie;

 private:
  const Complex* complex;
  Vertex_handle v;
  std::shared_ptr<Link> link_v;
  std::shared_ptr<Trie> trie;
  // TODO(DS): use a deque instead
  std::list<Trie*> nodes_to_be_seen;

 public:
  Simplex_around_vertex_iterator() : complex(0) { }

  Simplex_around_vertex_iterator(const Complex* complex_, Vertex_handle v_) :
      complex(complex_),
      v(v_),
      link_v(new Link(*complex_, v_)),
      trie(new Trie(v_)) {
    compute_trie_and_nodes_to_be_seen();
  }

  // TODO(DS): avoid useless copy
  // TODO(DS): currently just work if copy begin iterator
  Simplex_around_vertex_iterator(const Simplex_around_vertex_iterator& other) :
      complex(other.complex),
      v(other.v),
      link_v(other.link_v),
      trie(other.trie),
      nodes_to_be_seen(other.nodes_to_be_seen) {
    if (!other.is_end()) {
    }
  }

  /**
   * returns an iterator to the end
   */
  Simplex_around_vertex_iterator(const Complex* complex_, Vertex_handle v_, bool end) :
      complex(complex_),
      v(v_) {
    set_end();
  }

 private:
  void compute_trie_and_nodes_to_be_seen() {
    // once we go through every simplices passing through v0
    // we remove v0. That way, it prevents from passing a lot of times
    // though edges leaving v0.
    // another solution would have been to provides an adjacency iterator
    // to superior vertices that avoids lower ones.
    while (!link_v->empty()) {
      auto v0 = *(link_v->vertex_range().begin());
      trie->add_child(build_trie(v0, trie.get()));
      link_v->remove_vertex(v0);
    }
    nodes_to_be_seen.push_back(trie.get());
  }

  Trie* build_trie(Link_vertex_handle link_vh, Trie* parent) {
    Trie* res = new Trie(parent_vertex(link_vh), parent);
    for (Link_vertex_handle nv : link_v->vertex_range(link_vh)) {
      if (link_vh < nv) {
        Simplex simplex_node_plus_nv(res->simplex());
        simplex_node_plus_nv.add_vertex(parent_vertex(nv));
        if (complex->contains(simplex_node_plus_nv)) {
          res->add_child(build_trie(nv, res));
        }
      }
    }
    return res;
  }

  bool is_node_in_complex(Trie* trie) {
    return true;
  }

  Vertex_handle parent_vertex(Link_vertex_handle link_vh) const {
    return complex->convert_handle_from_another_complex(*link_v, link_vh);
  }

 public:
  friend std::ostream& operator<<(std::ostream& stream, const Simplex_around_vertex_iterator& savi) {
    stream << savi.trie << std::endl;
    stream << "(" << savi.nodes_to_be_seen.size() << ") nodes to see\n";
    return stream;
  }

  bool equal(const Simplex_around_vertex_iterator& other) const {
    bool same_complex = (complex == other.complex);
    if (!same_complex)
      return false;

    if (!(v == other.v))
      return false;

    bool both_empty = nodes_to_be_seen.empty() && other.nodes_to_be_seen.empty();
    if (both_empty)
      return true;

    bool both_non_empty = !nodes_to_be_seen.empty() && !other.nodes_to_be_seen.empty();

    // one is empty the other is not
    if (!both_non_empty) return false;

    bool same_node = (**(nodes_to_be_seen.begin()) == **(other.nodes_to_be_seen.begin()));
    return same_node;
  }

  void increment() {
    assert(!is_end());
    Trie* first_node = nodes_to_be_seen.front();

    nodes_to_be_seen.pop_front();

    for (auto childs : first_node->childs) {
      nodes_to_be_seen.push_back(childs.get());
    }
  }

  Simplex dereference() const {
    assert(!nodes_to_be_seen.empty());
    Trie* first_node = nodes_to_be_seen.front();
    return first_node->simplex();
  }

  Simplex get_trie_address() const {
    assert(!nodes_to_be_seen.empty());
    return nodes_to_be_seen.front();
  }

 private:
  void set_end() {
    nodes_to_be_seen.clear();
  }

  bool is_end() const {
    return nodes_to_be_seen.empty();
  }
};

template<typename SkeletonBlockerComplex>
class Simplex_iterator :
public boost::iterator_facade < Simplex_iterator<SkeletonBlockerComplex>
, typename SkeletonBlockerComplex::Simplex
, boost::forward_traversal_tag
, typename SkeletonBlockerComplex::Simplex
> {
  typedef Skeleton_blocker_link_superior<SkeletonBlockerComplex> Link;

  friend class boost::iterator_core_access;

  template<class SkBlDS> friend class Skeleton_blocker_complex;

  typedef SkeletonBlockerComplex Complex;
  typedef typename Complex::Vertex_handle Vertex_handle;
  typedef typename Complex::Edge_handle Edge_handle;
  typedef typename Complex::Simplex Simplex;
  typedef typename Complex::Complex_vertex_iterator Complex_vertex_iterator;
  typedef typename Link::Vertex_handle Link_vertex_handle;

 private:
  const Complex* complex_;
  Complex_vertex_iterator current_vertex_;

  typedef Simplex_around_vertex_iterator<SkeletonBlockerComplex, Link> SAVI;
  SAVI current_simplex_around_current_vertex_;
  SAVI simplices_around_current_vertex_end_;

 public:
  Simplex_iterator() : complex_(0) { }

  Simplex_iterator(const Complex* complex) :
      complex_(complex),
      current_vertex_(complex->vertex_range().begin()),
      current_simplex_around_current_vertex_(complex, *current_vertex_),
      simplices_around_current_vertex_end_(complex, *current_vertex_, true) {
    // should not be called if the complex is empty
    assert(!complex->empty());
  }

 private:
  Simplex_iterator(const Complex* complex, bool end) :
      complex_(complex) {
    set_end();
  }

 public:
  Simplex_iterator(const Simplex_iterator& other)
      :
      complex_(other.complex_),
      current_vertex_(other.current_vertex_),
      current_simplex_around_current_vertex_(other.current_simplex_around_current_vertex_),
      simplices_around_current_vertex_end_(other.simplices_around_current_vertex_end_) { }

  friend Simplex_iterator make_begin_iterator(const Complex* complex) {
    if (complex->empty())
      return make_end_simplex_iterator(complex);
    else
      return Simplex_iterator(complex);
  }

  friend Simplex_iterator make_end_simplex_iterator(const Complex* complex) {
    return Simplex_iterator(complex, true);
  }

  bool equal(const Simplex_iterator& other) const {
    if (complex_ != other.complex_) return false;
    if (current_vertex_ != other.current_vertex_) return false;
    if (is_end() && other.is_end()) return true;
    if (current_simplex_around_current_vertex_ != other.current_simplex_around_current_vertex_)
      return false;
    return true;
  }

  void increment() {
    if (current_simplex_around_current_vertex_ != simplices_around_current_vertex_end_) {
      current_simplex_around_current_vertex_.increment();
      if (current_simplex_around_current_vertex_ == simplices_around_current_vertex_end_)
        goto_next_vertex();
    } else {
      goto_next_vertex();
    }
  }

  void goto_next_vertex() {
    current_vertex_.increment();
    if (!is_end()) {
      current_simplex_around_current_vertex_ = SAVI(complex_, *current_vertex_);
      simplices_around_current_vertex_end_ = SAVI(complex_, *current_vertex_, true);
    }
  }

  Simplex dereference() const {
    return current_simplex_around_current_vertex_.dereference();
  }

 private:
  void set_end() {
    current_vertex_ = complex_->vertex_range().end();
  }

  bool is_end() const {
    return (current_vertex_ == complex_->vertex_range().end());
  }
};

/**
 * Iterator through the maximal faces of the coboundary of a simplex.
 */
template<typename SkeletonBlockerComplex, typename Link>
class Simplex_coboundary_iterator :
public boost::iterator_facade < Simplex_coboundary_iterator<SkeletonBlockerComplex, Link>
, typename SkeletonBlockerComplex::Simplex, boost::forward_traversal_tag, typename SkeletonBlockerComplex::Simplex> {
  friend class boost::iterator_core_access;
  typedef SkeletonBlockerComplex Complex;
  typedef typename Complex::Vertex_handle Vertex_handle;
  typedef typename Complex::Edge_handle Edge_handle;
  typedef typename Complex::Simplex Simplex;
  typedef typename Complex::Complex_vertex_iterator Complex_vertex_iterator;

  // Link_vertex_handle == Complex_Vertex_handle but this renaming helps avoiding confusion
  typedef typename Link::Vertex_handle Link_vertex_handle;

 private:
  const Complex* complex;
  const Simplex& sigma;
  std::shared_ptr<Link> link;
  Complex_vertex_iterator current_vertex;
  Complex_vertex_iterator link_vertex_end;

 public:
  Simplex_coboundary_iterator() : complex(0) { }

  Simplex_coboundary_iterator(const Complex* complex_, const Simplex& sigma_) :
      complex(complex_),
      sigma(sigma_),
      // need only vertices of the link hence the true flag
      link(new Link(*complex_, sigma_, false, true)) {
    auto link_vertex_range = link->vertex_range();
    current_vertex = link_vertex_range.begin();
    link_vertex_end = link_vertex_range.end();
  }

  Simplex_coboundary_iterator(const Simplex_coboundary_iterator& other) :
      complex(other.complex),
      sigma(other.sigma),
      link(other.link),
      current_vertex(other.current_vertex),
      link_vertex_end(other.link_vertex_end) { }

  // returns an iterator to the end
  Simplex_coboundary_iterator(const Complex* complex_, const Simplex& sigma_, bool end) :
      complex(complex_),
      sigma(sigma_) {
    // to represent an end iterator without having to build a useless link, we use
    // the convection that link is not initialized.
  }

 private:
  Vertex_handle parent_vertex(Link_vertex_handle link_vh) const {
    return complex->convert_handle_from_another_complex(*link, link_vh);
  }

 public:
  friend std::ostream& operator<<(std::ostream& stream, const Simplex_coboundary_iterator& sci) {
    return stream;
  }

  // assume that iterator points to the same complex and comes from the same simplex
  bool equal(const Simplex_coboundary_iterator& other) const {
    assert(complex == other.complex && sigma == other.sigma);
    if (is_end()) return other.is_end();
    if (other.is_end()) return is_end();
    return *current_vertex == *(other.current_vertex);
  }

  void increment() {
    ++current_vertex;
  }

  Simplex dereference() const {
    Simplex res(sigma);
    res.add_vertex(parent_vertex(*current_vertex));
    return res;
  }

 private:
  bool is_end() const {
    return !link || current_vertex == link_vertex_end;
  }
};

}  // namespace skeleton_blocker

namespace skbl = skeleton_blocker;

}  // namespace Gudhi

#endif  // SKELETON_BLOCKER_ITERATORS_SKELETON_BLOCKERS_SIMPLICES_ITERATORS_H_
