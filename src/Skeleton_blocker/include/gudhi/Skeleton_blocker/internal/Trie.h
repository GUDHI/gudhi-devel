/* This file is part of the Gudhi Library. The Gudhi library 
 *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
 *    library for computational topology.
 *
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014  INRIA Sophia Antipolis-Méditerranée (France)
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
 * 
 */


#ifndef SKELETON_BLOCKER_INTERNAL_TRIE_H_
#define SKELETON_BLOCKER_INTERNAL_TRIE_H_

#include <memory>
#include <vector>
#include <deque>
#include <set>

namespace Gudhi {

namespace skbl {

template<typename SimplexHandle>
struct Trie {
  typedef SimplexHandle Simplex_handle;
  typedef typename SimplexHandle::Vertex_handle Vertex_handle;

  Vertex_handle v;
  std::vector<std::shared_ptr<Trie> > childs;
  // std::vector<std::unique_ptr<Trie> > childs; -> use of deleted function
 private:
  const Trie* parent_;

 public:
  Trie() : parent_(0) { }

  Trie(Vertex_handle v_) : v(v_), parent_(0) { }

  Trie(Vertex_handle v_, Trie* parent) : v(v_), parent_(parent) { }

  bool operator==(const Trie& other) const {
    return (v == other.v);
  }

  void add_child(Trie* child) {
    if (child) {
      std::shared_ptr<Trie> ptr_to_add(child);
      childs.push_back(ptr_to_add);
      child->parent_ = this;
    }
  }

  typedef typename Simplex_handle::Simplex_vertex_const_iterator Simplex_vertex_const_iterator;

  Trie* make_trie(Simplex_vertex_const_iterator s_it, Simplex_vertex_const_iterator s_end) {
    if (s_it == s_end) {
      return 0;
    } else {
      Trie* res = new Trie(*s_it);
      Trie* child = make_trie(++s_it, s_end);
      res->add_child(child);
      return res;
    }
  }

 private:
  // go down recursively in the tree while advancing the simplex iterator.
  // when it reaches a leaf, it inserts the remaining that is not present
  void add_simplex_helper(Simplex_vertex_const_iterator s_it, Simplex_vertex_const_iterator s_end) {
    assert(*s_it == v);
    ++s_it;
    if (s_it == s_end) return;
    if (!is_leaf()) {
      for (auto child : childs) {
        if (child->v == *s_it)
          return child->add_simplex_helper(s_it, s_end);
      }
      // s_it is not found and needs to be inserted
    }
    // not leaf -> remaining of s needs to be inserted
    Trie * son_with_what_remains_of_s(make_trie(s_it, s_end));
    add_child(son_with_what_remains_of_s);
    return;
  }

  void maximal_faces_helper(std::vector<Simplex_handle>& res) const {
    if (is_leaf()) res.push_back(simplex());
    else
      for (auto child : childs)
        child->maximal_faces_helper(res);
  }

 public:
  /**
   * adds the simplex to the trie
   */
  void add_simplex(const Simplex_handle& s) {
    if (s.empty()) return;
    assert(v == s.first_vertex());
    add_simplex_helper(s.begin(), s.end());
  }

  std::vector<Simplex_handle> maximal_faces() const {
    std::vector<Simplex_handle> res;
    maximal_faces_helper(res);
    return res;
  }

  /**
   * Goes to the root in the trie to consitute simplex
   */
  void add_vertices_up_to_the_root(Simplex_handle& res) const {
    res.add_vertex(v);
    if (parent_)
      parent_->add_vertices_up_to_the_root(res);
  }

  Simplex_handle simplex() const {
    Simplex_handle res;
    add_vertices_up_to_the_root(res);
    return res;
  }

  bool is_leaf() const {
    return childs.empty();
  }

  bool is_root() const {
    return parent_ == 0;
  }

  const Trie* parent() {
    return parent_;
  }

  void remove_leaf() {
    assert(is_leaf);
    if (!is_root())
      parent_->childs.erase(this);
  }

  /**
   * true iff the simplex corresponds to one node in the trie
   */
  bool contains(const Simplex_handle& s) const {
    Trie const* current = this;
    if (s.empty()) return true;
    if (current->v != s.first_vertex()) return false;
    auto s_pos = s.begin();
    ++s_pos;
    while (s_pos != s.end() && current != 0) {
      bool found = false;
      for (const auto child : current->childs) {
        if (child->v == *s_pos) {
          ++s_pos;
          current = child.get();
          found = true;
          break;
        }
      }
      if (!found) return false;
    }
    return current != 0;
  }

  Trie* go_bottom_left() {
    if (is_leaf())
      return this;
    else
      return (*childs.begin())->go_bottom_left();
  }

  friend std::ostream& operator<<(std::ostream& stream, const Trie& trie) {
    stream << "T( " << trie.v << " ";
    for (auto t : trie.childs)
      stream << *t;
    stream << ")";
    return stream;
  }
};

template<typename SimplexHandle>
struct Tries {
  typedef typename SimplexHandle::Vertex_handle Vertex_handle;
  typedef SimplexHandle Simplex_handle;

  typedef Trie<Simplex_handle> STrie;

  template<typename SimpleHandleOutputIterator>
  Tries(unsigned num_vertices, SimpleHandleOutputIterator simplex_begin, SimpleHandleOutputIterator simplex_end) :
      cofaces_(num_vertices, 0) {
    for (auto i = 0u; i < num_vertices; ++i)
      cofaces_[i] = new STrie(Vertex_handle(i));
    for (auto s_it = simplex_begin; s_it != simplex_end; ++s_it) {
      if (s_it->dimension() >= 1)
        cofaces_[s_it->first_vertex()]->add_simplex(*s_it);
    }
  }

  ~Tries() {
    for (STrie* t : cofaces_)
      delete t;
  }

  // return a simplex that consists in all u such uv is an edge and u>v

  Simplex_handle positive_neighbors(Vertex_handle v) const {
    Simplex_handle res;
    for (auto child : cofaces_[v]->childs)
      res.add_vertex(child->v);
    return res;
  }

  bool contains(const Simplex_handle& s) const {
    auto first_v = s.first_vertex();
    return cofaces_[first_v]->contains(s);
  }

  friend std::ostream& operator<<(std::ostream& stream, const Tries& tries) {
    for (auto trie : tries.cofaces_)
      stream << *trie << std::endl;
    return stream;
  }

  // init_next_dimension must be called first

  std::vector<Simplex_handle> next_dimension_simplices() const {
    std::vector<Simplex_handle> res;
    while (!to_see_.empty() && to_see_.front()->simplex().dimension() == current_dimension_) {
      res.emplace_back(to_see_.front()->simplex());
      for (auto child : to_see_.front()->childs)
        to_see_.push_back(child.get());
      to_see_.pop_front();
    }
    ++current_dimension_;
    return res;
  }

  void init_next_dimension() const {
    for (auto trie : cofaces_)
      to_see_.push_back(trie);
  }

 private:
  mutable std::deque<STrie*> to_see_;
  mutable unsigned current_dimension_ = 0;
  std::vector<STrie*> cofaces_;
};

}  // namespace skbl

}  // namespace Gudhi

#endif  // SKELETON_BLOCKER_INTERNAL_TRIE_H_
