/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef SIMPLEX_TREE_SIMPLEX_TREE_STAR_SIMPLEX_ITERATORS_H_
#define SIMPLEX_TREE_SIMPLEX_TREE_STAR_SIMPLEX_ITERATORS_H_

#include <gudhi/Debug_utils.h>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/version.hpp>
#include <boost/iterator/filter_iterator.hpp>

#include <vector>
#include <stdexcept>
#include <utility>  // for std::move
#include <functional>  // for std::greater
#include <algorithm>  // for std::includes

namespace Gudhi {

/** \private
 * \brief Iterator over all the roots of subtrees containing cofaces of all
 * dimension of a given simplex.
 *
 * Specifically, consider a simplex \f$\sigma\f$ whose vertices have maximal label
 * u. A maximal subtree of cofaces for \f$\sigma\f$ is a subtree of the simplex tree
 * rooted at a Node of label u, for which the root represents itself a coface of
 * \f$\sigma\f$.
 *
 * \details All Nodes of the simplex must store intrusive list hooks (option
 * link_nodes_by_label == true) in order to connect all nodes with a fixed label u
 * into a list.
 *
 * Forward iterator, value_type is SimplexTree::Simplex_handle.
 *
 * The implementation uses a filter iterator on all Nodes of label u, checking
 * whether they represent a coface of \f$\sigma\f$ or not.
 */
template <class SimplexTree>
class Simplex_tree_optimized_cofaces_rooted_subtrees_simplex_iterator
    : public boost::iterator_facade<Simplex_tree_optimized_cofaces_rooted_subtrees_simplex_iterator<SimplexTree>,
                                    typename SimplexTree::Simplex_handle const, boost::forward_traversal_tag> {
 public:
  using Simplex_handle = typename SimplexTree::Simplex_handle;
  using Siblings = typename SimplexTree::Siblings;
  using Vertex_handle = typename SimplexTree::Vertex_handle;
  using Node = typename SimplexTree::Node;
  using Static_vertex_vector = typename SimplexTree::Static_vertex_vector;

  /** \brief Predicate to check whether an input SimplexTree::Node represents a
   * coface of a simplex simp_, stored as a std::vector of
   * SimplexTree::Vertex_handle sorted in decreasing Vertex_handle order.
   *
   * \details Given a SimplexHandle in a simplex tree cpx_, traverses the tree
   * upwards
   * to find the sequence of Vertex_handle of simp_ as a subsequence of labels
   * encountered.
   * Does not test sh itself.
   * Used for filter_iterator in the optimized algorithm for
   * star_simplex_range.
   */
  class is_coface {
   public:
    is_coface() : cpx_(nullptr) {}
    is_coface(SimplexTree* cpx, Static_vertex_vector&& simp) : cpx_(cpx), simp_(simp) {}

    // Return true iff traversing the Node upwards to the root reads a
    // coface of simp_
    bool operator()(typename SimplexTree::Hooks_simplex_base& curr_hooks) {
      Node& curr_node = static_cast<Node&>(curr_hooks);
      Simplex_handle sh = cpx_->simplex_handle_from_node(curr_node);
      // first Node must always have label simp_.begin(); we assume it is true
      auto&& rng = cpx_->simplex_vertex_range(sh);
      auto rng_it = rng.begin();
      GUDHI_CHECK(*rng_it == simp_.front(), std::invalid_argument("first Node must always have label simp_.begin()"));
      auto simp_it = simp_.begin();
      // is simp_ a face of the simplex defined by sh ?
      return std::includes(++rng_it, rng.end(), ++simp_it, simp_.end(), std::greater<Vertex_handle>());
    }

   private:
    SimplexTree* cpx_;
    Static_vertex_vector simp_;  // vertices of simplex, in reverse order
  };

  typedef boost::filter_iterator<is_coface, typename SimplexTree::List_max_vertex::iterator>
      Filtered_cofaces_simplex_iterator;
  // any end() iterator
  Simplex_tree_optimized_cofaces_rooted_subtrees_simplex_iterator() : predicate_(), st_(nullptr) {}

  Simplex_tree_optimized_cofaces_rooted_subtrees_simplex_iterator(SimplexTree* cpx,
                                                                  Static_vertex_vector&& simp)
      : predicate_(cpx, std::move(simp)), st_(cpx) {
    GUDHI_CHECK(!simp.empty(), std::invalid_argument("cannot call for cofaces of an empty simplex"));
    auto list_ptr = st_->nodes_by_label(simp.front());
    GUDHI_CHECK(list_ptr != nullptr, std::runtime_error("invalid call to cofaces forest"));

    it_ = boost::make_filter_iterator(predicate_, list_ptr->begin(), list_ptr->end());
    end_ = boost::make_filter_iterator(predicate_, list_ptr->end(), list_ptr->end());
    Node& curr_node = static_cast<Node&>(*it_);
    sh_ = st_->simplex_handle_from_node(curr_node);
  }

 private:
  friend class boost::iterator_core_access;

  // valid when iterating along the SAME list of max vertex.
  bool equal(Simplex_tree_optimized_cofaces_rooted_subtrees_simplex_iterator const& other) const {
    if (other.st_ == nullptr) {
      return (st_ == nullptr);
    }
    if (st_ == nullptr) {
      return false;
    }
    return (it_ == other.it_);
  }

  Simplex_handle const& dereference() const { return sh_; }

  void increment() {
    if (++it_ == end_) {
      st_ = nullptr;
    }       //== end
    else {  // update sh_
      Node& curr_node = static_cast<Node&>(*it_);
      sh_ = st_->simplex_handle_from_node(curr_node);
    }
  }

  // given a Node of label max_v, returns true if the associated simplex is a coface of the simplex {..., max_v}. The
  // predicate stores the vertices of the simplex whose star we compute.
  is_coface predicate_;
  SimplexTree* st_;
  // filtered iterators over Nodes, filtered with predicate_
  Filtered_cofaces_simplex_iterator it_;
  Filtered_cofaces_simplex_iterator end_;
  // current Simplex_handle corresponding to Node pointed at by it_
  Simplex_handle sh_;
};

/** \brief Iterator over the simplices of the star of a simplex.
 *
 * \details All Nodes of the simplex must store intrusive list hooks (option
 * link_nodes_by_label == true) in order to connect all nodes with a fixed label u
 * into a list.
 *
 * Forward iterator, value_type is SimplexTree::Simplex_handle.
 *
 * The implementation uses a
 * Simplex_tree_optimized_cofaces_rooted_subtrees_simplex_iterator to iterate through all
 * roots of cofaces Nodes, and traverses each such subtree of the simplex tree.
 */
template <class SimplexTree>
class Simplex_tree_optimized_star_simplex_iterator
    : public boost::iterator_facade<Simplex_tree_optimized_star_simplex_iterator<SimplexTree>,
                                    typename SimplexTree::Simplex_handle const, boost::forward_traversal_tag> {
 public:
  using Simplex_handle = typename SimplexTree::Simplex_handle;
  using Siblings = typename SimplexTree::Siblings;
  using Vertex_handle = typename SimplexTree::Vertex_handle;
  using Node = typename SimplexTree::Node;
  using Static_vertex_vector = typename SimplexTree::Static_vertex_vector;

  // any end() iterator
  Simplex_tree_optimized_star_simplex_iterator() : st_(nullptr) {}

  Simplex_tree_optimized_star_simplex_iterator(SimplexTree* cpx, Static_vertex_vector&& simp)
      : st_(cpx), it_(cpx, std::move(simp)), end_(), sh_(*it_), sib_(st_->self_siblings(sh_)), children_stack_() {
    if (it_ == end_) {
      st_ = nullptr;
      return;
    }  // no coface subtree => end()
    is_root_ = true;
    sh_ = *it_;                      // sh_ is the root
    sib_ = st_->self_siblings(sh_);  // Siblings containing sh_
    if (st_->has_children(sh_)) {
      children_stack_.push_back(st_->children(sh_));
    }
    return;  // first root of coface subtree
  }

 private:
  friend class boost::iterator_core_access;

  // valid when iterating along the SAME list of max vertex.
  bool equal(Simplex_tree_optimized_star_simplex_iterator const& other) const {
    if (other.st_ == nullptr) {
      return (st_ == nullptr);
    }
    if (st_ == nullptr) {
      return false;
    }
    return (&(sh_->second) == &(other.sh_->second));
  }

  Simplex_handle const& dereference() const { return sh_; }

  /* Go to the next valid Simplex_handle for a coface, using a breadth first
   * search approach.
   *
   * Invariant:
   *
   * sh_ is a coface,
   * sib_ the Siblings containing sh_,
   * it_ the root of a coface subtree, such that dim_root_ <= exact_dim_cofaces_.
   * bfs_queue contains Siblings inside the coface subtree.
   *
   * Additionally,
   *
   *  - computing all cofaces: sh_ points to a coface of any dimension. bfs_queue
   * contains a collection of Siblings* that must be considered (as well as there
   * children). These are all sets of children of simplices in
   * [ sib_->members().begin(), sh_ ]
   */
  void increment_all_cofaces() {
    ++sh_;  // next sibling
    // if no more sibling or root of coface tree, go down or to next subtree
    if (is_root_ || sh_ == sib_->members().end()) {
      is_root_ = false;
      if (!children_stack_.empty()) {
        sib_ = children_stack_.back();
        children_stack_.pop_back();
        sh_ = sib_->members().begin();  // don't track dimensions
        if (st_->has_children(sh_)) {
          children_stack_.push_back(st_->children(sh_));
        }
      } else {  // bfs_queue == empty, go to root of next coface subtree
        if (++it_ == end_) {
          st_ = nullptr;
          return;
        }  // no more subtree => end()
        is_root_ = true;
        sh_ = *it_;                      // sh_ is the root
        sib_ = st_->self_siblings(sh_);  // Siblings containing sh_
        if (st_->has_children(sh_)) {
          children_stack_.push_back(st_->children(sh_));
        }
        return;  // next root of coface
      }
    } else {  // sh_ is valid, simply add its children to the queue
      if (st_->has_children(sh_)) {
        children_stack_.push_back(st_->children(sh_));
      }
    }
  }

  // For cofaces, enumerating all the star and testing which simplices have the right dimension is suboptimal.
  // This could be optimized later.
  void increment() { increment_all_cofaces(); }

  // Let s be the simplex in a complex C whose star is
  // iterated through. Let max_v denote the maximal label of vertices in s.
  SimplexTree* st_;  // Simplex tree for complex C
  // The cofaces of s form a subforest of the simplex tree. The roots of trees in this
  // forest have label max_v.
  //[it_,end_) == range of Simplex_handles of the roots of the cofaces trees (any dim)
  Simplex_tree_optimized_cofaces_rooted_subtrees_simplex_iterator<SimplexTree> it_;
  Simplex_tree_optimized_cofaces_rooted_subtrees_simplex_iterator<SimplexTree> end_;
  // curr Simplex_handle, returned by operator*, pointing to a coface of s
  Simplex_handle sh_;
  // set of siblings containing sh_ in the Simplex_tree
  Siblings* sib_;  //
  // Save children in a list to avoid calling sib_->members().find(.)
  std::vector<Siblings*> children_stack_;
  // true iff sh_ points to the root of a coface subtree
  bool is_root_;
};

}  // namespace Gudhi

#endif  // SIMPLEX_TREE_SIMPLEX_TREE_STAR_SIMPLEX_ITERATORS_H_
