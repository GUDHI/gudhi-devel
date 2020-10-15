/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - 2020/09 Clément Maria: Add the class 
 *                               Simplex_tree_optimized_cofaces_simplex_iterator
 */

#ifndef SIMPLEX_TREE_COFACES_ITERATORS_
#define SIMPLEX_TREE_COFACES_ITERATORS_

#include <gudhi/Debug_utils.h>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/version.hpp>
#include <boost/container/static_vector.hpp>
#include <boost/iterator/filter_iterator.hpp>

#include <vector>

namespace Gudhi {

/* \addtogroup simplex_tree
 * Iterators and range types for the Simplex_tree.
 * @{
 */

/* \brief Iterator over all the roots of subtrees containing cofaces of all 
 * dimension of a given simplex. If the simplex as vertex of maximal label u, all
 * roots have label u.
 *
 * Forward iterator, value_type is SimplexTree::Simplex_handle.*/
template <class SimplexTree>
class Simplex_tree_opt_cofaces_rooted_subtrees_simplex_iterator
	: public boost::iterator_facade<Simplex_tree_opt_cofaces_rooted_subtrees_simplex_iterator<SimplexTree>,
	typename SimplexTree::Simplex_handle const, boost::forward_traversal_tag> {
public:
    typedef typename SimplexTree::Simplex_handle   Simplex_handle;
    typedef typename SimplexTree::Siblings         Siblings;
    typedef typename SimplexTree::Vertex_handle    Vertex_handle;
    typedef typename SimplexTree::Node             Node;

  /* Predicate to check whether an input SimplexTree::Node represents a
   * coface of a simplex simp_,
   * stored as a std::vector of SimplexTree::Vertex_handle.
   *
   * Given a SimplexHandle in a simplex tree cpx_, traverse the tree upwards
   * to find the sequence of Vertex_handle of simp_ as a subsequence of labels 
   * encountered. 
   * Does not test sh itself.
   * Used for filter_iterator in the optimized algorithm for
   * cofaces_simplex_range.
   */
  class is_coface {
    public:
    	is_coface() : cpx_(nullptr) {}
    	is_coface(SimplexTree* cpx, const std::vector<Vertex_handle> &simp) 
      : cpx_(cpx), simp_(simp) {}
    	// Return true iff traversing the Node upwards to the root reads a
    	// coface of simp_ of codimension codim_
    	bool operator()(typename SimplexTree::Hooks_simplex_base &curr_hooks) {
  	    Node& curr_node = static_cast<Node&>(curr_hooks);
  	    auto vertex_it = simp_.begin();//largest label is first
  	    // first Node must always have label simp_.begin(); we assume it is true
  	    auto curr_sib = cpx_->self_siblings(curr_node, *vertex_it);
  	    if (++vertex_it == simp_.end()) { return true; }
  	    while (curr_sib->oncles() != nullptr) {  
      		if (curr_sib->parent() == *vertex_it) {
    		    if (++vertex_it == simp_.end()) { return true; }//we found a coface
      		}
      		curr_sib = curr_sib->oncles();
    	  }
  	    return false;
    	}

    private:
    	SimplexTree                * cpx_;
    	std::vector<Vertex_handle>   simp_;//vertices of simplex, reverse order
  };

  typedef boost::filter_iterator<is_coface, typename SimplexTree::List_max_vertex::iterator>                                    Filtered_cofaces_simplex_iterator;
  // any end() iterator
  Simplex_tree_opt_cofaces_rooted_subtrees_simplex_iterator() 
  : predicate_(), st_(nullptr) {}

  Simplex_tree_opt_cofaces_rooted_subtrees_simplex_iterator(SimplexTree* cpx, 
                                             const std::vector<Vertex_handle>& simp)
	: predicate_(cpx, simp), st_(cpx) 
  {
    GUDHI_CHECK(!simp.empty(), "cannot call for cofaces of an empty simplex");
  	max_v_ = *(simp.begin());
  	auto list_ptr = st_->nodes_by_label(max_v_);
    GUDHI_CHECK(list_ptr != nullptr, "invalid call to cofaces forest");

  	it_ = boost::make_filter_iterator(predicate_, 
                                      list_ptr->begin(), list_ptr->end());
  	end_ = boost::make_filter_iterator(predicate_, 
                                       list_ptr->end(), list_ptr->end());
  	Node& curr_node = static_cast<Node&>(*it_);
  	auto curr_sib = st_->self_siblings(curr_node, max_v_);
  	sh_ = curr_sib->find(max_v_);
  }

private:
  friend class boost::iterator_core_access;

  // valid when iterating along the SAME list of max vertex.
  bool equal(Simplex_tree_opt_cofaces_rooted_subtrees_simplex_iterator const& other)
  const {
  	if (other.st_ == nullptr) { return (st_ == nullptr); }
  	if (st_ == nullptr) { return false; }
  	return (it_ == other.it_);//(&(*it_) == &(*(other.it_)));
  }

  Simplex_handle const& dereference() const { return sh_; }

  void increment() {
  	if (++it_ == end_) { st_ = nullptr; } //== end
  	else {  // update sh_
      Node& curr_node = static_cast<Node&>(*it_);
      auto curr_sib = st_->self_siblings(curr_node, max_v_);
      sh_ = curr_sib->find(max_v_);
  	}
  }

//given a Node of label max_v, returns true if the associated simplex is a coface of the simplex
  is_coface                         predicate_;
  SimplexTree                     * st_;
//filtered iterators over Nodes of same label max_v_, filtered with predicate_
  Filtered_cofaces_simplex_iterator it_;
  Filtered_cofaces_simplex_iterator end_;
//max label of the simplex whose cofaces are computed
  Vertex_handle                     max_v_;
// current Simplex_handle corresponding to Node pointed at by it_
  Simplex_handle                    sh_;  
};

/* \brief Iterator over all cofaces of dimension exactly exact_dim_cofaces_ for
 * a simplex.
 *
 * Uses hooks stored in the Node of the SimplexTree.*/
template <class SimplexTree>
class Simplex_tree_optimized_cofaces_simplex_iterator
	: public boost::iterator_facade<Simplex_tree_optimized_cofaces_simplex_iterator<SimplexTree>,
	typename SimplexTree::Simplex_handle const, boost::forward_traversal_tag> {
public:
  typedef typename SimplexTree::Simplex_handle Simplex_handle;
  typedef typename SimplexTree::Siblings       Siblings;
  typedef typename SimplexTree::Vertex_handle  Vertex_handle;
  typedef typename SimplexTree::Node           Node;

  // any end() iterator
  Simplex_tree_optimized_cofaces_simplex_iterator() : st_(nullptr) {}

  Simplex_tree_optimized_cofaces_simplex_iterator(SimplexTree* cpx, 
                                            const std::vector<Vertex_handle>& simp)
  : st_(cpx),
    it_(cpx, simp),
    end_(),
    sh_(*it_),
    sib_(st_->self_siblings(sh_)),
    bfs_queue_() 
  { 
    if(it_ == end_) { st_ = nullptr; return; }//no coface subtree => end()
    is_root_  = true;
    sh_       = *it_;  //sh_ is the root
    sib_      = st_->self_siblings(sh_);//Siblings containing sh_
    if(st_->has_children(sh_)) { bfs_queue_.push(st_->children(sh_)); }
    return;//first root of coface subtree
  }

private:
  friend class boost::iterator_core_access;

  // valid when iterating along the SAME list of max vertex.
  bool equal(Simplex_tree_optimized_cofaces_simplex_iterator const& other) const {
  	if (other.st_ == nullptr) { return (st_ == nullptr); }
  	if (st_ == nullptr) { return false; }
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
 *  - if all_cofaces == true: sh_ points to a coface of any dimension. bfs_queue 
 * contains a collection of Siblings* that must be considered (as well as there 
 * children). These are all sets of children of simplices in 
 * [ sib_->members().begin(), sh_ ]
 */
  void increment_all_cofaces() {
    ++sh_; //next sibling
    //if no more sibling or root of coface tree, go down or to next subtree
    if( is_root_ || sh_ == sib_->members().end() ) {
      is_root_ = false;
      if(!bfs_queue_.empty()) {
        sib_ = bfs_queue_.front();  bfs_queue_.pop();
        sh_  = sib_->members().begin();//don't track dimensions
        if(st_->has_children(sh_)) { bfs_queue_.push(st_->children(sh_)); }
      }
      else {//bfs_queue == empty, go to root of next coface subtree
        if (++it_ == end_) { st_ = nullptr; return; }//no more subtree => end()
        is_root_  = true;
        sh_       = *it_;  //sh_ is the root
        sib_      = st_->self_siblings(sh_);//Siblings containing sh_
        if(st_->has_children(sh_)) { bfs_queue_.push(st_->children(sh_)); }
        return;//next root of coface
      } 
    }
    else{ //sh_ is valid, simply add its children to the queue
      if(st_->has_children(sh_)) { bfs_queue_.push(st_->children(sh_)); }
    }
  }

  void increment() { increment_all_cofaces(); }

//Let s be the simplex in a complex C whose cofaces (of a certain codimension) are 
//iterated through. Let max_v denote the maximal label of vertices in s.
  SimplexTree          * st_;//Simplex tree for complex C
//The cofaces of s form a subforest of the simplex tree. The roots of trees in this
//forest have label max_v.
//[it_,end_) == range of Simplex_handles of the roots of the cofaces trees (any dim)
  Simplex_tree_opt_cofaces_rooted_subtrees_simplex_iterator<SimplexTree> it_;
  Simplex_tree_opt_cofaces_rooted_subtrees_simplex_iterator<SimplexTree> end_;
//curr Simplex_handle, returned by operator*. A coface of s of appropriate codim  
  Simplex_handle         sh_;               
//sh_ is always in a coface subtree
  // int            dim_root_; //dimension of the root       
  // int            dim_sh_;   //dimension of simplex sh_, inside the subtree
 //exact dimension of the cofaces computed, max_int otherwise if we want ALL cofaces
  // int            exact_dim_cofaces_;
  Siblings             * sib_;//
  // bool           all_cofaces;//true if we compute all cofaces, false otherwise
  //use a bfs search to avoid calling sib_->members().find(.)
  std::queue<Siblings *> bfs_queue_;
  bool                   is_root_;
};

/* @} */  // end addtogroup simplex_tree
}  // namespace Gudhi

#endif  // SIMPLEX_TREE_COFACES_ITERATORS_
