/*    This file is a prototype for the Gudhi Library.
 *    Author(s):       Clément Maria
 *    Copyright (C) 2021 Inria
 *    This version is under developement, please do not redistribute this software. 
 *    This program is for academic research use only. 
 */

#ifndef SIMPLEX_TREE_SIMPLEX_TREE_MORSE_ITERATORS_H_
#define SIMPLEX_TREE_SIMPLEX_TREE_MORSE_ITERATORS_H_

#include <boost/iterator/iterator_facade.hpp>
#include <gudhi/Discrete_morse_theory.h>

/*---------------------------------------------------------------------------*/
/* \brief Iterator over the simplices of the boundary of a
 *  critical simplex in a Morse complex.
 *
 * Forward iterator, value_type is SimplexTree::Simplex_handle.*/
template<class SimplexTree>
class Simplex_tree_morse_boundary_simplex_iterator : public boost::iterator_facade<
    Simplex_tree_morse_boundary_simplex_iterator<SimplexTree>,
    typename SimplexTree::Simplex_handle const, boost::forward_traversal_tag> {
 public:
  typedef typename SimplexTree::Simplex_handle Simplex_handle;
  typedef typename SimplexTree::Vertex_handle Vertex_handle;
  typedef typename SimplexTree::Siblings Siblings;

// any end() iterator
  explicit Simplex_tree_morse_boundary_simplex_iterator(SimplexTree * st)
      : st_(st), sh_it_(boundary_.end()), end_(true) {}

  template<class SimplexHandle>
  Simplex_tree_morse_boundary_simplex_iterator(SimplexTree * st, SimplexHandle sh)
    : st_(st) 
  {
    auto boundary_map = Gudhi::dmt::boundary_morse_complex(st_, sh);
    boundary_(boundary_map.begin(),boundary_map.end());
    sh_it_ = boundary_.begin();
    if(sh_it_ == boundary_.end()) { end_ = true; }
    else { end_ = false; }
  }

/** \brief Return the coefficient associated to the incidence 
 * cell - cell in the boundary.
 * 
 * \details For a simplicial complex, takes value +1 or -1.
 */ 
  int coefficient() { return sh_it_->second; }

 private:
  friend class boost::iterator_core_access;
// valid when iterating along the SAME boundary.
  bool equal(Simplex_tree_morse_boundary_simplex_iterator const& other) const {
    if(end_) { return end_ == other.end_; }
    return sh_it_ == other.sh_it_;
  }

  Simplex_handle const& dereference() const {
    return sh_it_->first;
  }

  void increment() {
    if(!end_) {
      ++sh_it_;
      if(sh_it_ == boundary_.end()) { end_ = true; }
    }
  }

  SimplexTree  * st_;//the simplex tree we are working in
  //precomputed set of simplex handles for the boundary
  std::vector< std::pair<Simplex_handle,int> > boundary_;
  //iterator in boundry_ pointing to the current simplex handle
  typename std::vector< Simplex_handle >::iterator sh_it_;
  bool                                         end_;//true iff the iterator it end()

};

#endif //SIMPLEX_TREE_SIMPLEX_TREE_MORSE_ITERATORS_H_