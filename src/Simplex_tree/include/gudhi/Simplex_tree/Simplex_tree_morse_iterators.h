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
      // : last_(st->null_vertex()),
      //   next_(st->null_vertex()),
      //   sib_(nullptr),
      //   sh_(st->null_simplex()),
      //   st_(st)  {
  }

  template<class SimplexHandle>
  Simplex_tree_morse_boundary_simplex_iterator(SimplexTree * st, SimplexHandle sh)
      // : last_(sh->first),
      //   next_(st->null_vertex()),
      //   sib_(nullptr),
      //   sh_(st->null_simplex()),
      //   st_(st) {
    // // Only check once at the beginning instead of for every increment, as this is expensive.
    // if (SimplexTree::Options::contiguous_vertices)
    //   GUDHI_CHECK(st_->contiguous_vertices(), "The set of vertices is not { 0, ..., n } without holes");
    // Siblings * sib = st->self_siblings(sh);
    // next_ = sib->parent();
    // sib_ = sib->oncles();
    // if (sib_ != nullptr) {
    //   if constexpr(SimplexTree::Options::contiguous_vertices) {
    //     if(sib_->oncles() == nullptr) { sh_ = sib_->members_.begin()+next_; }
    //     else { sh_ = sib_->find(next_); }
    //   }
    //   else { sh_ = sib_->find(next_); }
    // }
  }

 private:
  friend class boost::iterator_core_access;
// valid when iterating along the SAME boundary.
  bool equal(Simplex_tree_morse_boundary_simplex_iterator const& other) const {
    // return sh_ == other.sh_;
  }

  Simplex_handle const& dereference() const {
    // assert(sh_ != st_->null_simplex());
    // return sh_;
  }

  void increment() {
    // if (sib_ == nullptr) {
    //   sh_ = st_->null_simplex();
    //   return;
    // }

    // Siblings * for_sib = sib_;
    // Siblings * new_sib = sib_->oncles();
    // auto rit = suffix_.rbegin();
    // if (new_sib == nullptr) {
    //   // We reached the root, use a short-cut to find a vertex.
    //   if (rit == suffix_.rend()) {
    //     // Segment, this vertex is the last boundary simplex
    //     if constexpr(SimplexTree::Options::contiguous_vertices) {
    //       sh_ = for_sib->members_.begin()+last_;
    //     }
    //     else { sh_ = for_sib->members_.find(last_); }
    //     sib_ = nullptr;
    //     return;
    //   } else {
    //     // Dim >= 2, initial step of the descent
    //     if constexpr(SimplexTree::Options::contiguous_vertices) {
    //       sh_ = for_sib->members_.begin()+*rit;
    //     }
    //     else { sh_ = for_sib->members_.find(*rit); }
    //     for_sib = sh_->second.children();
    //     ++rit;
    //   }
    // }
    // for (; rit != suffix_.rend(); ++rit) {
    //   sh_ = for_sib->find(*rit);
    //   for_sib = sh_->second.children();
    // }
    // sh_ = for_sib->find(last_);  // sh_ points to the right simplex now
    // suffix_.push_back(next_);
    // next_ = sib_->parent();
    // sib_ = new_sib;
  }

  // Most of the storage should be moved to the range, iterators should be light.
  // Vertex_handle last_;  // last vertex of the simplex
  // Vertex_handle next_;  // next vertex to push in suffix_
  // // 40 seems a conservative bound on the dimension of a Simplex_tree for now,
  // // as it would not fit on the biggest hard-drive.
  // boost::container::static_vector<Vertex_handle, 40> suffix_;
  // // static_vector still has some overhead compared to a trivial hand-made
  // // version using std::aligned_storage, or compared to making suffix_ static.
  // Siblings * sib_;  // where the next search will start from
  // Simplex_handle sh_;  // current Simplex_handle in the boundary
  // SimplexTree * st_;  // simplex containing the simplicial complex
};

#endif //SIMPLEX_TREE_SIMPLEX_TREE_MORSE_ITERATORS_H_