/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef SIMPLEX_TREE_SIMPLEX_TREE_ITERATORS_H_
#define SIMPLEX_TREE_SIMPLEX_TREE_ITERATORS_H_

#include <gudhi/Debug_utils.h>

#include <boost/iterator/iterator_facade.hpp>
#include <boost/version.hpp>
#include <boost/container/static_vector.hpp>

#include <vector>

namespace Gudhi {

/* \addtogroup simplex_tree
 * Iterators and range types for the Simplex_tree.
 * @{
 */

/* \brief Iterator over the vertices of a simplex
 * in a SimplexTree.
 *
 * Forward iterator, 'value_type' is SimplexTree::Vertex_handle.*/
template<class SimplexTree>
class Simplex_tree_simplex_vertex_iterator : public boost::iterator_facade<
    Simplex_tree_simplex_vertex_iterator<SimplexTree>,
    typename SimplexTree::Vertex_handle const, boost::forward_traversal_tag,
    typename SimplexTree::Vertex_handle const> {
 public:
  typedef typename SimplexTree::Simplex_handle Simplex_handle;
  typedef typename SimplexTree::Siblings Siblings;
  typedef typename SimplexTree::Vertex_handle Vertex_handle;

  explicit Simplex_tree_simplex_vertex_iterator(SimplexTree const* st)
      :  // any end() iterator
        sib_(nullptr),
        v_(st->null_vertex()) {
  }

  Simplex_tree_simplex_vertex_iterator(SimplexTree const* st, Simplex_handle sh)
      : sib_(st->self_siblings(sh)),
        v_(sh->first) {
  }

 private:
  friend class boost::iterator_core_access;

  bool equal(Simplex_tree_simplex_vertex_iterator const &other) const {
    return sib_ == other.sib_ && v_ == other.v_;
  }

  Vertex_handle const& dereference() const {
    return v_;
  }

  void increment() {
    v_ = sib_->parent();
    sib_ = sib_->oncles();
  }

  Siblings * sib_;
  Vertex_handle v_;
};

/*---------------------------------------------------------------------------*/
/* \brief Iterator over the simplices of the boundary of a
 *  simplex.
 *
 * Forward iterator, value_type is SimplexTree::Simplex_handle.*/
template<class SimplexTree>
class Simplex_tree_boundary_simplex_iterator : public boost::iterator_facade<
    Simplex_tree_boundary_simplex_iterator<SimplexTree>,
    typename SimplexTree::Simplex_handle const, boost::forward_traversal_tag> {
 public:
  typedef typename SimplexTree::Simplex_handle Simplex_handle;
  typedef typename SimplexTree::Vertex_handle Vertex_handle;
  typedef typename SimplexTree::Siblings Siblings;

  // For cython purpose only. The object it initializes should be overwritten ASAP and never used before it is overwritten.
  Simplex_tree_boundary_simplex_iterator()
      : sib_(nullptr),
        st_(nullptr)  {
  }

// any end() iterator
  explicit Simplex_tree_boundary_simplex_iterator(SimplexTree * st)
      : last_(st->null_vertex()),
        next_(st->null_vertex()),
        sib_(nullptr),
        sh_(st->null_simplex()),
        st_(st)  {
  }

  template<class SimplexHandle>
  Simplex_tree_boundary_simplex_iterator(SimplexTree * st, SimplexHandle sh)
      : last_(sh->first),
        next_(st->null_vertex()),
        sib_(nullptr),
        sh_(st->null_simplex()),
        st_(st) {
    // Only check once at the beginning instead of for every increment, as this is expensive.
    if (SimplexTree::Options::contiguous_vertices)
      GUDHI_CHECK(st_->contiguous_vertices(), "The set of vertices is not { 0, ..., n } without holes");
    Siblings * sib = st->self_siblings(sh);
    next_ = sib->parent();
    sib_ = sib->oncles();
    if (sib_ != nullptr) {
      if (SimplexTree::Options::contiguous_vertices && sib_->oncles() == nullptr)
        // Only relevant for edges
        sh_ = sib_->members_.begin()+next_;
      else
        sh_ = sib_->find(next_);
    }
  }

 private:
  friend class boost::iterator_core_access;
// valid when iterating along the SAME boundary.
  bool equal(Simplex_tree_boundary_simplex_iterator const& other) const {
    return sh_ == other.sh_;
  }

  Simplex_handle const& dereference() const {
    assert(sh_ != st_->null_simplex());
    return sh_;
  }

  void increment() {
    if (sib_ == nullptr) {
      sh_ = st_->null_simplex();
      return;
    }

    Siblings * for_sib = sib_;
    Siblings * new_sib = sib_->oncles();
    auto rit = suffix_.rbegin();
    if (SimplexTree::Options::contiguous_vertices && new_sib == nullptr) {
      // We reached the root, use a short-cut to find a vertex.
      if (rit == suffix_.rend()) {
        // Segment, this vertex is the last boundary simplex
        sh_ = for_sib->members_.begin()+last_;
        sib_ = nullptr;
        return;
      } else {
        // Dim >= 2, initial step of the descent
        sh_ = for_sib->members_.begin()+*rit;
        for_sib = sh_->second.children();
        ++rit;
      }
    }
    for (; rit != suffix_.rend(); ++rit) {
      sh_ = for_sib->find(*rit);
      for_sib = sh_->second.children();
    }
    sh_ = for_sib->find(last_);  // sh_ points to the right simplex now
    suffix_.push_back(next_);
    next_ = sib_->parent();
    sib_ = new_sib;
  }

  // Most of the storage should be moved to the range, iterators should be light.
  Vertex_handle last_;  // last vertex of the simplex
  Vertex_handle next_;  // next vertex to push in suffix_
  // 40 seems a conservative bound on the dimension of a Simplex_tree for now,
  // as it would not fit on the biggest hard-drive.
  boost::container::static_vector<Vertex_handle, 40> suffix_;
  // static_vector still has some overhead compared to a trivial hand-made
  // version using std::aligned_storage, or compared to making suffix_ static.
  Siblings * sib_;  // where the next search will start from
  Simplex_handle sh_;  // current Simplex_handle in the boundary
  SimplexTree * st_;  // simplex containing the simplicial complex
};
/*---------------------------------------------------------------------------*/
/* \brief Iterator over the simplices of a simplicial complex.
 *
 * Forward iterator, value_type is SimplexTree::Simplex_handle.*/
template<class SimplexTree>
class Simplex_tree_complex_simplex_iterator : public boost::iterator_facade<
    Simplex_tree_complex_simplex_iterator<SimplexTree>,
    typename SimplexTree::Simplex_handle const, boost::forward_traversal_tag> {
 public:
  typedef typename SimplexTree::Simplex_handle Simplex_handle;
  typedef typename SimplexTree::Siblings Siblings;
  typedef typename SimplexTree::Vertex_handle Vertex_handle;

// any end() iterator
  Simplex_tree_complex_simplex_iterator()
      : sib_(nullptr),
        st_(nullptr) {
  }

  explicit Simplex_tree_complex_simplex_iterator(SimplexTree * st)
      : sib_(nullptr),
        st_(st) {
    if (st == nullptr || st->root() == nullptr || st->root()->members().empty()) {
      st_ = nullptr;
    } else {
      sh_ = st->root()->members().begin();
      sib_ = st->root();
      while (st->has_children(sh_)) {
        sib_ = sh_->second.children();
        sh_ = sib_->members().begin();
      }
    }
  }
 private:
  friend class boost::iterator_core_access;

// valid when iterating along the SAME boundary.
  bool equal(Simplex_tree_complex_simplex_iterator const& other) const {
    if (other.st_ == nullptr) {
      return (st_ == nullptr);
    }
    if (st_ == nullptr) {
      return false;
    }
    return (&(sh_->second) == &(other.sh_->second));
  }

  Simplex_handle const& dereference() const {
    return sh_;
  }

// Depth first traversal.
  void increment() {
    ++sh_;
    if (sh_ == sib_->members().end()) {
      if (sib_->oncles() == nullptr) {
        st_ = nullptr;
        return;
      }  // reach the end
      sh_ = sib_->oncles()->members().find(sib_->parent());
      sib_ = sib_->oncles();
      return;
    }
    while (st_->has_children(sh_)) {
      sib_ = sh_->second.children();
      sh_ = sib_->members().begin();
    }
  }

  Simplex_handle sh_;
  Siblings * sib_;
  SimplexTree * st_;
};

/* \brief Iterator over the simplices of the skeleton of a given
 * dimension of the simplicial complex.
 *
 * Forward iterator, value_type is SimplexTree::Simplex_handle.*/
template<class SimplexTree>
class Simplex_tree_skeleton_simplex_iterator : public boost::iterator_facade<
    Simplex_tree_skeleton_simplex_iterator<SimplexTree>,
    typename SimplexTree::Simplex_handle const, boost::forward_traversal_tag> {
 public:
  typedef typename SimplexTree::Simplex_handle Simplex_handle;
  typedef typename SimplexTree::Siblings Siblings;
  typedef typename SimplexTree::Vertex_handle Vertex_handle;

// any end() iterator
  Simplex_tree_skeleton_simplex_iterator()
      : sib_(nullptr),
        st_(nullptr),
        dim_skel_(0),
        curr_dim_(0) {
  }

  Simplex_tree_skeleton_simplex_iterator(SimplexTree * st, int dim_skel)
      : sib_(nullptr),
        st_(st),
        dim_skel_(dim_skel),
        curr_dim_(0) {
    if (st == nullptr || st->root() == nullptr || st->root()->members().empty()) {
      st_ = nullptr;
    } else {
      sh_ = st->root()->members().begin();
      sib_ = st->root();
      while (st->has_children(sh_) && curr_dim_ < dim_skel_) {
        sib_ = sh_->second.children();
        sh_ = sib_->members().begin();
        ++curr_dim_;
      }
    }
  }
 private:
  friend class boost::iterator_core_access;

// valid when iterating along the SAME boundary.
  bool equal(Simplex_tree_skeleton_simplex_iterator const& other) const {
    if (other.st_ == nullptr) {
      return (st_ == nullptr);
    }
    if (st_ == nullptr) {
      return false;
    }
    return (&(sh_->second) == &(other.sh_->second));
  }

  Simplex_handle const& dereference() const {
    return sh_;
  }

// Depth first traversal of the skeleton.
  void increment() {
    ++sh_;
    if (sh_ == sib_->members().end()) {
      if (sib_->oncles() == nullptr) {
        st_ = nullptr;
        return;
      }  // reach the end
      sh_ = sib_->oncles()->members().find(sib_->parent());
      sib_ = sib_->oncles();
      --curr_dim_;
      return;
    }
    while (st_->has_children(sh_) && curr_dim_ < dim_skel_) {
      sib_ = sh_->second.children();
      sh_ = sib_->members().begin();
      ++curr_dim_;
    }
  }

  Simplex_handle sh_;
  Siblings * sib_;
  SimplexTree * st_;
  int dim_skel_;
  int curr_dim_;
};

/* @} */  // end addtogroup simplex_tree
}  // namespace Gudhi

#endif  // SIMPLEX_TREE_SIMPLEX_TREE_ITERATORS_H_
