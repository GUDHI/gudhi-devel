/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef SIMPLEX_TREE_H_
#define SIMPLEX_TREE_H_

#include <gudhi/Simplex_tree/Simplex_tree_node_explicit_storage.h>
#include <gudhi/Simplex_tree/Simplex_tree_siblings.h>
#include <gudhi/Simplex_tree/Simplex_tree_iterators.h>
#include <gudhi/Simplex_tree/indexing_tag.h>

#include <gudhi/reader_utils.h>
#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/Debug_utils.h>

#include <boost/container/flat_map.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/range/adaptor/reversed.hpp>

#ifdef GUDHI_USE_TBB
#include <tbb/parallel_sort.h>
#endif

#include <utility>
#include <vector>
#include <functional>  // for greater<>
#include <stdexcept>
#include <limits>  // Inf
#include <initializer_list>
#include <algorithm>  // for std::max
#include <cstdint>  // for std::uint32_t
#include <iterator>  // for std::distance

namespace Gudhi {

struct Simplex_tree_options_full_featured;

/**
 * \class Simplex_tree Simplex_tree.h gudhi/Simplex_tree.h
 * \brief Simplex Tree data structure for representing simplicial complexes.
 *
 * \details Every simplex \f$[v_0, \cdots ,v_d]\f$ admits a canonical orientation
 * induced by the order relation on vertices \f$ v_0 < \cdots < v_d \f$.
 *
 * Details may be found in \cite boissonnatmariasimplextreealgorithmica.
 *
 * \implements FilteredComplex
 *
 */

template<typename SimplexTreeOptions = Simplex_tree_options_full_featured>
class Simplex_tree {
 public:
  typedef SimplexTreeOptions Options;
  typedef typename Options::Indexing_tag Indexing_tag;
  /** \brief Type for the value of the filtration function.
   *
   * Must be comparable with <. */
  typedef typename Options::Filtration_value Filtration_value;
  /** \brief Key associated to each simplex.
   *
   * Must be an integer type. */
  typedef typename Options::Simplex_key Simplex_key;
  /** \brief Type for the vertex handle.
   *
   * Must be a signed integer type. It admits a total order <. */
  typedef typename Options::Vertex_handle Vertex_handle;

  /* Type of node in the simplex tree. */
  typedef Simplex_tree_node_explicit_storage<Simplex_tree> Node;
  /* Type of dictionary Vertex_handle -> Node for traversing the simplex tree. */
  // Note: this wastes space when Vertex_handle is 32 bits and Node is aligned on 64 bits. It would be better to use a
  // flat_set (with our own comparator) where we can control the layout of the struct (put Vertex_handle and
  // Simplex_key next to each other).
  typedef typename boost::container::flat_map<Vertex_handle, Node> Dictionary;

  /* \brief Set of nodes sharing a same parent in the simplex tree. */
  /* \brief Set of nodes sharing a same parent in the simplex tree. */
  typedef Simplex_tree_siblings<Simplex_tree, Dictionary> Siblings;

  struct Key_simplex_base_real {
    Key_simplex_base_real() : key_(-1) {}
    void assign_key(Simplex_key k) { key_ = k; }
    Simplex_key key() const { return key_; }
   private:
    Simplex_key key_;
  };
  struct Key_simplex_base_dummy {
    Key_simplex_base_dummy() {}
    // Undefined so it will not link
    void assign_key(Simplex_key);
    Simplex_key key() const;
  };
  typedef typename std::conditional<Options::store_key, Key_simplex_base_real, Key_simplex_base_dummy>::type
      Key_simplex_base;

  struct Filtration_simplex_base_real {
    Filtration_simplex_base_real() : filt_(0) {}
    void assign_filtration(Filtration_value f) { filt_ = f; }
    Filtration_value filtration() const { return filt_; }
   private:
    Filtration_value filt_;
  };
  struct Filtration_simplex_base_dummy {
    Filtration_simplex_base_dummy() {}
    void assign_filtration(Filtration_value GUDHI_CHECK_code(f)) { GUDHI_CHECK(f == 0, "filtration value specified for a complex that does not store them"); }
    Filtration_value filtration() const { return 0; }
  };
  typedef typename std::conditional<Options::store_filtration, Filtration_simplex_base_real,
    Filtration_simplex_base_dummy>::type Filtration_simplex_base;

 public:
  /** \brief Handle type to a simplex contained in the simplicial complex represented
   * by the simplex tree. */
  typedef typename Dictionary::iterator Simplex_handle;

 private:
  typedef typename Dictionary::iterator Dictionary_it;
  typedef typename Dictionary_it::value_type Dit_value_t;

  struct return_first {
    Vertex_handle operator()(const Dit_value_t& p_sh) const {
      return p_sh.first;
    }
  };

 public:
  /** \name Range and iterator types
   *
   * The naming convention is Container_content_(iterator/range). A Container_content_range is
   * essentially an object on which the methods begin() and end() can be called. They both return
   * an object of type Container_content_iterator, and allow the traversal of the range
   * [ begin();end() ).
   * @{ */

  /** \brief Iterator over the vertices of the simplicial complex.
   *
   * 'value_type' is Vertex_handle. */
  typedef boost::transform_iterator<return_first, Dictionary_it> Complex_vertex_iterator;
  /** \brief Range over the vertices of the simplicial complex. */
  typedef boost::iterator_range<Complex_vertex_iterator> Complex_vertex_range;
  /** \brief Iterator over the vertices of a simplex.
   *
   * 'value_type' is Vertex_handle. */
  typedef Simplex_tree_simplex_vertex_iterator<Simplex_tree> Simplex_vertex_iterator;
  /** \brief Range over the vertices of a simplex. */
  typedef boost::iterator_range<Simplex_vertex_iterator> Simplex_vertex_range;
  /** \brief Range over the cofaces of a simplex. */
  typedef std::vector<Simplex_handle> Cofaces_simplex_range;
  /** \brief Iterator over the simplices of the boundary of a simplex.
   *
   * 'value_type' is Simplex_handle. */
  typedef Simplex_tree_boundary_simplex_iterator<Simplex_tree> Boundary_simplex_iterator;
  /** \brief Range over the simplices of the boundary of a simplex. */
  typedef boost::iterator_range<Boundary_simplex_iterator> Boundary_simplex_range;
  /** \brief Iterator over the simplices of the simplicial complex.
   *
   * 'value_type' is Simplex_handle. */
  typedef Simplex_tree_complex_simplex_iterator<Simplex_tree> Complex_simplex_iterator;
  /** \brief Range over the simplices of the simplicial complex. */
  typedef boost::iterator_range<Complex_simplex_iterator> Complex_simplex_range;
  /** \brief Iterator over the simplices of the skeleton of the simplicial complex, for a given 
   * dimension.
   *
   * 'value_type' is Simplex_handle. */
  typedef Simplex_tree_skeleton_simplex_iterator<Simplex_tree> Skeleton_simplex_iterator;
  /** \brief Range over the simplices of the skeleton of the simplicial complex, for a given 
   * dimension. */
  typedef boost::iterator_range<Skeleton_simplex_iterator> Skeleton_simplex_range;
  /** \brief Range over the simplices of the simplicial complex, ordered by the filtration. */
  typedef std::vector<Simplex_handle> Filtration_simplex_range;
  /** \brief Iterator over the simplices of the simplicial complex, ordered by the filtration.
   *
   * 'value_type' is Simplex_handle. */
  typedef typename Filtration_simplex_range::const_iterator Filtration_simplex_iterator;

  /* @} */  // end name range and iterator types
  /** \name Range and iterator methods
   * @{ */

  /** \brief Returns a range over the vertices of the simplicial complex. 
   * The order is increasing according to < on Vertex_handles.*/
  Complex_vertex_range complex_vertex_range() {
    return Complex_vertex_range(
                                boost::make_transform_iterator(root_.members_.begin(), return_first()),
                                boost::make_transform_iterator(root_.members_.end(), return_first()));
  }

  /** \brief Returns a range over the simplices of the simplicial complex.
   *
   * In the Simplex_tree, the tree is traverse in a depth-first fashion.
   * Consequently, simplices are ordered according to lexicographic order on the list of
   * Vertex_handles of a simplex, read in increasing < order for Vertex_handles. */
  Complex_simplex_range complex_simplex_range() {
    return Complex_simplex_range(Complex_simplex_iterator(this),
                                 Complex_simplex_iterator());
  }

  /** \brief Returns a range over the simplices of the dim-skeleton of the simplicial complex.
   *
   * The \f$d\f$-skeleton of a simplicial complex \f$\mathbf{K}\f$ is the simplicial complex containing the
   * simplices of \f$\mathbf{K}\f$ of dimension at most \f$d\f$.
   *
   * @param[in] dim The maximal dimension of the simplices in the skeleton.
   *
   * The simplices are ordered according to lexicographic order on the list of
   * Vertex_handles of a simplex, read in increasing < order for Vertex_handles. */
  Skeleton_simplex_range skeleton_simplex_range(int dim) {
    return Skeleton_simplex_range(Skeleton_simplex_iterator(this, dim),
                                  Skeleton_simplex_iterator());
  }

  /** \brief Returns a range over the simplices of the simplicial complex,
   * in the order of the filtration.
   *
   * The filtration is a monotonic function \f$ f: \mathbf{K} \rightarrow \mathbb{R} \f$, i.e. if two simplices
   * \f$\tau\f$ and \f$\sigma\f$ satisfy \f$\tau \subseteq \sigma\f$ then
   * \f$f(\tau) \leq f(\sigma)\f$.
   *
   * The method returns simplices ordered according to increasing filtration values. Ties are
   * resolved by considering inclusion relation (subsimplices appear before their cofaces). If two
   * simplices have same filtration value but are not comparable w.r.t. inclusion, lexicographic
   * order is used.
   *
   * The filtration must be valid. If the filtration has not been initialized yet, the
   * method initializes it (i.e. order the simplices). If the complex has changed since the last time the filtration
   * was initialized, please call `initialize_filtration()` to recompute it. */
  Filtration_simplex_range const& filtration_simplex_range(Indexing_tag = Indexing_tag()) {
    if (filtration_vect_.empty()) {
      initialize_filtration();
    }
    return filtration_vect_;
  }

  /** \brief Returns a range over the vertices of a simplex.
   *
   * The order in which the vertices are visited is the decreasing order for < on Vertex_handles,
   * which is consequenlty
   * equal to \f$(-1)^{\text{dim} \sigma}\f$ the canonical orientation on the simplex.
   */
  Simplex_vertex_range simplex_vertex_range(Simplex_handle sh) {
    assert(sh != null_simplex());  // Empty simplex
    return Simplex_vertex_range(Simplex_vertex_iterator(this, sh),
                                Simplex_vertex_iterator(this));
  }

  /** \brief Returns a range over the simplices of the boundary of a simplex.
   *
   * The boundary of a simplex is the set of codimension \f$1\f$ subsimplices of the simplex.
   * If the simplex is \f$[v_0, \cdots ,v_d]\f$, with canonical orientation
   * induced by \f$ v_0 < \cdots < v_d \f$, the iterator enumerates the
   * simplices of the boundary in the order:
   * \f$[v_0,\cdots,\widehat{v_i},\cdots,v_d]\f$ for \f$i\f$ from \f$0\f$ to \f$d\f$,
   * where \f$\widehat{v_i}\f$ means that the vertex \f$v_i\f$ is omitted.
   *
   * We note that the alternate sum of the simplices given by the iterator
   * gives \f$(-1)^{\text{dim} \sigma}\f$ the chains corresponding to the boundary
   * of the simplex.
   *
   * @param[in] sh Simplex for which the boundary is computed. */
  template<class SimplexHandle>
  Boundary_simplex_range boundary_simplex_range(SimplexHandle sh) {
    return Boundary_simplex_range(Boundary_simplex_iterator(this, sh),
                                  Boundary_simplex_iterator(this));
  }

  /** @} */  // end range and iterator methods
  /** \name Constructor/Destructor
   * @{ */

  /** \brief Constructs an empty simplex tree. */
  Simplex_tree()
      : null_vertex_(-1),
      root_(nullptr, null_vertex_),
      filtration_vect_(),
      dimension_(-1) { }

  /** \brief User-defined copy constructor reproduces the whole tree structure. */
  Simplex_tree(const Simplex_tree& complex_source) {
#ifdef DEBUG_TRACES
    std::cout << "Simplex_tree copy constructor" << std::endl;
#endif  // DEBUG_TRACES
    copy_from(complex_source);
  }

  /** \brief User-defined move constructor relocates the whole tree structure.
   *  \exception std::invalid_argument In debug mode, if the complex_source is invalid.
   */
  Simplex_tree(Simplex_tree && complex_source) {
#ifdef DEBUG_TRACES
    std::cout << "Simplex_tree move constructor" << std::endl;
#endif  // DEBUG_TRACES
    move_from(complex_source);

    // just need to set dimension_ on source to make it available again
    // (filtration_vect_ and members are already set from the move)
    complex_source.dimension_ = -1;
  }

  /** \brief Destructor; deallocates the whole tree structure. */
  ~Simplex_tree() {
    root_members_recursive_deletion();
  }

  /** \brief User-defined copy assignment reproduces the whole tree structure. */
  Simplex_tree& operator= (const Simplex_tree& complex_source) {
#ifdef DEBUG_TRACES
    std::cout << "Simplex_tree copy assignment" << std::endl;
#endif  // DEBUG_TRACES
    // Self-assignment detection
    if (&complex_source != this) {
      // We start by deleting root_ if not empty
      root_members_recursive_deletion();

      copy_from(complex_source);
    }
    return *this;
  }

  /** \brief User-defined move assignment relocates the whole tree structure.
   *  \exception std::invalid_argument In debug mode, if the complex_source is invalid.
   */
  Simplex_tree& operator=(Simplex_tree&& complex_source) {
#ifdef DEBUG_TRACES
    std::cout << "Simplex_tree move assignment" << std::endl;
#endif  // DEBUG_TRACES
    // Self-assignment detection
    if (&complex_source != this) {
      // root_ deletion in case it was not empty
      root_members_recursive_deletion();

      move_from(complex_source);
    }
    return *this;
  }
  /** @} */  // end constructor/destructor

 private:
  // Copy from complex_source to "this"
  void copy_from(const Simplex_tree& complex_source) {
    null_vertex_ = complex_source.null_vertex_;
    filtration_vect_.clear();
    dimension_ = complex_source.dimension_;
    auto root_source = complex_source.root_;

    // root members copy
    root_.members() = Dictionary(boost::container::ordered_unique_range, root_source.members().begin(), root_source.members().end());
    // Needs to reassign children
    for (auto& map_el : root_.members()) {
      map_el.second.assign_children(&root_);
    }
    rec_copy(&root_, &root_source);
  }

  /** \brief depth first search, inserts simplices when reaching a leaf. */
  void rec_copy(Siblings *sib, Siblings *sib_source) {
    for (auto sh = sib->members().begin(), sh_source = sib_source->members().begin();
         sh != sib->members().end(); ++sh, ++sh_source) {
      if (has_children(sh_source)) {
        Siblings * newsib = new Siblings(sib, sh_source->first);
        newsib->members_.reserve(sh_source->second.children()->members().size());
        for (auto & child : sh_source->second.children()->members())
          newsib->members_.emplace_hint(newsib->members_.end(), child.first, Node(newsib, child.second.filtration()));
        rec_copy(newsib, sh_source->second.children());
        sh->second.assign_children(newsib);
      }
    }
  }

  // Move from complex_source to "this"
  void move_from(Simplex_tree& complex_source) {
    null_vertex_ = std::move(complex_source.null_vertex_);
    root_ = std::move(complex_source.root_);
    filtration_vect_ = std::move(complex_source.filtration_vect_);
    dimension_ = std::move(complex_source.dimension_);

    // Need to update root members (children->oncles and children need to point on the new root pointer)
    for (auto& map_el : root_.members()) {
      if (map_el.second.children() != &(complex_source.root_)) {
        // reset children->oncles with the moved root_ pointer value
        map_el.second.children()->oncles_ = &root_;
      } else {
        // if simplex is of dimension 0, oncles_ shall be nullptr
        GUDHI_CHECK(map_el.second.children()->oncles_ == nullptr,
                    std::invalid_argument("Simplex_tree move constructor from an invalid Simplex_tree"));
        // and children points on root_ - to be moved
        map_el.second.assign_children(&root_);
      }
    }
  }

  // delete all root_.members() recursively
  void root_members_recursive_deletion() {
    for (auto sh = root_.members().begin(); sh != root_.members().end(); ++sh) {
      if (has_children(sh)) {
        rec_delete(sh->second.children());
      }
    }
    root_.members().clear();
  }

  // Recursive deletion
  void rec_delete(Siblings * sib) {
    for (auto sh = sib->members().begin(); sh != sib->members().end(); ++sh) {
      if (has_children(sh)) {
        rec_delete(sh->second.children());
      }
    }
    delete sib;
  }

 public:
  /** \brief Checks if two simplex trees are equal. */
  bool operator==(Simplex_tree& st2) {
    if ((null_vertex_ != st2.null_vertex_) ||
        (dimension_ != st2.dimension_))
      return false;
    return rec_equal(&root_, &st2.root_);
  }

  /** \brief Checks if two simplex trees are different. */
  bool operator!=(Simplex_tree& st2) {
    return (!(*this == st2));
  }

 private:
  /** rec_equal: Checks recursively whether or not two simplex trees are equal, using depth first search. */
  bool rec_equal(Siblings* s1, Siblings* s2) {
    if (s1->members().size() != s2->members().size())
      return false;
    for (auto sh1 = s1->members().begin(), sh2 = s2->members().begin();
         (sh1 != s1->members().end() && sh2 != s2->members().end()); ++sh1, ++sh2) {
      if (sh1->first != sh2->first || sh1->second.filtration() != sh2->second.filtration())
        return false;
      if (has_children(sh1) != has_children(sh2))
        return false;
      // Recursivity on children only if both have children
      else if (has_children(sh1))
        if (!rec_equal(sh1->second.children(), sh2->second.children()))
          return false;
    }
    return true;
  }

 public:
  /** \brief Returns the key associated to a simplex.
   *
   * The filtration must be initialized.
   * \pre SimplexTreeOptions::store_key
   */
  static Simplex_key key(Simplex_handle sh) {
    return sh->second.key();
  }

  /** \brief Returns the simplex that has index idx in the filtration.
   *
   * The filtration must be initialized.
   * \pre SimplexTreeOptions::store_key
   */
  Simplex_handle simplex(Simplex_key idx) const {
    return filtration_vect_[idx];
  }

  /** \brief Returns the filtration value of a simplex.
   *
   * Called on the null_simplex, it returns infinity.
   * If SimplexTreeOptions::store_filtration is false, returns 0.
   */
  static Filtration_value filtration(Simplex_handle sh) {
    if (sh != null_simplex()) {
      return sh->second.filtration();
    } else {
      return std::numeric_limits<Filtration_value>::infinity();
    }
  }

  /** \brief Sets the filtration value of a simplex.
   * \exception std::invalid_argument In debug mode, if sh is a null_simplex.
   */
  void assign_filtration(Simplex_handle sh, Filtration_value fv) {
    GUDHI_CHECK(sh != null_simplex(),
                std::invalid_argument("Simplex_tree::assign_filtration - cannot assign filtration on null_simplex"));
    sh->second.assign_filtration(fv);
  }

  /** \brief Returns a Simplex_handle different from all Simplex_handles
   * associated to the simplices in the simplicial complex.
   *
   * One can call filtration(null_simplex()). */
  static Simplex_handle null_simplex() {
    return Dictionary_it(nullptr);
  }

  /** \brief Returns a key different for all keys associated to the
   * simplices of the simplicial complex. */
  static Simplex_key null_key() {
    return -1;
  }

  /** \brief Returns a Vertex_handle different from all Vertex_handles associated
   * to the vertices of the simplicial complex. */
  Vertex_handle null_vertex() const {
    return null_vertex_;
  }

  /** \brief Returns the number of vertices in the complex. */
  size_t num_vertices() const {
    return root_.members_.size();
  }

 public:
  /** \brief returns the number of simplices in the simplex_tree. */
  size_t num_simplices() {
    return num_simplices(&root_);
  }

 private:
  /** \brief returns the number of simplices in the simplex_tree. */
  size_t num_simplices(Siblings * sib) {
    auto sib_begin = sib->members().begin();
    auto sib_end = sib->members().end();
    size_t simplices_number = sib_end - sib_begin;
    for (auto sh = sib_begin; sh != sib_end; ++sh) {
      if (has_children(sh)) {
        simplices_number += num_simplices(sh->second.children());
      }
    }
    return simplices_number;
  }

 public:
  /** \brief Returns the dimension of a simplex.
   *
   * Must be different from null_simplex().*/
  int dimension(Simplex_handle sh) {
    Siblings * curr_sib = self_siblings(sh);
    int dim = 0;
    while (curr_sib != nullptr) {
      ++dim;
      curr_sib = curr_sib->oncles();
    }
    return dim - 1;
  }

  /** \brief Returns an upper bound on the dimension of the simplicial complex. */
  int upper_bound_dimension() const {
    return dimension_;
  }

  /** \brief Returns the dimension of the simplicial complex.
      \details This function is not constant time because it can recompute dimension if required (can be triggered by
      `remove_maximal_simplex()` or `prune_above_filtration()`).
  */
  int dimension() {
    if (dimension_to_be_lowered_)
      lower_upper_bound_dimension();
    return dimension_;
  }

  /** \brief Returns true if the node in the simplex tree pointed by
   * sh has children.*/
  template<class SimplexHandle>
  bool has_children(SimplexHandle sh) const {
    // Here we rely on the root using null_vertex(), which cannot match any real vertex.
    return (sh->second.children()->parent() == sh->first);
  }

    /** \brief Given a range of Vertex_handles, returns the Simplex_handle
   * of the simplex in the simplicial complex containing the corresponding
   * vertices. Return null_simplex() if the simplex is not in the complex.
   *
   * The type InputVertexRange must be a range of <CODE>Vertex_handle</CODE>
   * on which we can call std::begin() function
   */
  template<class InputVertexRange = std::initializer_list<Vertex_handle>>
  Simplex_handle find(const InputVertexRange & s) {
    auto first = std::begin(s);
    auto last = std::end(s);

    if (first == last)
      return null_simplex();  // ----->>

    // Copy before sorting
    std::vector<Vertex_handle> copy(first, last);
    std::sort(std::begin(copy), std::end(copy));
    return find_simplex(copy);
  }

 private:
  /** Find function, with a sorted range of vertices. */
  Simplex_handle find_simplex(const std::vector<Vertex_handle> & simplex) {
    Siblings * tmp_sib = &root_;
    Dictionary_it tmp_dit;
    auto vi = simplex.begin();
    if (Options::contiguous_vertices) {
      // Equivalent to the first iteration of the normal loop
      GUDHI_CHECK(contiguous_vertices(), "non-contiguous vertices");
      Vertex_handle v = *vi++;
      if(v < 0 || v >= static_cast<Vertex_handle>(root_.members_.size()))
        return null_simplex();
      tmp_dit = root_.members_.begin() + v;
      if (vi == simplex.end())
        return tmp_dit;
      if (!has_children(tmp_dit))
        return null_simplex();
      tmp_sib = tmp_dit->second.children();
    }
    for (;;) {
      tmp_dit = tmp_sib->members_.find(*vi++);
      if (tmp_dit == tmp_sib->members_.end())
        return null_simplex();
      if (vi == simplex.end())
        return tmp_dit;
      if (!has_children(tmp_dit))
        return null_simplex();
      tmp_sib = tmp_dit->second.children();
    }
  }

  /** \brief Returns the Simplex_handle corresponding to the 0-simplex
   * representing the vertex with Vertex_handle v. */
  Simplex_handle find_vertex(Vertex_handle v) {
    if (Options::contiguous_vertices) {
      assert(contiguous_vertices());
      return root_.members_.begin() + v;
    } else {
      return root_.members_.find(v);
    }
  }

 public:
  /** \private \brief Test if the vertices have contiguous numbering: 0, 1, etc.  */
  bool contiguous_vertices() const {
    if (root_.members_.empty()) return true;
    if (root_.members_.begin()->first != 0) return false;
    if (std::prev(root_.members_.end())->first != static_cast<Vertex_handle>(root_.members_.size() - 1)) return false;
    return true;
  }

 private:
  /** \brief Inserts a simplex represented by a vector of vertex.
   * @param[in]  simplex    vector of Vertex_handles, representing the vertices of the new simplex. The vector must be
   * sorted by increasing vertex handle order.
   * @param[in]  filtration the filtration value assigned to the new simplex.
   * @return If the new simplex is inserted successfully (i.e. it was not in the
   * simplicial complex yet) the bool is set to true and the Simplex_handle is the handle assigned
   * to the new simplex.
   * If the insertion fails (the simplex is already there), the bool is set to false. If the insertion
   * fails and the simplex already in the complex has a filtration value strictly bigger than 'filtration',
   * we assign this simplex with the new value 'filtration', and set the Simplex_handle field of the
   * output pair to the Simplex_handle of the simplex. Otherwise, we set the Simplex_handle part to
   * null_simplex.
   * 
  */
  std::pair<Simplex_handle, bool> insert_vertex_vector(const std::vector<Vertex_handle>& simplex,
                                                     Filtration_value filtration) {
    Siblings * curr_sib = &root_;
    std::pair<Simplex_handle, bool> res_insert;
    auto vi = simplex.begin();
    for (; vi != simplex.end() - 1; ++vi) {
      GUDHI_CHECK(*vi != null_vertex(), "cannot use the dummy null_vertex() as a real vertex");
      res_insert = curr_sib->members_.emplace(*vi, Node(curr_sib, filtration));
      if (!(has_children(res_insert.first))) {
        res_insert.first->second.assign_children(new Siblings(curr_sib, *vi));
      }
      curr_sib = res_insert.first->second.children();
    }
    GUDHI_CHECK(*vi != null_vertex(), "cannot use the dummy null_vertex() as a real vertex");
    res_insert = curr_sib->members_.emplace(*vi, Node(curr_sib, filtration));
    if (!res_insert.second) {
      // if already in the complex
      if (res_insert.first->second.filtration() > filtration) {
        // if filtration value modified
        res_insert.first->second.assign_filtration(filtration);
        return res_insert;
      }
      // if filtration value unchanged
      return std::pair<Simplex_handle, bool>(null_simplex(), false);
    }
    // otherwise the insertion has succeeded - size is a size_type
    if (static_cast<int>(simplex.size()) - 1 > dimension_) {
      // Update dimension if needed
      dimension_ = static_cast<int>(simplex.size()) - 1;
    }
    return res_insert;
  }

 public:
  /** \brief Insert a simplex, represented by a range of Vertex_handles, in the simplicial complex.
   *
   * @param[in]  simplex    range of Vertex_handles, representing the vertices of the new simplex
   * @param[in]  filtration the filtration value assigned to the new simplex.
   * @return If the new simplex is inserted successfully (i.e. it was not in the
   * simplicial complex yet) the bool is set to true and the Simplex_handle is the handle assigned
   * to the new simplex.
   * If the insertion fails (the simplex is already there), the bool is set to false. If the insertion
   * fails and the simplex already in the complex has a filtration value strictly bigger than 'filtration',
   * we assign this simplex with the new value 'filtration', and set the Simplex_handle field of the
   * output pair to the Simplex_handle of the simplex. Otherwise, we set the Simplex_handle part to
   * null_simplex.
   *
   * All subsimplices do not necessary need to be already in the simplex tree to proceed to an
   * insertion. However, the property of being a simplicial complex will be violated. This allows
   * us to insert a stream of simplices contained in a simplicial complex without considering any
   * order on them.
   *
   * The filtration value
   * assigned to the new simplex must preserve the monotonicity of the filtration.
   *
   * The type InputVertexRange must be a range for which .begin() and
   * .end() return input iterators, with 'value_type' Vertex_handle. */
  template<class InputVertexRange = std::initializer_list<Vertex_handle>>
  std::pair<Simplex_handle, bool> insert_simplex(const InputVertexRange & simplex,
                                                 Filtration_value filtration = 0) {
    auto first = std::begin(simplex);
    auto last = std::end(simplex);

    if (first == last)
      return std::pair<Simplex_handle, bool>(null_simplex(), true);  // ----->>

    // Copy before sorting
    std::vector<Vertex_handle> copy(first, last);
    std::sort(std::begin(copy), std::end(copy));
    return insert_vertex_vector(copy, filtration);
  }

    /** \brief Insert a N-simplex and all his subfaces, from a N-simplex represented by a range of
   * Vertex_handles, in the simplicial complex.
   *
   * @param[in]  Nsimplex   range of Vertex_handles, representing the vertices of the new N-simplex
   * @param[in]  filtration the filtration value assigned to the new N-simplex.
   * @return If the new simplex is inserted successfully (i.e. it was not in the
   * simplicial complex yet) the bool is set to true and the Simplex_handle is the handle assigned
   * to the new simplex.
   * If the insertion fails (the simplex is already there), the bool is set to false. If the insertion
   * fails and the simplex already in the complex has a filtration value strictly bigger than 'filtration',
   * we assign this simplex with the new value 'filtration', and set the Simplex_handle field of the
   * output pair to the Simplex_handle of the simplex. Otherwise, we set the Simplex_handle part to
   * null_simplex.
   */
  template<class InputVertexRange = std::initializer_list<Vertex_handle>>
  std::pair<Simplex_handle, bool> insert_simplex_and_subfaces(const InputVertexRange& Nsimplex,
				    Filtration_value filtration = 0) {
    auto first = std::begin(Nsimplex);
    auto last = std::end(Nsimplex);

    if (first == last)
      return { null_simplex(), true }; // FIXME: false would make more sense to me.

    // Copy before sorting
    // Thread local is not available on XCode version < V.8 - It will slow down computation
#ifdef GUDHI_CAN_USE_CXX11_THREAD_LOCAL
    thread_local
#endif  // GUDHI_CAN_USE_CXX11_THREAD_LOCAL
    std::vector<Vertex_handle> copy;
    copy.clear();
    copy.insert(copy.end(), first, last);
    std::sort(copy.begin(), copy.end());
    auto last_unique = std::unique(copy.begin(), copy.end());
    copy.erase(last_unique, copy.end());
    GUDHI_CHECK_code(
      for (Vertex_handle v : copy)
        GUDHI_CHECK(v != null_vertex(), "cannot use the dummy null_vertex() as a real vertex");
    )
    // Update dimension if needed. We could wait to see if the insertion succeeds, but I doubt there is much to gain.
    dimension_ = (std::max)(dimension_, static_cast<int>(std::distance(copy.begin(), copy.end())) - 1);

    return rec_insert_simplex_and_subfaces_sorted(root(), copy.begin(), copy.end(), filtration);
  }

 private:
  // To insert {1,2,3,4}, we insert {2,3,4} twice, once at the root, and once below 1.
  template<class ForwardVertexIterator>
  std::pair<Simplex_handle, bool> rec_insert_simplex_and_subfaces_sorted(Siblings* sib,
  	                                                                     ForwardVertexIterator first,
  	                                                                     ForwardVertexIterator last,
  	                                                                     Filtration_value filt) {
    // An alternative strategy would be:
    // - try to find the complete simplex, if found (and low filtration) exit
    // - insert all the vertices at once in sib
    // - loop over those (new or not) simplices, with a recursive call(++first, last)
    Vertex_handle vertex_one = *first;
    auto&& dict = sib->members();
    auto insertion_result = dict.emplace(vertex_one, Node(sib, filt));
    Simplex_handle simplex_one = insertion_result.first;
    bool one_is_new = insertion_result.second;
    if (!one_is_new) {
      if (filtration(simplex_one) > filt) {
        assign_filtration(simplex_one, filt);
      } else {
        // FIXME: this interface makes no sense, and it doesn't seem to be tested.
        insertion_result.first = null_simplex();
      }
    }
    if (++first == last) return insertion_result;
    if (!has_children(simplex_one))
      // TODO: have special code here, we know we are building the whole subtree from scratch.
      simplex_one->second.assign_children(new Siblings(sib, vertex_one));
    auto res = rec_insert_simplex_and_subfaces_sorted(simplex_one->second.children(), first, last, filt);
    // No need to continue if the full simplex was already there with a low enough filtration value.
    if (res.first != null_simplex()) rec_insert_simplex_and_subfaces_sorted(sib, first, last, filt);
    return res;
  }

 public:
  /** \brief Assign a value 'key' to the key of the simplex
   * represented by the Simplex_handle 'sh'. */
  void assign_key(Simplex_handle sh, Simplex_key key) {
    sh->second.assign_key(key);
  }

  /** Returns the two Simplex_handle corresponding to the endpoints of
   * and edge. sh must point to a 1-dimensional simplex. This is an
   * optimized version of the boundary computation. */
  std::pair<Simplex_handle, Simplex_handle> endpoints(Simplex_handle sh) {
    assert(dimension(sh) == 1);
    return { find_vertex(sh->first), find_vertex(self_siblings(sh)->parent()) };
  }

  /** Returns the Siblings containing a simplex.*/
  template<class SimplexHandle>
  Siblings* self_siblings(SimplexHandle sh) {
    if (sh->second.children()->parent() == sh->first)
      return sh->second.children()->oncles();
    else
      return sh->second.children();
  }

 public:
  /** Returns a pointer to the root nodes of the simplex tree. */
  Siblings * root() {
    return &root_;
  }

  /** \brief Set a dimension for the simplicial complex.
   *  \details This function must be used with caution because it disables dimension recomputation when required
   * (this recomputation can be triggered by `remove_maximal_simplex()` or `prune_above_filtration()`).
   */
  void set_dimension(int dimension) {
    dimension_to_be_lowered_ = false;
    dimension_ = dimension;
  }

 public:
  /** \brief Initializes the filtrations, i.e. sort the
   * simplices according to their order in the filtration and initializes all Simplex_keys.
   *
   * After calling this method, filtration_simplex_range() becomes valid, and each simplex is
   * assigned a Simplex_key corresponding to its order in the filtration (from 0 to m-1 for a
   * simplicial complex with m simplices).
   *
   * Will be automatically called when calling filtration_simplex_range()
   * if the filtration has never been initialized yet. */
  void initialize_filtration() {
    filtration_vect_.clear();
    filtration_vect_.reserve(num_simplices());
    for (Simplex_handle sh : complex_simplex_range())
      filtration_vect_.push_back(sh);

    /* We use stable_sort here because with libstdc++ it is faster than sort.
     * is_before_in_filtration is now a total order, but we used to call
     * stable_sort for the following heuristic:
     * The use of a depth-first traversal of the simplex tree, provided by
     * complex_simplex_range(), combined with a stable sort is meant to
     * optimize the order of simplices with same filtration value. The
     * heuristic consists in inserting the cofaces of a simplex as soon as
     * possible.
     */
#ifdef GUDHI_USE_TBB
    tbb::parallel_sort(filtration_vect_.begin(), filtration_vect_.end(), is_before_in_filtration(this));
#else
    std::stable_sort(filtration_vect_.begin(), filtration_vect_.end(), is_before_in_filtration(this));
#endif
  }

 private:
  /** Recursive search of cofaces
   * This function uses DFS
   *\param vertices contains a list of vertices, which represent the vertices of the simplex not found yet.
   *\param curr_nbVertices represents the number of vertices of the simplex we reached by going through the tree.
   *\param cofaces contains a list of Simplex_handle, representing all the cofaces asked.
   *\param star true if we need the star of the simplex
   *\param nbVertices number of vertices of the cofaces we search
   * Prefix actions : When the bottom vertex matches with the current vertex in the tree, we remove the bottom vertex from vertices.
   * Infix actions : Then we call or not the recursion.
   * Postfix actions : Finally, we add back the removed vertex into vertices, and remove this vertex from curr_nbVertices so that we didn't change the parameters.
   * If the vertices list is empty, we need to check if curr_nbVertices matches with the dimension of the cofaces asked.
   */
  void rec_coface(std::vector<Vertex_handle> &vertices, Siblings *curr_sib, int curr_nbVertices,
                  std::vector<Simplex_handle>& cofaces, bool star, int nbVertices) {
    if (!(star || curr_nbVertices <= nbVertices))  // dimension of actual simplex <= nbVertices
      return;
    for (Simplex_handle simplex = curr_sib->members().begin(); simplex != curr_sib->members().end(); ++simplex) {
      if (vertices.empty()) {
        // If we reached the end of the vertices, and the simplex has more vertices than the given simplex
        // => we found a coface

        // Add a coface if we wan't the star or if the number of vertices of the current simplex matches with nbVertices
        bool addCoface = (star || curr_nbVertices == nbVertices);
        if (addCoface)
          cofaces.push_back(simplex);
        if ((!addCoface || star) && has_children(simplex))  // Rec call
          rec_coface(vertices, simplex->second.children(), curr_nbVertices + 1, cofaces, star, nbVertices);
      } else {
        if (simplex->first == vertices.back()) {
          // If curr_sib matches with the top vertex
          bool equalDim = (star || curr_nbVertices == nbVertices);  // dimension of actual simplex == nbVertices
          bool addCoface = vertices.size() == 1 && equalDim;
          if (addCoface)
            cofaces.push_back(simplex);
          if ((!addCoface || star) && has_children(simplex)) {
            // Rec call
            Vertex_handle tmp = vertices.back();
            vertices.pop_back();
            rec_coface(vertices, simplex->second.children(), curr_nbVertices + 1, cofaces, star, nbVertices);
            vertices.push_back(tmp);
          }
        } else if (simplex->first > vertices.back()) {
          return;
        } else {
          // (simplex->first < vertices.back()
          if (has_children(simplex))
            rec_coface(vertices, simplex->second.children(), curr_nbVertices + 1, cofaces, star, nbVertices);
        }
      }
    }
  }

 public:
  /** \brief Compute the star of a n simplex
   * \param simplex represent the simplex of which we search the star
   * \return Vector of Simplex_handle, empty vector if no cofaces found.
   */

  Cofaces_simplex_range star_simplex_range(const Simplex_handle simplex) {
    return cofaces_simplex_range(simplex, 0);
  }

  /** \brief Compute the cofaces of a n simplex
   * \param simplex represent the n-simplex of which we search the n+codimension cofaces
   * \param codimension The function returns the n+codimension-cofaces of the n-simplex. If codimension = 0, 
   * return all cofaces (equivalent of star function)
   * \return Vector of Simplex_handle, empty vector if no cofaces found.
   */

  Cofaces_simplex_range cofaces_simplex_range(const Simplex_handle simplex, int codimension) {
    Cofaces_simplex_range cofaces;
    // codimension must be positive or null integer
    assert(codimension >= 0);
    Simplex_vertex_range rg = simplex_vertex_range(simplex);
    std::vector<Vertex_handle> copy(rg.begin(), rg.end());
    if (codimension + static_cast<int>(copy.size()) > dimension_ + 1 ||
        (codimension == 0 && static_cast<int>(copy.size()) > dimension_))  // n+codimension greater than dimension_
      return cofaces;
    // must be sorted in decreasing order
    assert(std::is_sorted(copy.begin(), copy.end(), std::greater<Vertex_handle>()));
    bool star = codimension == 0;
    rec_coface(copy, &root_, 1, cofaces, star, codimension + static_cast<int>(copy.size()));
    return cofaces;
  }

 private:
  /** \brief Returns true iff the list of vertices of sh1
   * is smaller than the list of vertices of sh2 w.r.t.
   * lexicographic order on the lists read in reverse.
   *
   * It defines a StrictWeakOrdering on simplices. The Simplex_vertex_iterators
   * must traverse the Vertex_handle in decreasing order. Reverse lexicographic order satisfy
   * the property that a subsimplex of a simplex is always strictly smaller with this order. */
  bool reverse_lexicographic_order(Simplex_handle sh1, Simplex_handle sh2) {
    Simplex_vertex_range rg1 = simplex_vertex_range(sh1);
    Simplex_vertex_range rg2 = simplex_vertex_range(sh2);
    Simplex_vertex_iterator it1 = rg1.begin();
    Simplex_vertex_iterator it2 = rg2.begin();
    while (it1 != rg1.end() && it2 != rg2.end()) {
      if (*it1 == *it2) {
        ++it1;
        ++it2;
      } else {
        return *it1 < *it2;
      }
    }
    return ((it1 == rg1.end()) && (it2 != rg2.end()));
  }

  /** \brief StrictWeakOrdering, for the simplices, defined by the filtration.
   *
   * It corresponds to the partial order
   * induced by the filtration values, with ties resolved using reverse lexicographic order.
   * Reverse lexicographic order has the property to always consider the subsimplex of a simplex
   * to be smaller. The filtration function must be monotonic. */
  struct is_before_in_filtration {
    explicit is_before_in_filtration(Simplex_tree * st)
        : st_(st) { }

    bool operator()(const Simplex_handle sh1, const Simplex_handle sh2) const {
      // Not using st_->filtration(sh1) because it uselessly tests for null_simplex.
      if (sh1->second.filtration() != sh2->second.filtration()) {
        return sh1->second.filtration() < sh2->second.filtration();
      }
      // is sh1 a proper subface of sh2
      return st_->reverse_lexicographic_order(sh1, sh2);
    }

    Simplex_tree * st_;
  };

 public:
  /** \brief Inserts a 1-skeleton in an empty Simplex_tree.
   *
   * The Simplex_tree must contain no simplex when the method is
   * called.
   *
   * Inserts all vertices and edges given by a OneSkeletonGraph.
   * OneSkeletonGraph must be a model of
   * <a href="http://www.boost.org/doc/libs/1_65_1/libs/graph/doc/EdgeListGraph.html">boost::EdgeListGraph</a>
   * and <a href="http://www.boost.org/doc/libs/1_65_1/libs/graph/doc/PropertyGraph.html">boost::PropertyGraph</a>.
   *
   * The vertex filtration value is accessible through the property tag
   * vertex_filtration_t.
   * The edge filtration value is accessible through the property tag
   * edge_filtration_t.
   *
   * boost::graph_traits<OneSkeletonGraph>::vertex_descriptor
   *                                    must be Vertex_handle.
   * boost::graph_traits<OneSkeletonGraph>::directed_category
   *                                    can be directed_tag (the fastest, the least RAM use), undirected_tag or even
   *                                    bidirected_tag.
   *
   * If an edge appears with multiplicity, the function will arbitrarily pick
   * one representative to read the filtration value.  */
  template<class OneSkeletonGraph>
  void insert_graph(const OneSkeletonGraph& skel_graph) {
    // the simplex tree must be empty
    assert(num_simplices() == 0);

    if (boost::num_vertices(skel_graph) == 0) {
      return;
    }
    if (num_edges(skel_graph) == 0) {
      dimension_ = 0;
    } else {
      dimension_ = 1;
    }

    root_.members_.reserve(boost::num_vertices(skel_graph));

    typename boost::graph_traits<OneSkeletonGraph>::vertex_iterator v_it,
        v_it_end;
    for (std::tie(v_it, v_it_end) = boost::vertices(skel_graph); v_it != v_it_end;
         ++v_it) {
      root_.members_.emplace_hint(
                                  root_.members_.end(), *v_it,
                                  Node(&root_, boost::get(vertex_filtration_t(), skel_graph, *v_it)));
    }
    std::pair<typename boost::graph_traits<OneSkeletonGraph>::edge_iterator,
              typename boost::graph_traits<OneSkeletonGraph>::edge_iterator> boost_edges = boost::edges(skel_graph);
    // boost_edges.first is the equivalent to boost_edges.begin()
    // boost_edges.second is the equivalent to boost_edges.end()
    for (; boost_edges.first != boost_edges.second; boost_edges.first++) {
      auto edge = *(boost_edges.first);
      auto u = source(edge, skel_graph);
      auto v = target(edge, skel_graph);
      if (u == v) throw "Self-loops are not simplicial";
      // We cannot skip edges with the wrong orientation and expect them to
      // come a second time with the right orientation, that does not always
      // happen in practice. emplace() should be a NOP when an element with the
      // same key is already there, so seeing the same edge multiple times is
      // ok.
      // Should we actually forbid multiple edges? That would be consistent
      // with rejecting self-loops.
      if (v < u) std::swap(u, v);
      auto sh = find_vertex(u);
      if (!has_children(sh)) {
        sh->second.assign_children(new Siblings(&root_, sh->first));
      }

      sh->second.children()->members().emplace(v,
          Node(sh->second.children(), boost::get(edge_filtration_t(), skel_graph, edge)));
    }
  }

  /** \brief Expands the Simplex_tree containing only its one skeleton
   * until dimension max_dim.
   *
   * The expanded simplicial complex until dimension \f$d\f$
   * attached to a graph \f$G\f$ is the maximal simplicial complex of
   * dimension at most \f$d\f$ admitting the graph \f$G\f$ as \f$1\f$-skeleton.
   * The filtration value assigned to a simplex is the maximal filtration
   * value of one of its edges.
   *
   * The Simplex_tree must contain no simplex of dimension bigger than
   * 1 when calling the method. */
  void expansion(int max_dim) {
    if (max_dim <= 1) return;
    dimension_ = max_dim;
    for (Dictionary_it root_it = root_.members_.begin();
         root_it != root_.members_.end(); ++root_it) {
      if (has_children(root_it)) {
        siblings_expansion(root_it->second.children(), max_dim - 1);
      }
    }
    dimension_ = max_dim - dimension_;
  }

 private:
  /** \brief Recursive expansion of the simplex tree.*/
  void siblings_expansion(Siblings * siblings,  // must contain elements
                          int k) {
    if (dimension_ > k) {
      dimension_ = k;
    }
    if (k == 0)
      return;
    Dictionary_it next = siblings->members().begin();
    ++next;

#ifdef GUDHI_CAN_USE_CXX11_THREAD_LOCAL
    thread_local
#endif  // GUDHI_CAN_USE_CXX11_THREAD_LOCAL
    std::vector<std::pair<Vertex_handle, Node> > inter;
    for (Dictionary_it s_h = siblings->members().begin();
         s_h != siblings->members().end(); ++s_h, ++next) {
      Simplex_handle root_sh = find_vertex(s_h->first);
      if (has_children(root_sh)) {
        intersection(
                     inter,  // output intersection
                     next,  // begin
                     siblings->members().end(),  // end
                     root_sh->second.children()->members().begin(),
                     root_sh->second.children()->members().end(),
                     s_h->second.filtration());
        if (inter.size() != 0) {
          Siblings * new_sib = new Siblings(siblings,  // oncles
                                            s_h->first,  // parent
                                            inter);  // boost::container::ordered_unique_range_t
          inter.clear();
          s_h->second.assign_children(new_sib);
          siblings_expansion(new_sib, k - 1);
        } else {
          // ensure the children property
          s_h->second.assign_children(siblings);
          inter.clear();
        }
      }
    }
  }

  /** \brief Intersects Dictionary 1 [begin1;end1) with Dictionary 2 [begin2,end2)
   * and assigns the maximal possible Filtration_value to the Nodes. */
  static void intersection(std::vector<std::pair<Vertex_handle, Node> >& intersection,
                           Dictionary_it begin1, Dictionary_it end1,
                           Dictionary_it begin2, Dictionary_it end2,
                           Filtration_value filtration_) {
    if (begin1 == end1 || begin2 == end2)
      return;  // ----->>
    while (true) {
      if (begin1->first == begin2->first) {
        Filtration_value filt = (std::max)({begin1->second.filtration(), begin2->second.filtration(), filtration_});
        intersection.emplace_back(begin1->first, Node(nullptr, filt));
        if (++begin1 == end1 || ++begin2 == end2)
          return;  // ----->>
      } else if (begin1->first < begin2->first) {
        if (++begin1 == end1)
          return;
      } else /* begin1->first > begin2->first */ {
        if (++begin2 == end2)
          return;  // ----->>
      }
    }
  }

 public:
  /** \brief Expands a simplex tree containing only a graph. Simplices corresponding to cliques in the graph are added
   * incrementally, faces before cofaces, unless the simplex has dimension larger than `max_dim` or `block_simplex`
   * returns true for this simplex.
   *
   * @param[in] max_dim Expansion maximal dimension value.
   * @param[in] block_simplex Blocker oracle. Its concept is <CODE>bool block_simplex(Simplex_handle sh)</CODE>
   *
   * The function identifies a candidate simplex whose faces are all already in the complex, inserts
   * it with a filtration value corresponding to the maximum of the filtration values of the faces, then calls
   * `block_simplex` on a `Simplex_handle` for this new simplex. If `block_simplex` returns true, the simplex is
   * removed, otherwise it is kept. Note that the evaluation of `block_simplex` is a good time to update the
   * filtration value of the simplex if you want a customized value. The algorithm then proceeds with the next
   * candidate.
   *
   * @warning several candidates of the same dimension may be inserted simultaneously before calling `block_simplex`,
   * so if you examine the complex in `block_simplex`, you may hit a few simplices of the same dimension that have not
   * been vetted by `block_simplex` yet, or have already been rejected but not yet removed.
   */
  template< typename Blocker >
  void expansion_with_blockers(int max_dim, Blocker block_simplex) {
    // Loop must be from the end to the beginning, as higher dimension simplex are always on the left part of the tree
    for (auto& simplex : boost::adaptors::reverse(root_.members())) {
      if (has_children(&simplex)) {
        siblings_expansion_with_blockers(simplex.second.children(), max_dim, max_dim - 1, block_simplex);
      }
    }
  }

 private:
  /** \brief Recursive expansion with blockers of the simplex tree.*/
  template< typename Blocker >
  void siblings_expansion_with_blockers(Siblings* siblings, int max_dim, int k, Blocker block_simplex) {
    if (dimension_ < max_dim - k) {
      dimension_ = max_dim - k;
    }
    if (k == 0)
      return;
    // No need to go deeper
    if (siblings->members().size() < 2)
      return;
    // Reverse loop starting before the last one for 'next' to be the last one
    for (auto simplex = siblings->members().rbegin() + 1; simplex != siblings->members().rend(); simplex++) {
      std::vector<std::pair<Vertex_handle, Node> > intersection;
      for(auto next = siblings->members().rbegin(); next != simplex; next++) {
        bool to_be_inserted = true;
        Filtration_value filt = simplex->second.filtration();
        // If all the boundaries are present, 'next' needs to be inserted
        for (Simplex_handle border : boundary_simplex_range(simplex)) {
          Simplex_handle border_child = find_child(border, next->first);
          if (border_child == null_simplex()) {
            to_be_inserted=false;
            break;
          }
          filt = (std::max)(filt, filtration(border_child));
        }
        if (to_be_inserted) {
          intersection.emplace_back(next->first, Node(nullptr, filt));
        }
      }
      if (intersection.size() != 0) {
        // Reverse the order to insert
        Siblings * new_sib = new Siblings(siblings,  // oncles
                                          simplex->first,  // parent
                                          boost::adaptors::reverse(intersection));  // boost::container::ordered_unique_range_t
        std::vector<Vertex_handle> blocked_new_sib_vertex_list;
        // As all intersections are inserted, we can call the blocker function on all new_sib members
        for (auto new_sib_member = new_sib->members().begin();
             new_sib_member != new_sib->members().end();
             new_sib_member++) {
           bool blocker_result = block_simplex(new_sib_member);
           // new_sib member has been blocked by the blocker function
           // add it to the list to be removed - do not perform it while looping on it
           if (blocker_result) {
             blocked_new_sib_vertex_list.push_back(new_sib_member->first);
           }
        }
        if (blocked_new_sib_vertex_list.size() == new_sib->members().size()) {
          // Specific case where all have to be deleted
          delete new_sib;
          // ensure the children property
          simplex->second.assign_children(siblings);
        } else {
          for (auto& blocked_new_sib_member : blocked_new_sib_vertex_list) {
            new_sib->members().erase(blocked_new_sib_member);
          }
          // ensure recursive call
          simplex->second.assign_children(new_sib);
          siblings_expansion_with_blockers(new_sib, max_dim, k - 1, block_simplex);
        }
      } else {
        // ensure the children property
        simplex->second.assign_children(siblings);
      }
    }
  }

  /* \private Returns the Simplex_handle composed of the vertex list (from the Simplex_handle), plus the given
   * Vertex_handle if the Vertex_handle is found in the Simplex_handle children list.
   * Returns null_simplex() if it does not exist
  */
  Simplex_handle find_child(Simplex_handle sh, Vertex_handle vh) const {
    if (!has_children(sh))
      return null_simplex();

    Simplex_handle child = sh->second.children()->find(vh);
    // Specific case of boost::flat_map does not find, returns boost::flat_map::end()
    // in simplex tree we want a null_simplex()
    if (child == sh->second.children()->members().end())
      return null_simplex();

    return child;
  }

 public:
  /** \brief Write the hasse diagram of the simplicial complex in os.
   *
   * Each row in the file correspond to a simplex. A line is written:
   * dim idx_1 ... idx_k fil   where dim is the dimension of the simplex,
   * idx_1 ... idx_k are the row index (starting from 0) of the simplices of the boundary
   * of the simplex, and fil is its filtration value. */
  void print_hasse(std::ostream& os) {
    os << num_simplices() << " " << std::endl;
    for (auto sh : filtration_simplex_range()) {
      os << dimension(sh) << " ";
      for (auto b_sh : boundary_simplex_range(sh)) {
        os << key(b_sh) << " ";
      }
      os << filtration(sh) << " \n";
    }
  }

 public:
  /** \brief This function ensures that each simplex has a higher filtration value than its faces by increasing the
   * filtration values.
   * @return True if any filtration value was modified, false if the filtration was already non-decreasing.
   * \post Some simplex tree functions require the filtration to be valid. `make_filtration_non_decreasing()`
   * function is not launching `initialize_filtration()` but returns the filtration modification information. If the
   * complex has changed , please call `initialize_filtration()` to recompute it.
   */
  bool make_filtration_non_decreasing() {
    bool modified = false;
    // Loop must be from the end to the beginning, as higher dimension simplex are always on the left part of the tree
    for (auto& simplex : boost::adaptors::reverse(root_.members())) {
      if (has_children(&simplex)) {
        modified |= rec_make_filtration_non_decreasing(simplex.second.children());
      }
    }
    return modified;
  }

 private:
  /** \brief Recursively Browse the simplex tree to ensure the filtration is not decreasing.
   * @param[in] sib Siblings to be parsed.
   * @return The filtration modification information in order to trigger initialize_filtration.
   */
  bool rec_make_filtration_non_decreasing(Siblings * sib) {
    bool modified = false;

    // Loop must be from the end to the beginning, as higher dimension simplex are always on the left part of the tree
    for (auto& simplex : boost::adaptors::reverse(sib->members())) {
      // Find the maximum filtration value in the border
      Boundary_simplex_range boundary = boundary_simplex_range(&simplex);
      Boundary_simplex_iterator max_border = std::max_element(std::begin(boundary), std::end(boundary),
                                                              [](Simplex_handle sh1, Simplex_handle sh2) {
                                                                return filtration(sh1) < filtration(sh2);
                                                              });

      Filtration_value max_filt_border_value = filtration(*max_border);
      if (simplex.second.filtration() < max_filt_border_value) {
        // Store the filtration modification information
        modified = true;
        simplex.second.assign_filtration(max_filt_border_value);
      }
      if (has_children(&simplex)) {
        modified |= rec_make_filtration_non_decreasing(simplex.second.children());
      }
    }
    // Make the modified information to be traced by upper call
    return modified;
  }

 public:
  /** \brief Prune above filtration value given as parameter.
   * @param[in] filtration Maximum threshold value.
   * @return The filtration modification information.
   * \post Some simplex tree functions require the filtration to be valid. `prune_above_filtration()`
   * function is not launching `initialize_filtration()` but returns the filtration modification information. If the
   * complex has changed , please call `initialize_filtration()` to recompute it.
   * \post Note that the dimension of the simplicial complex may be lower after calling `prune_above_filtration()`
   * than it was before. However, `upper_bound_dimension()` will return the old value, which remains a valid upper
   * bound. If you care, you can call `dimension()` to recompute the exact dimension.
   */
  bool prune_above_filtration(Filtration_value filtration) {
    return rec_prune_above_filtration(root(), filtration);
  }

 private:
  bool rec_prune_above_filtration(Siblings* sib, Filtration_value filt) {
    auto&& list = sib->members();
    auto last = std::remove_if(list.begin(), list.end(), [=](Dit_value_t& simplex) {
        if (simplex.second.filtration() <= filt) return false;
        if (has_children(&simplex)) rec_delete(simplex.second.children());
        // dimension may need to be lowered
        dimension_to_be_lowered_ = true;
        return true;
      });

    bool modified = (last != list.end());
    if (last == list.begin() && sib != root()) {
      // Removing the whole siblings, parent becomes a leaf.
      sib->oncles()->members()[sib->parent()].assign_children(sib->oncles());
      delete sib;
      // dimension may need to be lowered
      dimension_to_be_lowered_ = true;
      return true;
    } else {
      // Keeping some elements of siblings. Remove the others, and recurse in the remaining ones.
      list.erase(last, list.end());
      for (auto&& simplex : list)
        if (has_children(&simplex))
          modified |= rec_prune_above_filtration(simplex.second.children(), filt);
    }
    return modified;
  }

 private:
  /** \brief Deep search simplex tree dimension recompute.
   * @return True if the dimension was modified, false otherwise.
   * \pre Be sure the simplex tree has not a too low dimension value as the deep search stops when the former dimension
   * has been reached (cf. `upper_bound_dimension()` and `set_dimension()` methods).
   */
  bool lower_upper_bound_dimension() {
    // reset automatic detection to recompute
    dimension_to_be_lowered_ = false;
    int new_dimension = -1;
    // Browse the tree from the left to the right as higher dimension cells are more likely on the left part of the tree
    for (Simplex_handle sh : complex_simplex_range()) {
#ifdef DEBUG_TRACES
      for (auto vertex : simplex_vertex_range(sh)) {
        std::cout << " " << vertex;
      }
      std::cout << std::endl;
#endif  // DEBUG_TRACES

      int sh_dimension = dimension(sh);
      if (sh_dimension >= dimension_)
        // Stop browsing as soon as the dimension is reached, no need to go furter
        return false;
      new_dimension = (std::max)(new_dimension, sh_dimension);
    }
    dimension_ = new_dimension;
    return true;
  }


 public:
  /** \brief Remove a maximal simplex.
   * @param[in] sh Simplex handle on the maximal simplex to remove.
   * \pre Please check the simplex has no coface before removing it.
   * \exception std::invalid_argument In debug mode, if sh has children.
   * \post Be aware that removing is shifting data in a flat_map (initialize_filtration to be done).
   * \post Note that the dimension of the simplicial complex may be lower after calling `remove_maximal_simplex()`
   * than it was before. However, `upper_bound_dimension()` will return the old value, which remains a valid upper
   * bound. If you care, you can call `dimension()` to recompute the exact dimension.
   */
  void remove_maximal_simplex(Simplex_handle sh) {
    // Guarantee the simplex has no children
    GUDHI_CHECK(!has_children(sh),
                std::invalid_argument("Simplex_tree::remove_maximal_simplex - argument has children"));

    // Simplex is a leaf, it means the child is the Siblings owning the leaf
    Siblings* child = sh->second.children();

    if ((child->size() > 1) || (child == root())) {
      // Not alone, just remove it from members
      // Special case when child is the root of the simplex tree, just remove it from members
      child->erase(sh);
    } else {
      // Sibling is emptied : must be deleted, and its parent must point on his own Sibling
      child->oncles()->members().at(child->parent()).assign_children(child->oncles());
      delete child;
      // dimension may need to be lowered
      dimension_to_be_lowered_ = true;
    }
  }

 private:
  Vertex_handle null_vertex_;
  /** \brief Total number of simplices in the complex, without the empty simplex.*/
  /** \brief Set of simplex tree Nodes representing the vertices.*/
  Siblings root_;
  /** \brief Simplices ordered according to a filtration.*/
  std::vector<Simplex_handle> filtration_vect_;
  /** \brief Upper bound on the dimension of the simplicial complex.*/
  int dimension_;
  bool dimension_to_be_lowered_ = false;
};

// Print a Simplex_tree in os.
template<typename...T>
std::ostream& operator<<(std::ostream & os, Simplex_tree<T...> & st) {
  for (auto sh : st.filtration_simplex_range()) {
    os << st.dimension(sh) << " ";
    for (auto v : st.simplex_vertex_range(sh)) {
      os << v << " ";
    }
    os << st.filtration(sh) << "\n";  // TODO(VR): why adding the key ?? not read ?? << "     " << st.key(sh) << " \n";
  }
  return os;
}

template<typename...T>
std::istream& operator>>(std::istream & is, Simplex_tree<T...> & st) {
  typedef Simplex_tree<T...> ST;
  std::vector<typename ST::Vertex_handle> simplex;
  typename ST::Filtration_value fil;
  int max_dim = -1;
  while (read_simplex(is, simplex, fil)) {
    // read all simplices in the file as a list of vertices
    // Warning : simplex_size needs to be casted in int - Can be 0
    int dim = static_cast<int> (simplex.size() - 1);
    if (max_dim < dim) {
      max_dim = dim;
    }
    // insert every simplex in the simplex tree
    st.insert_simplex(simplex, fil);
    simplex.clear();
  }
  st.set_dimension(max_dim);

  return is;
}

/** Model of SimplexTreeOptions.
 * 
 * Maximum number of simplices to compute persistence is <CODE>std::numeric_limits<std::uint32_t>::max()</CODE>
 * (about 4 billions of simplices). */
struct Simplex_tree_options_full_featured {
  typedef linear_indexing_tag Indexing_tag;
  typedef int Vertex_handle;
  typedef double Filtration_value;
  typedef std::uint32_t Simplex_key;
  static const bool store_key = true;
  static const bool store_filtration = true;
  static const bool contiguous_vertices = false;
};

/** Model of SimplexTreeOptions, faster than `Simplex_tree_options_full_featured` but note the unsafe
 * `contiguous_vertices` option.
 * 
 * Maximum number of simplices to compute persistence is <CODE>std::numeric_limits<std::uint32_t>::max()</CODE>
 * (about 4 billions of simplices). */

struct Simplex_tree_options_fast_persistence {
  typedef linear_indexing_tag Indexing_tag;
  typedef int Vertex_handle;
  typedef float Filtration_value;
  typedef std::uint32_t Simplex_key;
  static const bool store_key = true;
  static const bool store_filtration = true;
  static const bool contiguous_vertices = true;
};

/** @} */  // end defgroup simplex_tree

}  // namespace Gudhi

#endif  // SIMPLEX_TREE_H_
