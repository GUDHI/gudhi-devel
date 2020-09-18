/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
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
#include <gudhi/Simplex_tree/Simplex_tree_zigzag_iterators.h>
#include <gudhi/Simplex_tree/indexing_tag.h>

#include <gudhi/reader_utils.h>
#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/Debug_utils.h>

#include <boost/container/map.hpp>
#include <boost/container/flat_map.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/container/static_vector.hpp>

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

/**
 * \class Extended_simplex_type Simplex_tree.h gudhi/Simplex_tree.h
 * \brief Extended simplex type data structure for representing the type of simplices in an extended filtration.
 *
 * \details The extended simplex type can be either UP (which means
 * that the simplex was present originally, and is thus part of the ascending extended filtration), DOWN (which means
 * that the simplex is the cone of an original simplex, and is thus part of the descending extended filtration) or
 * EXTRA (which means the simplex is the cone point).
 *
 * Details may be found in \cite Cohen-Steiner2009 and section 2.2 in \cite Carriere16.
 *
 */
enum class Extended_simplex_type {UP, DOWN, EXTRA};

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
  // typedef typename boost::container::flat_map<Vertex_handle, Node> Dictionary;

  typedef typename boost::container::flat_map<Vertex_handle, Node> flat_map;
  //Dictionary::iterator remain valid under insertions and deletions
  typedef typename boost::container::map<Vertex_handle, Node> map;
  typedef typename std::conditional<Options::simplex_handle_strong_validity, map, 
                                    flat_map>::type Dictionary;

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
  struct Extended_filtration_data {
    Filtration_value minval;
    Filtration_value maxval;
    Extended_filtration_data(){}
    Extended_filtration_data(Filtration_value vmin, Filtration_value vmax): minval(vmin), maxval(vmax) {}
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
   * by the simplex tree.
   *
   * They are essentially pointers into internal vectors, and any insertion or removal
   * of a simplex may invalidate any other Simplex_handle in the complex. */
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
    /** Forward iterator on the simplices (insertion and deletion) of a zigzag
      * filtration.
      *
      * 'value_type' is Simplex_handle.
      */
    typedef Flagzigzag_simplex_iterator< Simplex_tree >
                                                  Zigzagfiltration_simplex_iterator;
    /** Range for the flag zigzag filtration.*/
    typedef boost::iterator_range< Zigzagfiltration_simplex_iterator >
                                                     Zigzagfiltration_simplex_range;
    /** \brief Range over the simplices of the simplicial complex, ordered by the 
      * filtration. */
    typedef typename std::conditional< Options::is_zigzag, 
                      Zigzagfiltration_simplex_range,
                      std::vector<Simplex_handle> >::type  Filtration_simplex_range;
    /** \brief Iterator over the simplices of the simplicial complex, ordered by 
      * the filtration.
      *
      * 'value_type' is Simplex_handle. */
    typedef typename std::conditional< Options::is_zigzag,
          Zigzagfiltration_simplex_iterator,
          typename std::vector<Simplex_handle>::const_iterator >::type
                                                        Filtration_simplex_iterator;

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
   * was initialized, please call `clear_filtration()` or `initialize_filtration()` to recompute it. */
  // Filtration_simplex_range const& filtration_simplex_range(Indexing_tag = Indexing_tag()) {
  Filtration_simplex_range const& filtration_simplex_range(linear_indexing_tag) {
    maybe_initialize_filtration();
    return filtration_vect_;
  }

  //Call the appropriate filtration_simplex_range depending on whether the indexing 
  //tag is linear_indexing_tag or zigzag_indexing_tag  
  Filtration_simplex_range & filtration_simplex_range()
  { return filtration_simplex_range(Indexing_tag()); }

  /** \brief Returns a range over the vertices of a simplex.
   *
   * The order in which the vertices are visited is the decreasing order for < on Vertex_handles,
   * which is consequenlty
   * equal to \f$(-1)^{\text{dim} \sigma}\f$ the canonical orientation on the simplex.
   */
  Simplex_vertex_range simplex_vertex_range(Simplex_handle sh) const {
    GUDHI_CHECK(sh != null_simplex(), "empty simplex");
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
      dimension_(-1),
      zigzag_simplex_range_initialized_(false) { }

  /** \brief User-defined copy constructor reproduces the whole tree structure. */
  Simplex_tree(const Simplex_tree& complex_source) {
#ifdef DEBUG_TRACES
    std::clog << "Simplex_tree copy constructor" << std::endl;
#endif  // DEBUG_TRACES
    copy_from(complex_source);
  }

  /** \brief User-defined move constructor relocates the whole tree structure.
   *  \exception std::invalid_argument In debug mode, if the complex_source is invalid.
   */
  Simplex_tree(Simplex_tree && complex_source) {
#ifdef DEBUG_TRACES
    std::clog << "Simplex_tree move constructor" << std::endl;
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
    std::clog << "Simplex_tree copy assignment" << std::endl;
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
    std::clog << "Simplex_tree move assignment" << std::endl;
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

  /** \brief Returns the filtration value of a simplex.
   *
   * Same as `filtration()`, but does not handle `null_simplex()`.
   */
  static Filtration_value filtration_(Simplex_handle sh) {
    GUDHI_CHECK (sh != null_simplex(), "null simplex");
    return sh->second.filtration();
  }

 public:
  /** \brief Returns the key associated to a simplex.
   *
   * If no key has been assigned, returns `null_key()`.
   * \pre SimplexTreeOptions::store_key
   */
  static Simplex_key key(Simplex_handle sh) {
    return sh->second.key();
  }

  /** \brief Returns the simplex that has index idx in the filtration.
   *
   * The filtration must be initialized.
   */
  Simplex_handle simplex(Simplex_key idx) const {
    simplex(idx, Indexing_tag());
  }
  Simplex_handle simplex(Simplex_key idx, linear_indexing_tag) const {
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

  /** \brief Returns a fixed number not in the interval [0, `num_simplices()`).  */
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

    thread_local std::vector<Vertex_handle> copy;
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
  static Siblings* self_siblings(SimplexHandle sh) {
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
  /** \brief Initializes the filtration cache, i.e. sorts the
   * simplices according to their order in the filtration.
   *
   * It always recomputes the cache, even if one already exists.
   *
   * Any insertion, deletion or change of filtration value invalidates this cache,
   * which can be cleared with clear_filtration().  */
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
  /** \brief Initializes the filtration cache if it isn't initialized yet.
   *
   * Automatically called by filtration_simplex_range(). */
  void maybe_initialize_filtration() {
    if (filtration_vect_.empty()) {
      initialize_filtration();
    }
  }
  /** \brief Clears the filtration cache produced by initialize_filtration().
   *
   * Useful when initialize_filtration() has already been called and we perform an operation
   * (say an insertion) that invalidates the cache. */
  void clear_filtration() {
    filtration_vect_.clear();
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
    clear_filtration(); // Drop the cache.
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

    thread_local std::vector<std::pair<Vertex_handle, Node> > inter;
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
   * 
   * If a simplex has a `NaN` filtration value, it is considered lower than any other defined filtration value.
   */
  bool make_filtration_non_decreasing() {
    bool modified = false;
    // Loop must be from the end to the beginning, as higher dimension simplex are always on the left part of the tree
    for (auto& simplex : boost::adaptors::reverse(root_.members())) {
      if (has_children(&simplex)) {
        modified |= rec_make_filtration_non_decreasing(simplex.second.children());
      }
    }
    if(modified)
      clear_filtration(); // Drop the cache.
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
      // Replacing if(f<max) with if(!(f>=max)) would mean that if f is NaN, we replace it with the max of the children.
      // That seems more useful than keeping NaN.
      if (!(simplex.second.filtration() >= max_filt_border_value)) {
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
   * @return True if any simplex was removed, false if all simplices already had a value below the threshold.
   * \post Note that the dimension of the simplicial complex may be lower after calling `prune_above_filtration()`
   * than it was before. However, `upper_bound_dimension()` will return the old value, which remains a valid upper
   * bound. If you care, you can call `dimension()` to recompute the exact dimension.
   */
  bool prune_above_filtration(Filtration_value filtration) {
    bool modified = rec_prune_above_filtration(root(), filtration);
    if(modified)
      clear_filtration(); // Drop the cache.
    return modified;
  }

 private:
  bool rec_prune_above_filtration(Siblings* sib, Filtration_value filt) {
    auto&& list = sib->members();
    auto last = std::remove_if(list.begin(), list.end(), [this,filt](Dit_value_t& simplex) {
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
        std::clog << " " << vertex;
      }
      std::clog << std::endl;
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

  /** \brief Retrieve the original filtration value for a given simplex in the Simplex_tree. Since the 
   * computation of extended persistence requires modifying the filtration values, this function can be used
   * to recover the original values. Moreover, computing extended persistence requires adding new simplices
   * in the Simplex_tree. Hence, this function also outputs the type of each simplex. It can be either UP (which means
   * that the simplex was present originally, and is thus part of the ascending extended filtration), DOWN (which means
   * that the simplex is the cone of an original simplex, and is thus part of the descending extended filtration) or
   * EXTRA (which means the simplex is the cone point). See the definition of Extended_simplex_type. Note that if the simplex type is DOWN, the original filtration value
   * is set to be the original filtration value of the corresponding (not coned) original simplex. 
   * \pre This function should be called only if `extend_filtration()` has been called first!
   * \post The output filtration value is supposed to be the same, but might be a little different, than the
   * original filtration value, due to the internal transformation (scaling to [-2,-1]) that is 
   * performed on the original filtration values during the computation of extended persistence.
   * @param[in] f Filtration value of the simplex in the extended (i.e., modified) filtration.
   * @param[in] efd Structure containing the minimum and maximum values of the original filtration. This the output of `extend_filtration()`.
   * @return A pair containing the original filtration value of the simplex as well as the simplex type.
   */
  std::pair<Filtration_value, Extended_simplex_type> decode_extended_filtration(Filtration_value f, const Extended_filtration_data& efd){
    std::pair<Filtration_value, Extended_simplex_type> p;
    Filtration_value minval = efd.minval;
    Filtration_value maxval = efd.maxval;
    if (f >= -2 && f <= -1){
      p.first = minval + (maxval-minval)*(f + 2); p.second = Extended_simplex_type::UP;
    }
    else if (f >= 1 && f <= 2){
      p.first = minval - (maxval-minval)*(f - 2); p.second = Extended_simplex_type::DOWN;
    }
    else{
      p.first = std::numeric_limits<Filtration_value>::quiet_NaN(); p.second = Extended_simplex_type::EXTRA;
    }
    return p;
  };

  /** \brief Extend filtration for computing extended persistence. 
   * This function only uses the filtration values at the 0-dimensional simplices, 
   * and computes the extended persistence diagram induced by the lower-star filtration 
   * computed with these values. 
   * \post Note that after calling this function, the filtration 
   * values are actually modified. The function `decode_extended_filtration()` 
   * retrieves the original values and outputs the extended simplex type.
   * \pre Note that this code creates an extra vertex internally, so you should make sure that
   * the Simplex tree does not contain a vertex with the largest Vertex_handle.
   * @return A data structure containing the maximum and minimum values of the original filtration.
   * It is meant to be provided as input to `decode_extended_filtration()` in order to retrieve
   * the original filtration values for each simplex.
   */
  Extended_filtration_data extend_filtration() {
    clear_filtration(); // Drop the cache.

    // Compute maximum and minimum of filtration values
    Vertex_handle maxvert = std::numeric_limits<Vertex_handle>::min();
    Filtration_value minval = std::numeric_limits<Filtration_value>::infinity();
    Filtration_value maxval = -std::numeric_limits<Filtration_value>::infinity();
    for (auto sh = root_.members().begin(); sh != root_.members().end(); ++sh){
      Filtration_value f = this->filtration(sh);
      minval = std::min(minval, f);
      maxval = std::max(maxval, f);
      maxvert = std::max(sh->first, maxvert);
    }
    
    GUDHI_CHECK(maxvert < std::numeric_limits<Vertex_handle>::max(), std::invalid_argument("Simplex_tree contains a vertex with the largest Vertex_handle"));
    maxvert += 1;

    Simplex_tree st_copy = *this;

    // Add point for coning the simplicial complex
    this->insert_simplex({maxvert}, -3);

    // For each simplex
    std::vector<Vertex_handle> vr;
    for (auto sh_copy : st_copy.complex_simplex_range()){

      // Locate simplex
      vr.clear();
      for (auto vh : st_copy.simplex_vertex_range(sh_copy)){
        vr.push_back(vh);
      }
      auto sh = this->find(vr);

      // Create cone on simplex
      vr.push_back(maxvert);
      if (this->dimension(sh) == 0){
        Filtration_value v = this->filtration(sh);
        Filtration_value scaled_v = (v-minval)/(maxval-minval);
        // Assign ascending value between -2 and -1 to vertex
        this->assign_filtration(sh, -2 + scaled_v);
        // Assign descending value between 1 and 2 to cone on vertex
        this->insert_simplex(vr, 2 - scaled_v);
      }
      else{
        // Assign value -3 to simplex and cone on simplex
        this->assign_filtration(sh, -3);
        this->insert_simplex(vr, -3);
      }
    }

    // Automatically assign good values for simplices
    this->make_filtration_non_decreasing();

    // Return the filtration data 
    Extended_filtration_data efd(minval, maxval);
    return efd;
  }

  /** \brief Returns a vertex of `sh` that has the same filtration value as `sh` if it exists, and `null_vertex()` otherwise.
   *
   * For a lower-star filtration built with `make_filtration_non_decreasing()`, this is a way to invert the process and find out which vertex had its filtration value propagated to `sh`.
   * If several vertices have the same filtration value, the one it returns is arbitrary. */
  Vertex_handle vertex_with_same_filtration(Simplex_handle sh) {
    auto filt = filtration_(sh);
    for(auto v : simplex_vertex_range(sh))
      if(filtration_(find_vertex(v)) == filt)
        return v;
    return null_vertex();
  }

  /** \brief Returns an edge of `sh` that has the same filtration value as `sh` if it exists, and `null_simplex()` otherwise.
   *
   * For a flag-complex built with `expansion()`, this is a way to invert the process and find out which edge had its filtration value propagated to `sh`.
   * If several edges have the same filtration value, the one it returns is arbitrary.
   *
   * \pre `sh` must have dimension at least 1. */
  Simplex_handle edge_with_same_filtration(Simplex_handle sh) {
    // See issue #251 for potential speed improvements.
    auto&& vertices = simplex_vertex_range(sh); // vertices in decreasing order
    auto end = std::end(vertices);
    auto vi = std::begin(vertices);
    GUDHI_CHECK(vi != end, "empty simplex");
    auto v0 = *vi;
    ++vi;
    GUDHI_CHECK(vi != end, "simplex of dimension 0");
    if(std::next(vi) == end) return sh; // shortcut for dimension 1
    boost::container::static_vector<Vertex_handle, 40> suffix;
    suffix.push_back(v0);
    auto filt = filtration_(sh);
    do
    {
      Vertex_handle v = *vi;
      auto&& children1 = find_vertex(v)->second.children()->members_;
      for(auto w : suffix){
        // Can we take advantage of the fact that suffix is ordered?
        Simplex_handle s = children1.find(w);
        if(filtration_(s) == filt)
          return s;
      }
      suffix.push_back(v);
    }
    while(++vi != end);
    return null_simplex();
  }

  /** \brief Returns a minimal face of `sh` that has the same filtration value as `sh`.
   *
   * For a filtration built with `make_filtration_non_decreasing()`, this is a way to invert the process and find out which simplex had its filtration value propagated to `sh`.
   * If several minimal (for inclusion) simplices have the same filtration value, the one it returns is arbitrary, and it is not guaranteed to be the one with smallest dimension. */
  Simplex_handle minimal_simplex_with_same_filtration(Simplex_handle sh) {
    auto filt = filtration_(sh);
    // Naive implementation, it can be sped up.
    for(auto b : boundary_simplex_range(sh))
      if(filtration_(b) == filt)
        return minimal_simplex_with_same_filtration(b);
    return sh; // None of its faces has the same filtration.
  }



public:

/* Dynamic flag complex.
 *
 * The following methods are dedicated to a stream-like construction of
 * of flag complexes, with addition and removal of vertices and edges as
 * atomic operations.
 *
 * This can be used for a stream-like construction of zigzag filtrations (of flag
 * complexes) when knowing only the sequence of insertions and deletions of
 * vertices and edges.
 *
 * Options::store_key
 * ????, Options::contiguous_vertices,
 * ????  Options::link_simplices_through_max_vertex              must all be true.
 */

private:
/**
  * Returns a range over the simplices of the flag zigzag filtration encoded
  * for a vector of insertion and deletion of edges.
  *
  * @param[in] zz_edge_fil         range of Zigzag_edge< Simplex_tree >,
  *                                containing a valid sequence of insertions and
  *                                removal of edges.
  *
  * @param[in] dim_max             maximal dimension of the dynamic complex
  *                                constructed.
  *
  * A ZigzagEdgeRange must be a range of ZigzagEdge. A model of ZigzagEdge
  * must contain operations:
  * Vertex_handle u() return the endpoint with smaller label,
  * Vertex_handle v() return the endpoint with bigger label,
  * Filtration_value fil() return the filtration value in the zigzag
  * filtration,
  * bool type() return true if the edge is inserted, false if it is removed.
  * if u() == v(), this is a vertex.
  *
  * ZigzagEdge must be Zigzag_edge< Simplex_tree >.
  */
public:
/** Returns the zigzag simplex range for the (dim_max)-skeleton of a flag zigzag 
  * filtration represented implicitely by the sequence of insertions and 
  * deletions of vertices and edges in zz_edge_fil.
  */
  template< class ZigzagEdgeRange >
  Zigzagfiltration_simplex_range
  zigzag_simplex_range( ZigzagEdgeRange &zz_edge_fil, int dim_max )
  {
    return
      Zigzagfiltration_simplex_range(
          Zigzagfiltration_simplex_iterator(this, &zz_edge_fil, dim_max)
          , Zigzagfiltration_simplex_iterator()  );
  }
/** Returns the zigzag simplex range for the (dim_max)-skeleton of the 
  * oscillating Rips zigzag filtration with paramters mu and nu, built on top 
  * of the point cloud points with distance function distance. 
  */
  template< typename Kernel, typename Point_container, typename Distance >
  Zigzagfiltration_simplex_range
  zigzag_simplex_range( Point_container const        & points,
                        Distance        const          distance,
                        Filtration_value const          nu,
                        Filtration_value const          mu,
                        int                            dim_max )
  {
    return
      Zigzagfiltration_simplex_range(
          Zigzagfiltration_simplex_iterator(this, points, distance, nu, mu, dim_max) , Zigzagfiltration_simplex_iterator()  );
  }

public:
  //Initialize a Flag_zigzag_simplex_range
  template< class ZigzagEdgeRange >
  void initialize_filtration( ZigzagEdgeRange & zz_edge_fil, int dim_max )
  { //empty complex
    //todo
    // #ifdef COMPLEX_TIME
    //   start_ = std::chrono::system_clock::now();
    // #endif
    zigzag_simplex_range_ = zigzag_simplex_range(zz_edge_fil, dim_max);
    zigzag_simplex_range_initialized_ = true;
    // #ifdef COMPLEX_TIME
    //   end_ = std::chrono::system_clock::now();
    //   enlapsed_sec_ = end_ - start_;
    //   complexOperationTime_ += enlapsed_sec_.count();
    // #endif
  }

  template<typename Point_container, typename Distance>
  void initialize_filtration( Point_container const        & points,
                              Distance        const          distance,
                              Filtration_value const          nu,
                              Filtration_value const          mu,
                              int                            dim_max ) 
  {
    zigzag_simplex_range_ = zigzag_simplex_range(points, distance, nu, mu, dim_max);
    zigzag_simplex_range_initialized_ = true;
  }

  //must call initialize_filtration beforehand
  Zigzagfiltration_simplex_range & filtration_simplex_range(zigzag_indexing_tag)
  {
    assert(zigzag_simplex_range_initialized_);
    zigzag_simplex_range_initialized_ = false;
    return zigzag_simplex_range_;
  }

//CONTINUE FROM HERE

public:

/**
* Add an edge in the complex, its two vertices (if missing)
* and all the missing
* simplices of its star, defined to maintain the property
* of the complex to be a flag complex. The edge {u,v} must not be in the complex.
*
* In term of edges in the graph, inserting edge u,v only affects N^+(u),
* even if u and v were missing.
*
* For a new node with label v, we first do a local expansion for
* computing the
* children of this new node, and then a standard expansion for its children.
* Nodes with label v (and their subtrees) already in the tree
* do not get affected.
*
* Nodes with label u get affected only if a Node with label v is in their same
* siblings set.
* We then try to insert "ponctually" v all over the subtree rooted
* at Node(u). Each
* insertion of a Node with v label induces a local expansion at this
* Node (as explained above)
* and a sequence of "ponctual" insertion of Node(v) in the subtree
* rooted at sibling nodes of the new node, on its left.
*
* @param[in] u,v              Vertex_handle representing the new edge
* @param[in] zz_filtration    Must be empty. Contains at the end all new
*                             simplices induced by the insertion of the edge.
*
* SimplexTreeOptions::link_simplices_through_max_vertex must be true.
* Simplex_tree::Dictionary must sort Vertex_handles w/ increasing natural order <
*/
  void flag_add_edge( Vertex_handle                   u
                    , Vertex_handle                   v
                    , Filtration_value                fil
                    , int                             dim_max
                    , std::vector< Simplex_handle > & zz_filtration )
  {
// #ifdef COMPLEX_TIME
//   start_ = std::chrono::system_clock::now();
// #endif
  if(u == v) { // Are we inserting a vertex?
      auto res_ins = root_.members().emplace(u,Node(&root_,fil));
      if(res_ins.second) { //if the vertex was not in the complex
    update_simplex_tree_after_node_insertion(res_ins.first);
    zz_filtration.push_back(res_ins.first); //no more insert in root_.members()
      }
// #ifdef COMPLEX_TIME
//       end_ = std::chrono::system_clock::now();
//       enlapsed_sec_ = end_ - start_;
//       complexOperationTime_ += enlapsed_sec_.count();
// #endif
      return; //because the vertex is isolated, no more insertions.
  }
  // else, we are inserting an edge: ensure that u < v
  if(v < u) { std::swap(u,v); }

  //Note that we copy Simplex_handle (aka map iterators) in zz_filtration
  //whereas we are still modifying the Simplex_tree. Insertions in siblings may
  //invalidate Simplex_handles; we take care of this fact by first doing all
  //insertion in a Sibling, then inserting all handles in zz_filtration.

  //check whether vertices u and v are in the tree, insert them if necessary
  auto res_ins_v = root_.members().emplace(v,Node(&root_,fil));
  auto res_ins_u = root_.members().emplace(u,Node(&root_,fil));
  if(res_ins_v.second) {
      // update_simplex_tree_after_node_insertion(res_ins_v.first);
      zz_filtration.push_back(res_ins_v.first); //no more inserts in root_.members()
  }
  if(res_ins_u.second) {
      // update_simplex_tree_after_node_insertion(res_ins_u.first);
      zz_filtration.push_back(res_ins_u.first); //no more inserts in root_.members()
  } //check if the edge {u,v} is already in the complex, if true, nothing to do.
  if(has_children(res_ins_u.first) && res_ins_u.first->second.children()->members().find(v) != res_ins_u.first->second.children()->members().end()) {
// #ifdef COMPLEX_TIME
//       end_ = std::chrono::system_clock::now();
//       enlapsed_sec_ = end_ - start_;
//       complexOperationTime_ += enlapsed_sec_.count();
// #endif
      return;
  }

  //upper bound on dimension
  dimension_ = dim_max; dimension_to_be_lowered_ = true;

  //for all siblings containing a Node labeled with u (including the root), run
  //a zz_punctual_expansion       //todo parallelise
  auto list_u_ptr = cofaces_data_structure_.access(u);//list of all u Nodes
  for( auto hook_u_it =  list_u_ptr->begin();
       hook_u_it != list_u_ptr->end();    ++hook_u_it )
  {
      Node & node_u    = static_cast<Node&>(*hook_u_it);//corresponding node
      Siblings * sib_u = self_siblings(node_u, u);//Siblings containing node
      if(sib_u->members().find(v) != sib_u->members().end()) {
    int curr_dim = dimension(node_u,u);
    if(curr_dim < dim_max)
    {
        if(!has_children(node_u, u)) //now has a new child Node labeled v
        { node_u.assign_children(new Siblings(sib_u, u)); }
        zz_punctual_expansion( v
             , node_u.children()
             , fil
             , dim_max - curr_dim -1 //>= 0
             , zz_filtration ); //u on top
    }
      }
  }
  //todo sort zz_filtration appropriately
  // sort( zz_filtration.begin(), zz_filtration.end(),
  //       is_before_in_filtration(this)
  //       );
#ifdef GUDHI_USE_TBB
  tbb::parallel_sort(zz_filtration.begin(), zz_filtration.end(), is_before_in_filtration(this));
#else
  sort(zz_filtration.begin(), zz_filtration.end(), is_before_in_filtration(this));
#endif
  //update cofaces data structures
  for(auto sh : zz_filtration) { update_simplex_tree_after_node_insertion(sh); }
// #ifdef COMPLEX_TIME
//   end_ = std::chrono::system_clock::now();
//   enlapsed_sec_ = end_ - start_;
//   complexOperationTime_ += enlapsed_sec_.count();
// #endif
    }

private:
    /*
 * Insert a Node with label v in the set of siblings sib, and percolate the
 * expansion on the subtree rooted at sib. Sibling sib must not contain
 * v.
 * The percolation of the expansion is twofold:
 * 1- the newly inserted Node labeled v in sib has a subtree computed
 * via zz_local_expansion.
 * 2- All Node in the members of sib, with label x and x < v,
 * need in turn a local_expansion by v iff N^+(x) contains v.
 */
    void zz_punctual_expansion( Vertex_handle    v
        , Siblings *       sib
        , Filtration_value fil
        , int              k //k == dim_max - dimension simplices in sib
        , std::vector<Simplex_handle> & zz_filtration )
    { //insertion always succeeds because the edge {u,v} used to not be here.
  auto res_ins_v = sib->members().emplace(v, Node(sib,fil));
  // update_simplex_tree_after_node_insertion(res_ins_v.first);//for cofaces hooks
  zz_filtration.push_back(res_ins_v.first); //no more insertion in sib

  if(k == 0) { return; } //reached the maximal dimension

  //create the subtree of new Node(v)
  zz_local_expansion( res_ins_v.first
          , sib
          , fil
          , k
          , zz_filtration );

  //punctual expansion in nodes on the left of v, i.e. with label x < v
  for( auto sh = sib->members().begin(); sh != res_ins_v.first; ++sh )
  { //if v belongs to N^+(x), punctual expansion
      Simplex_handle root_sh = find_vertex(sh->first); //Node(x), x < v
      if( has_children(root_sh) &&
        root_sh->second.children()->members().find(v)
        != root_sh->second.children()->members().end()  )
      { //edge {x,v} is in the complex
    if(!has_children(sh))
    { sh->second.assign_children(new Siblings(sib, sh->first)); }
    //insert v in the children of sh, and expand.
    zz_punctual_expansion( v
               , sh->second.children()
               , fil
               , k-1
               , zz_filtration );
      }
  }
    }

    /* After the insertion of edge {u,v}, expansion of a subtree rooted at v, where the
 * Node with label v has just been inserted, and its parent is a Node labeled with
 * u. sh has no children here.
 *
 * k must be > 0
 */
    void zz_local_expansion(
      Simplex_handle   sh_v    //Node with label v which has just been inserted
      , Siblings       * curr_sib //Siblings containing the node sh_v
      , Filtration_value fil_uv //Fil value of the edge uv in the zz filtration
      , int              k //Stopping condition for recursion based on max dim
      , std::vector<Simplex_handle> &zz_filtration) //range of all new simplices
    { //pick N^+(v)
  Simplex_handle root_sh_v = find_vertex(sh_v->first);
  if(!has_children(root_sh_v)) { return; }
  //intersect N^+(v) with labels y > v in curr_sib
  Simplex_handle next_it = sh_v;    ++next_it;
  thread_local std::vector< std::pair<Vertex_handle, Node> > inter;

  zz_intersection( inter
       , next_it
       , curr_sib->members().end()
       , root_sh_v->second.children()->members().begin()
       , root_sh_v->second.children()->members().end()
       , fil_uv );

  if(!inter.empty())
  { //the construction assign the self_siblings as children to all nodes
      Siblings * new_sib = new Siblings(curr_sib, sh_v->first, inter);
      sh_v->second.assign_children(new_sib);
      //update new Nodes and cofaces data structure
      for( auto new_sh = new_sib->members().begin();
     new_sh != new_sib->members().end(); ++new_sh )
      {
    // update_simplex_tree_after_node_insertion(new_sh);
    zz_filtration.push_back(new_sh);//new_sib does not change anymore
      }
      inter.clear();
      //recursive standard expansion for the rest of the subtree
      zz_siblings_expansion(new_sib, fil_uv, k-1, zz_filtration );
  }
  else { sh_v->second.assign_children(curr_sib); inter.clear(); }
    }


    //TODO boost::container::ordered_unique_range_t in the creation of a Siblings

    /* Global expansion of a subtree in the simplex tree.
   *
   * The filtration value is absolute and defined by "Filtration_value fil".
   * The new Node are also connected appropriately in the coface
   * data structure.
   */
    void zz_siblings_expansion(
      Siblings       * siblings  // must contain elements
      , Filtration_value fil
      , int              k  //==max_dim expansion - dimension curr siblings
      , std::vector<Simplex_handle> & zz_filtration )
    {
  if (k == 0) { return; } //max dimension
  Dictionary_it next = ++(siblings->members().begin());

  thread_local std::vector< std::pair<Vertex_handle, Node> > inter;
  for( Dictionary_it s_h = siblings->members().begin();
       next != siblings->members().end(); ++s_h, ++next)
  { //find N^+(s_h)
      Simplex_handle root_sh = find_vertex(s_h->first);
      if( has_children(root_sh) )
      {
    zz_intersection( inter                      // output intersection
         , next                       // begin
         , siblings->members().end()  // end
         , root_sh->second.children()->members().begin()
         , root_sh->second.children()->members().end()
         , fil   );

    if ( !inter.empty() )
    { //inter is of type boost::container::ordered_unique_range_t
        Siblings * new_sib = new Siblings( siblings    // oncles
                   , s_h->first  // parent
                   , inter);
        s_h->second.assign_children(new_sib);
        for( auto new_sh = new_sib->members().begin();
       new_sh != new_sib->members().end(); ++new_sh )
        {
      // update_simplex_tree_after_node_insertion(new_sh);//cofaces hooks
      zz_filtration.push_back(new_sh); //new_sib does not change anymore
        }
        inter.clear();
        //recursive standard expansion for the rest of the subtree
        zz_siblings_expansion(new_sib, fil, k - 1, zz_filtration);
    }     // ensure the children property
    else { s_h->second.assign_children(siblings); inter.clear();}
      }
  }
    }

    /* \brief Intersects Dictionary 1 [begin1;end1) with Dictionary 2 [begin2,end2)
   * and assigns Filtration_value fil to the Nodes.
   *
   * The function is identical to Simplex_tree::intersection(...) except that it
   * forces the filtration value fil for the new Nodes.
   *
   * todo merge zz_intersection and intersection with a
   * "filtration_strategy predicate".
   */
    static void zz_intersection(
      std::vector<std::pair<Vertex_handle, Node> > & intersection
      , Dictionary_it                                  begin1
      , Dictionary_it                                  end1
      , Dictionary_it                                  begin2
      , Dictionary_it                                  end2
      , Filtration_value                               fil )
    {
  if (begin1 == end1 || begin2 == end2) { return; }
  while (true) {
      if (begin1->first < begin2->first) {++begin1; if(begin1 == end1) {return;} }
      else
      {
    if (begin1->first > begin2->first) {++begin2; if(begin2 == end2) {return;} }
    else // begin1->first == begin2->first
    {
        intersection.emplace_back( begin1->first, Node( nullptr, fil ) );
        ++begin1; ++begin2;
        if (begin1 == end1 || begin2 == end2) { return; }
    }
      }
  }
    }

    //basic methods implemented for Nodes, and not Simplex_handle. The hooks in
    //cofaces_data_structure_ gives access to Nodes.
private:
    int dimension(Node & node, Vertex_handle u) {
  Siblings * curr_sib = self_siblings(node, u);
  int dim = 0;
  while (curr_sib != nullptr) {
      ++dim;
      curr_sib = curr_sib->oncles();
  }
  return dim - 1;
    }
    /* \brief Returns true if the node in the simplex tree pointed by
   * sh has children. node must have label u*/
    bool has_children(Node & node, Vertex_handle u) const {
  return (node.children()->parent() == u);
    }

public:
    /* Computes all simplices that ought to be removed
  * if the edge {u,v} were to disappear (puts them in zz_filtration).
  * This method does NOT modify the simplex tree.
  *
  * zz_filtration must be empty.
  */
    void flag_lazy_remove_edge(
      Vertex_handle u
      , Vertex_handle v
      , std::vector< Simplex_handle > & zz_filtration )
    {
// #ifdef COMPLEX_TIME
//   start_ = std::chrono::system_clock::now();
// #endif
  Simplex_handle sh_uv; //The simplex that get removed (edge or vertex)

  if(v < u) { std::swap(u,v); } //so as u <= v

  auto root_it_u = root_.members().find(u);
  if(root_it_u == root_.members().end()) {
// #ifdef COMPLEX_TIME
//       end_ = std::chrono::system_clock::now();
//       enlapsed_sec_ = end_ - start_;
//       complexOperationTime_ += enlapsed_sec_.count();
// #endif
      return;
  }//u not in Simplex_tree

  if( u == v ) {
      sh_uv = root_it_u;
      //keep track of all cofaces of the simplex removed, including simplex itself
      for(auto sh : star_simplex_range(sh_uv)) {zz_filtration.push_back(sh);}
// #ifdef COMPLEX_TIME
//       end_ = std::chrono::system_clock::now();
//       enlapsed_sec_ = end_ - start_;
//       complexOperationTime_ += enlapsed_sec_.count();
// #endif
      return;
  } //vertex u
  else { //edge {u,v}, u < v
      if(!(has_children(root_it_u))) { return; }// N^+(u) = \emptyset
      //Simplex_handle for edge {u,v}
      Simplex_handle sh_uv = root_it_u->second.children()->members().find(v);
      if(sh_uv == root_it_u->second.children()->members().end()) { return; }//edge not here
      //keep track of all cofaces of the simplex removed, including simplex itself
      for(auto sh : star_simplex_range(sh_uv)) {zz_filtration.push_back(sh);}
// #ifdef COMPLEX_TIME
//       end_ = std::chrono::system_clock::now();
//       enlapsed_sec_ = end_ - start_;
//       complexOperationTime_ += enlapsed_sec_.count();
// #endif
      return;
  }
    }
    /* Put all remaining simplices in the complex into zz_filtration, in order to
   * empty it. Does NOT modify the complex.
   */
    void flag_lazy_empty_complex(
      std::vector< Simplex_handle > & zz_filtration)
    {
// #ifdef COMPLEX_TIME
//   start_ = std::chrono::system_clock::now();
// #endif
  for(auto sh : complex_simplex_range()) { zz_filtration.push_back(sh); }
// #ifdef COMPLEX_TIME
//   end_ = std::chrono::system_clock::now();
//   enlapsed_sec_ = end_ - start_;
//   complexOperationTime_ += enlapsed_sec_.count();
// #endif
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
  /** Range for zigzag filtrations.*/
  Zigzagfiltration_simplex_range zigzag_simplex_range_;
  /** True iff the zigzag_simplex_range_ has been initialized */
  bool                           zigzag_simplex_range_initialized_;
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
  static const bool is_zigzag = false;
  typedef int Vertex_handle;
  typedef double Filtration_value;
  typedef std::uint32_t Simplex_key;
  static const bool store_key = true;
  static const bool store_filtration = true;
  static const bool contiguous_vertices = false;
  static const bool simplex_handle_strong_validity = false;
};

/** Model of SimplexTreeOptions, faster than `Simplex_tree_options_full_featured` but note the unsafe
 * `contiguous_vertices` option.
 * 
 * Maximum number of simplices to compute persistence is <CODE>std::numeric_limits<std::uint32_t>::max()</CODE>
 * (about 4 billions of simplices). */
struct Simplex_tree_options_fast_persistence {
  typedef linear_indexing_tag Indexing_tag;
  static const bool is_zigzag = false;
  typedef int Vertex_handle;
  typedef float Filtration_value;
  typedef std::uint32_t Simplex_key;
  static const bool store_key = true;
  static const bool store_filtration = true;
  static const bool contiguous_vertices = true;
  static const bool simplex_handle_strong_validity = false;
};

/** Model of SimplexTreeOptions.
 * 
 * Maximum number of simplices to compute persistence is <CODE>std::numeric_limits<std::uint32_t>::max()</CODE>
 * (about 4 billions of simplices). */
struct Simplex_tree_options_zigzag_persistence {
  typedef zigzag_indexing_tag Indexing_tag;
  static const bool is_zigzag = true;
  typedef int Vertex_handle;
  typedef double Filtration_value;
  typedef std::uint32_t Simplex_key;
  static const bool store_key = true;
  static const bool store_filtration = true;
  static const bool contiguous_vertices = false;
  //allows deletion without invalidating iterators (Simplex_handles)
  static const bool simplex_handle_strong_validity = true;
};

/** @} */  // end defgroup simplex_tree

}  // namespace Gudhi

#endif  // SIMPLEX_TREE_H_
