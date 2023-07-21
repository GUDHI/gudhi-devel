/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2021 Inria
 *
 *    Modification(s):
 *      - 2023/05 Hannah Schreiber: Rework of the interface, reorganization and debug
 *      - 2023/05 Hannah Schreiber: Addition of infinit bars
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file Zigzag_persistence.h
 * @author Clément Maria, Hannah Schreiber
 * @brief Contains the implementation of the @ref Gudhi::zigzag_persistence::Zigzag_persistence class.
 */

#ifndef ZIGZAG_PERSISTENCE_H_
#define ZIGZAG_PERSISTENCE_H_

#include <boost/tuple/tuple.hpp>
#include <boost/intrusive/list.hpp>
#include <boost/intrusive/set.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/pool/object_pool.hpp>
#include <boost/timer/progress_display.hpp>

#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <list>
#include <ostream>
#include <set>
#include <unordered_set>
#include <utility>
#include <vector>
#include <functional>

#include <gudhi/Debug_utils.h>
#include <gudhi/matrix.h>
#include <gudhi/options.h>
#include <gudhi/utilities/utilities.h>

namespace Gudhi {
namespace zigzag_persistence {

/** \class Zigzag_persistence Zigzag_persistence.h gudhi/Zigzag_persistence.h
 * \brief Class computating the zigzag persistent homology of a zigzag
 * filtration. Algorithm based on \cite zigzag_reflection.
 *
 * \details The type ZigzagFilteredComplex::Simplex_key counts the number of
 * insertions and
 * deletions of simplices, which may be large in zigzag persistence and require
 * more than 32 bits of storage. The type used (int, long, etc) should be chosen in
 * consequence. Simplex_key must be signed.
 *
 * Over all insertions, the Simplex_key must be positive and strictly increasing
 * when forward iterating along the zigzag filtration.
 *
 * \tparam ZigzagFilteredComplex Complex storing the current simplices.
 * \tparam ZigzagPersistenceOptions Options for the matrix used to compute the persistence.
 */
template <typename ZigzagFilteredComplex,
          typename ZigzagPersistenceOptions = Gudhi::persistence_matrix::Zigzag_options<> >
class Zigzag_persistence 
{
 public:
  using Complex = ZigzagFilteredComplex;						/**< Complex type. */
  using Options = ZigzagPersistenceOptions;						/**< Matrix options */
  /*** Types defined in the complex ***/
  using Simplex_key = typename Complex::Simplex_key;			/**< Key type, must be signed. */
  using Simplex_handle = typename Complex::Simplex_handle;		/**< Simplex ID type in the complex. */
  using Vertex_handle = typename Complex::Vertex_handle;		/**< Vertex ID type in the complex. */
  using Filtration_value = typename Complex::Filtration_value;	/**< Filtration value type. */

  /** \brief Structure to store persistence intervals by their index values.
   *
   * \details By convention, interval [b;d] are
   * closed for finite indices b and d, and open for left-infinite and/or
   * right-infinite endpoints.
   */
  template<typename value_type>
  struct interval {
    interval() {}
    interval(int dim, value_type b, value_type d) : dim_(dim), b_(b), d_(d) {}
    /** Returns the dimension of the homological feature corresponding to the
     * interval. */
    int dim() const { return dim_; }  // return the homological dimension of the interval
    /** Returns the birth index of the interval.*/
    value_type birth() const { return b_; }  // return the birth value
    /** Returns the death index of the interval.*/
    value_type death() const { return d_; }  // return the death value

   protected:          // note that we don't assume b_ <= d_
    int dim_;        // homological dimension
    value_type b_;  // filtration value associated to birth index
    value_type d_;  // filtration value associated to death index
  };
  using index_interval = interval<Simplex_key>;

  /** \brief Structure to store persistence intervals by their filtration values.
   *
   * \details By convention, interval \f$[b;d]\f$ are
   * closed for finite indices b and d, and open for left-infinite and/or
   * right-infinite endpoints.
   */
  struct filtration_value_interval : interval<Filtration_value> 
  {
   private:
    using Base = interval<Filtration_value>;

   public:
    /**
     * @brief Default constructor
     */
    filtration_value_interval() : Base() {}
    /**
     * @brief Construct a new interval with given parameters
     *
     * @param dim Dimension of the interval.
     * @param b Start value of the interval.
     * @param d End value of the interval.
     */
    filtration_value_interval(int dim, Filtration_value b, Filtration_value d)
        : Base(dim, b, d) {}

    /**
     * @brief Returns the absolute length of the interval \f$|d-b|\f$.
     */
    Filtration_value length() const {
      if (Base::b_ == Base::d_) {
        return 0;
      }  // otherwise inf - inf would return nan.
      return std::abs(Base::b_ - Base::d_);
    }
    /**
     * @brief Returns the absolute length of the log values of birth and death, i.e.  \f$|\log d - \log b|\f$.
     */
    Filtration_value log_length() const {
      if (Base::b_ == Base::d_) {
        return 0;
      }  // otherwise inf - inf would return nan.
      return std::abs(std::log2(static_cast<double>(Base::b_)) - std::log2(static_cast<double>(Base::d_)));
    }
  };

 private:
  /** \brief Maintains the birth ordering \f$\leq_b\f$.
   *
   * \details Contains an std::map of size the number of
   * non-zero rows of the homology matrix, at any time during the computation of
   * zigzag persistence.
   *
   * By construction, we maintain the map satisfying
   * 'birth_to_pos_[i] < birth_to_pos_[j]', with \f$0 <= i,j <= k\f$ indices in the quiver 
   * '\f$0 \leftrightarrow ... \leftrightarrow i \leftrightarrow ... \leftrightarrow k\f$' visited at time 
   * \f$k\f$ of the algorithm (prefix of length \f$k\f$ of the full zigzag filtration 
   * '\f$0 \leftrightarrow ... \leftrightarrow i \leftrightarrow ... 
   * \leftrightarrow k \leftrightarrow ... \leftrightarrow n\f$'
   * that is studied), iff \f$i <_b j\f$ for the birth ordering.
   *
   * By construction, when adding index \f$k+1\f$ to '\f$0 \leftrightarrow ... \leftrightarrow i \leftrightarrow ..
   * \leftrightarrow k \leftrightarrow k+1'\f$, we have:
   * - if \f$k \rightarrow k+1\f$ forward, then \f$j <_b k+1\f$ for all indices \f$j < k+1\f$, otherwise
   * - if \f$k \leftarrow k+1\f$ backward, then \f$k+1 <_b j\f$ for all indices \f$j < k+1\f$.
   */
  struct birth_ordering {
    /**
     * @brief Default constructor
     */
    birth_ordering() : birth_to_pos_(), max_birth_pos_(0), min_birth_pos_(-1) {}

    /**
     * @brief Inserts arrow number in the ordering after an insertion.
     * When the arrow key-1 -> key is forward, key is larger than any other index
     * i < key in the birth ordering <b. We give key the largest value max_birth_pos_.
     *
     * @param arrow_number Forward arrow number.
     */
    void add_birth_forward(Simplex_key arrow_number) {  // amortized constant time
      birth_to_pos_.emplace_hint(birth_to_pos_.end(), arrow_number, max_birth_pos_);
      ++max_birth_pos_;
    }
    /**
     * @brief Inserts arrow number in the ordering after a removal.
     * When the arrow key-1 <- key is backward, key is smaller than any other index
     * i < key in the birth ordering <b. We give key the smallest value min_birth_pos_
     *
     * @param arrow_number Backward arrow number.
     */
    void add_birth_backward(Simplex_key arrow_number) {  // amortized constant time
      birth_to_pos_.emplace_hint(birth_to_pos_.end(), arrow_number, min_birth_pos_);
      --min_birth_pos_;
    }

    /**
     * @brief Removes the birth from the ordering.
     * When the row at index @a birth is removed from the homology matrix, we do not need
     * to maintain its position in <b anymore.
     *
     * @param birth Birth to remove.
     */
    void remove_birth(Simplex_key birth) { birth_to_pos_.erase(birth); }
    /**
     * @brief Increasing birth order <=b, true iff k1 <b k2.
     *
     * @param k1
     * @param k2
     * @return true if k1 <b k2, false otherwise.
     */
    bool birth_order(Simplex_key k1, Simplex_key k2) const { return birth_to_pos_.at(k1) < birth_to_pos_.at(k2); }
    /**
     * @brief Decreasing birth order <=b, true iff k1 >b k2.
     *
     * @param k1
     * @param k2
     * @return true if k1 >b k2, false otherwise.
     */
    bool reverse_birth_order(Simplex_key k1, Simplex_key k2) const { return birth_to_pos_.at(k1) > birth_to_pos_.at(k2); }

   private:
    std::unordered_map<Simplex_key, Simplex_key> birth_to_pos_;	/**< birth_to_pos_[i] < birth_to_pos_[j] iff i <b j */
    Simplex_key max_birth_pos_;							/**< is strictly larger than any other birth so far */
    Simplex_key min_birth_pos_;							/**< is strictly smaller than any other birth so far */
  };

  using matrix_type = Gudhi::persistence_matrix::Matrix<ZigzagPersistenceOptions>;
  using index = Gudhi::persistence_matrix::index;

 public:
  /**
   * @brief Constructor of the Zigzag_persistence class.
   * @details After construction of the class, the zigzag filtration should be given in a streaming like way, i.e., 
   * call @ref insert_simplex or @ref remove_simplex for each step of the filtration in order of the filtration. 
   * If simplices are added (resp. removed) continuously, they can be inserted (resp. removed) in batches by using 
   * @ref insert_simplices_contiguously (resp. @ref remove_simplices_contiguously).
   * To retrieve the current persistence diagram at any moment of the filtration, 
   * use @ref get_persistence_diagram or @ref get_index_persistence_diagram.
   *
   * @param min_number_of_simplices Minimum number of simplices that will be inserted at some point in the filtration.
   * If the total number of simplices is known in advance, the memory allocation can be better optimized. 
   * Default value: 0.
   * @param ignore_cycles_above_dim Ignores cycles in dimension larger or equal in the final diagram. 
   * If -1, no cycles are ignored. Default value: -1.
   */
  Zigzag_persistence(unsigned int min_number_of_simplices = 0, int ignore_cycles_above_dim = -1)
      : dim_max_(ignore_cycles_above_dim),
        matrix_(min_number_of_simplices,
                [this](index columnIndex1, index columnIndex2) {
                  return birth_ordering_.birth_order(births_.at(columnIndex1), births_.at(columnIndex2));
                }),
        num_arrow_(-1),
        previous_filtration_value_(std::numeric_limits<Filtration_value>::infinity()) {}

  /**
   * @brief Updates the zigzag persistence diagram after the insertion of the given simplex.
   *
   * @tparam VertexRange Range type needing begin and end members.
   * @param simplex Simplex to insert, represented by its vertices.
   * @param filtration_value Filtration value associated to the simplex. 
   * Assumed to be larger or equal to previously used filtration values.
   */
  template <class VertexRange = std::initializer_list<Vertex_handle>>
  void insert_simplex(const VertexRange& simplex, Filtration_value filtration_value) {
    if (dim_max_ != -1 && simplex.size() > static_cast<unsigned int>(dim_max_) + 1) return;

    ++num_arrow_;

    if (filtration_value != previous_filtration_value_)  // check whether the filt value has changed
    {  // consecutive pairs (i,f), (j,f') mean simplices of index k in [i,j-1] have
      previous_filtration_value_ = filtration_value;  // filtration value f
      filtration_values_.emplace_back(num_arrow_, previous_filtration_value_);
    }

    std::pair<Simplex_handle, bool> res = cpx_.insert_simplex(simplex, filtration_value);
    GUDHI_CHECK(res.second, "Zigzag_persistence::insert_simplex - insertion of a simplex already in the complex");
    cpx_.assign_key(res.first, num_arrow_);
    _process_forward_arrow(res.first);
  }

  /**
   * @brief Updates the zigzag persistence diagram after the removal of the given simplex.
   *
   * @tparam VertexRange Range type needing begin and end members.
   * @param simplex Simplex to remove, represented by its vertices.
   * @param filtration_value Filtration value associated to the removal. 
   * Assumed to be larger or equal to previously used filtration values.
   */
  template <class VertexRange = std::initializer_list<Vertex_handle>>
  void remove_simplex(const VertexRange& simplex, Filtration_value filtration_value) {
    if (dim_max_ != -1 && simplex.size() > static_cast<unsigned int>(dim_max_) + 1) return;

    ++num_arrow_;

    Simplex_handle sh = cpx_.find(simplex);
    GUDHI_CHECK(sh != cpx_.null_simplex(),
                "Zigzag_persistence::remove_simplex - removal of a simplex not in the complex");

    if (filtration_value != previous_filtration_value_)  // check whether the filt value has changed
    {  // consecutive pairs (i,f), (j,f') mean simplices of index k in [i,j-1] have
      previous_filtration_value_ = filtration_value;  // filtration value f
      filtration_values_.emplace_back(num_arrow_, previous_filtration_value_);
    }

    _process_backward_arrow(sh);
    cpx_.remove_maximal_simplex(sh);
  }

  /**
   * @brief Updates the zigzag persistence diagram after the insertion of the given simplices.
   *
   * @tparam SimplexRange Range type needing begin and end members.
   * @tparam FiltrationRange Range type needing begin and end members.
   * @param simplices Simplices which are inserted, represented by their vertices. They have to be in the order they 
   * are inserted in the filtration and ``contiguous'' in the filtration, that is, no other simplex 
   * which is not in the range is inserted or removed between two simplices in the range.
   * @param filtration_values Filtration values associated to the insertion of the given simplices. 
   * The order has to correspond to the order in @a simplices. Their values have to ascending in this order and 
   * they are assumed to be larger or equal to previously used filtration values.
   */
  template <class SimplexRange = std::initializer_list<std::initializer_list<Vertex_handle>>,
            class FiltrationRange = std::initializer_list<Filtration_value>>
  void insert_simplices_contiguously(const SimplexRange& simplices, const FiltrationRange& filtration_values) {
    auto simplexIt = simplices.begin();
    auto filIt = filtration_values.begin();
    for (; simplexIt != simplices.end(); ++simplexIt, ++filIt) {
      insert_simplex(*simplexIt, *filIt);
    }
  }

  /**
   * @brief Updates the zigzag persistence diagram after the removal of the given simplices.
   *
   * @tparam SimplexRange Range type needing begin and end members.
   * @tparam FiltrationRange Range type needing begin and end members.
   * @param simplices Simplices which are removed, represented by their vertices. They have to be in the order they
   * are removed in the filtration and ``contiguous'' in the filtration, that is, no other simplex
   * which is not in the range is inserted or removed between two simplices in the range.
   * @param filtration_values Filtration values associated to the removal of the given simplices. Has therefore the 
   * same size as @a simplices. The order has to correspond to the order in @a simplices. Their values have to 
   * ascending in this order and they are assumed to be larger or equal to previously used filtration values.
   */
  template <class SimplexRange = std::initializer_list<std::initializer_list<Vertex_handle>>,
            class FiltrationRange = std::initializer_list<Filtration_value>>
  void remove_simplices_contiguously(const SimplexRange& simplices, const FiltrationRange& filtration_values) {
    auto simplexIt = simplices.begin();
    auto filIt = filtration_values.begin();
    for (; simplexIt != simplices.end(); ++simplexIt, ++filIt) {
      remove_simplex(*simplexIt, *filIt);
    }
  }

  /**
   * @brief Updates the zigzag persistence diagram after the insertion of the given simplices.
   *
   * @tparam SimplexRangeIterators Forward iterator of a range.
   * @tparam FiltrationRangeIterators Forward iterator of a range.
   * @param simplex_range_start Iterator pointing to the begining of the range of simplices to insert.
   * The simplices should be represented by their vertices. They have to be in the order they 
   * are inserted in the filtration and ``contiguous'' in the filtration, that is, no other simplex 
   * which is not in the range is inserted or removed between two simplices in the range.
   * @param simplex_range_end Iterator pointing to the end of the range of simplices to insert.
   * @param filtration_range_start Iterator pointing to the begining of the range of filtration values. The range is
   * assumed to end at the same time than the simplices range and has the same order. The filtration values should be 
   * ascending in this order and they are assumed to be larger or equal to previously used filtration values.
   */
  template <class SimplexRangeIterators, class FiltrationRangeIterators>
  void insert_simplices_contiguously(SimplexRangeIterators simplex_range_start, SimplexRangeIterators simplex_range_end,
                                     FiltrationRangeIterators filtration_range_start) {
    for (; simplex_range_start != simplex_range_end; ++simplex_range_start, ++filtration_range_start) {
      insert_simplex(*simplex_range_start, *filtration_range_start);
    }
  }

  /**
   * @brief Updates the zigzag persistence diagram after the removal of the given simplices.
   *
   * @tparam SimplexRangeIterators Forward iterator of a range.
   * @tparam FiltrationRangeIterators Forward iterator of a range.
   * @param simplex_range_start Iterator pointing to the begining of the range of simplices to remove.
   * The simplices should be represented by their vertices. They have to be in the order they 
   * are removed in the filtration and ``contiguous'' in the filtration, that is, no other simplex 
   * which is not in the range is inserted or removed between two simplices in the range.
   * @param simplex_range_end Iterator pointing to the end of the range of simplices to remove.
   * @param filtration_range_start Iterator pointing to the begining of the range of filtration values. The range is
   * assumed to end at the same time than the simplices range and has the same order. The filtration values should be 
   * ascending in this order and they are assumed to be larger or equal to previously used filtration values.
   */
  template <class SimplexRangeIterators, class FiltrationRangeIterators>
  void remove_simplices_contiguously(SimplexRangeIterators simplex_range_start, SimplexRangeIterators simplex_range_end,
                                     FiltrationRangeIterators filtration_range_start) {
    for (; simplex_range_start != simplex_range_end; ++simplex_range_start, ++filtration_range_start) {
      remove_simplex(*simplex_range_start, *filtration_range_start);
    }
  }

  /**
   * @brief Returns the ``index persistence diagram'' of the current filtration, that is, the pairs of atomic arrow 
   * numbers corresponding to a birth-death pair. Does not contain points at infinity, only the cycle classes which 
   * already died are represented.
   *
   * @return Reference to the list of intervals.
   */
  const std::list<index_interval>& get_index_persistence_diagram() const { return persistence_diagram_; }

  /**
   * @brief Returns the filtration values \f$[f(b),f(d)]\f$ associated to the indices \f$[b,d]\f$ which are retrieved
   * by @ref get_index_persistence_diagram.
   *
   * @param b_key Birth index
   * @param d_key Death index
   * @return A pair of filtration values associated to the given indices.
   */
  std::pair<Filtration_value, Filtration_value> map_index_to_filtration_value(
    Simplex_key b_key, Simplex_key d_key) const 
  {
    // filtration_values_ must be sorted by increasing keys.
    auto it_b =  // lower_bound(x) returns leftmost y s.t. x <= y
        std::lower_bound(
            filtration_values_.begin(), filtration_values_.end(),
            std::pair<Simplex_key, Filtration_value>(b_key, std::numeric_limits<Filtration_value>::infinity()),
            [](std::pair<Simplex_key, Filtration_value> p1, std::pair<Simplex_key, Filtration_value> p2) {
              return p1.first < p2.first;
            });
    if (it_b == filtration_values_.end() || it_b->first > b_key) {
      --it_b;
    }
    // it points to the rightmost z such that z <= x

    auto it_d =  //
        std::lower_bound(
            filtration_values_.begin(), filtration_values_.end(),
            std::pair<Simplex_key, Filtration_value>(d_key, std::numeric_limits<Filtration_value>::infinity()),
            [](std::pair<Simplex_key, Filtration_value> p1, std::pair<Simplex_key, Filtration_value> p2) {
              return p1.first < p2.first;
            });
    if (it_d == filtration_values_.end() || it_d->first > d_key) {
      --it_d;
    }

    return std::make_pair(it_b->second, it_d->second);
  }

  /**
   * @brief Returns the current persistence diagram ordered first by length, than by dimension, 
   * than by birth value and finally by death value.
   * 
   * @param shortest_interval Threshold. Every bar shorter than the given value will be ignored. Default value: 0.
   * @param include_infinit_bars If set to true, infinit bars are included in the diagram. Default value: false.
   * @return A vector of pairs of filtration values representing the persistence diagram.
   */
  std::vector<filtration_value_interval> get_persistence_diagram(Filtration_value shortest_interval = 0.,
                                                                 bool include_infinit_bars = false) {
    auto comp = [](filtration_value_interval p, filtration_value_interval q) {
      if (p.length() != q.length()) {
        return p.length() > q.length();
      }  // longest 1st
      if (p.dim() != q.dim()) {
        return p.dim() < q.dim();
      }  // lower dimension first
      if (p.birth() != q.birth()) {
        return p.birth() < q.birth();
      }  // lex order
      return p.death() < q.death();
    };

    std::vector<filtration_value_interval> diag = _get_persistence_diagram(shortest_interval);

    if (include_infinit_bars) {
      _retrieve_infinit_bars(diag);
    }

    std::stable_sort(diag.begin(), diag.end(), comp);

    return diag;
  }

  /**
   * @brief Returns a reference to the complex storing the simplices. 
   * A simplex is added in a call of @ref insert_simplex and is removed in a call of @ref remove_simplex.
   *
   * @return Const reference to the complex.
   */
//   const ZigzagFilteredComplex& get_complex() const{
//     return cpx_;
//   }

  /**
   * @brief Returns a reference to the complex storing the simplices. 
   * A simplex is added in a call of @ref insert_simplex and is removed in a call of @ref remove_simplex.
   * @warning The complex is not const for now for technical reasons, but DO NOT modify it.
   *
   * @return Reference to the complex.
   */
  ZigzagFilteredComplex& get_complex() const{
    return cpx_;
  }

  /**
   * @brief For debug purposes, to remove.
   */
  void print_current_complex() const {
    for (auto& sh : cpx_.complex_simplex_range()) {
      for (auto v : cpx_.simplex_vertex_range(sh)) {
        std::cout << v << " ";
      }
      std::cout << " - " << cpx_.filtration(sh) << "\n";
    }
  }

 private:
  /**
   * @brief Computes the boundary cycle of the new simplex zzsh, and express it as a
   * sum of cycles in a matrix. If some cycles are not boundary cycles, i.e., columns with F-index
   * in the matrix, it applies a surjective diamond to the zigzag module.
   *
   * @param zzsh Simplex handle of the inserted simplex.
   */
  void _process_forward_arrow(Simplex_handle zzsh) {  // maintain the <=b order
    // Reduce the boundary of zzsh in the basis of cycles.
    // Compute the simplex keys of the simplices of the boundary of zzsh.
    std::set<Simplex_key> col_bsh;  // set maintains the natural order on indices
    for (auto b_sh : cpx_.boundary_simplex_range(zzsh)) {
      col_bsh.insert(cpx_.key(b_sh));
    }

    std::vector<index> chains_in_F = matrix_.insert_boundary(num_arrow_, col_bsh);

    if (!chains_in_F.empty()) {
      _apply_surjective_reflection_diamond(zzsh, chains_in_F);
    } else {
      birth_ordering_.add_birth_forward(num_arrow_);
      births_.emplace_hint(births_.end(), matrix_.get_column_with_pivot(num_arrow_), num_arrow_);
    }
  }

  /**
   * @brief Applies the surjective reflection diamond principle to the current filtration.
   * 
   * @details The vector chains_in_F is sorted by decreasing lowest index values in the
   * columns corresponding to the chains, due to its computation in the reduction of
   * \partial zzsh in _process_forward_arrow(...). It is equivalent to decreasing death index
   * order w.r.t. the <d ordering.
   *
   * @param zzsh Simplex handle of the inserted simplex.
   * @param chains_in_F Indices of the non paired columns in the matrix.
   */
  void _apply_surjective_reflection_diamond(Simplex_handle zzsh,
                                            const std::vector<index>& chains_in_F) { 
    // fp is the largest death index for <=d
    // Set col_fp: col_fp <- col_f1+...+col_fp (now in G); preserves lowest idx
    auto chain_fp = chains_in_F[0];  // col_fp, with largest death <d index.

    // chains_in_F is ordered, from .begin() to end(), by decreasing lowest_idx_. The
    // lowest_idx_ is also the death of the chain in the right suffix of the
    // filtration (all backward arrows). Consequently, the chains in F are ordered by
    // decreasing death for <d.
    // Pair the col_fi, i = 1 ... p-1, according to the reflection diamond principle
    // Order the fi by reverse birth ordering <=_b
    auto cmp_birth = [this](Simplex_key k1, Simplex_key k2) -> bool {
      return birth_ordering_.reverse_birth_order(k1, k2);
    };  // true iff b(k1) >b b(k2)

    // available_birth: for all i by >d value of the d_i,
    // contains at step i all b_j, j > i, and maybe b_i if not stolen
    std::set<Simplex_key, decltype(cmp_birth)> available_birth(cmp_birth);
    // for f1 to f_{p} (i by <=d), insertion in available_birth_to_fidx sorts by >=b
    for (auto& chain_f : chains_in_F) {
      available_birth.insert(births_.at(chain_f));
    }

    auto maxb_it = available_birth.begin();  // max birth cycle
    auto maxb = *maxb_it;                    // max birth value, for persistence diagram
    available_birth.erase(maxb_it);          // remove max birth cycle (stolen)

    auto last_modified_chain_it = chains_in_F.rbegin();

    // consider all death indices by increasing <d order i.e., increasing lowest_idx_
    for (auto chain_f_it = chains_in_F.rbegin();  // by increasing death order <d
         *chain_f_it != chain_fp; ++chain_f_it)   // chain_fp=*begin() has max death
    {                                             // find which reduced col has this birth
      auto birth_it = available_birth.find(births_.at(*chain_f_it));
      if (birth_it == available_birth.end())  // birth is not available. *chain_f_it
      {                                       // must become the sum of all chains in F with smaller death index.
        // this gives as birth the maximal birth of all chains with strictly larger
        // death <=> the maximal availabe death.
        // Let c_1 ... c_f be the chains s.t. <[c_1+...+c_f]> is the kernel and
        //  death(c_i) >d death(c_i-1). If the birth of c_i is not available, we set
        // c_i <- c_i + c_i-1 + ... + c_1, which is [c_i + c_i-1 + ... + c_1] on
        // the right (of death the maximal<d death(c_i)), and is [c_i + c_i-1 + ... +
        // c_1] + kernel = [c_f + c_f-1 + ... + c_i+1] on the left (of birth the max<b
        // of the birth of the c_j, j>i  <=> the max<b available birth).
        // N.B. some of the c_k, k<i, ahve already been modified to be equal to
        // c_k + c_k-1 + ... + c_1. The largest k with this property is maintained in
        // last_modified_chain_it (no need to compute from scratch the full sum).

        // last_modified is equal to c_k+...+c_1, all c_j, i>j>k, are indeed c_j
        // set c_i <- c_i + (c_i-1) + ... + (c_k+1) + (c_k + ... + c_1)
        for (auto chain_passed_it = last_modified_chain_it;  // all with smaller <d death
             chain_passed_it != chain_f_it; ++chain_passed_it) {
          matrix_.add_to(*chain_passed_it, *chain_f_it);
        }
        last_modified_chain_it = chain_f_it;  // new cumulated c_i+...+c_1
        // remove the max available death
        auto max_avail_b_it = available_birth.begin();  // max because order by deacr <b
        Simplex_key max_avail_b = *max_avail_b_it;      // max available birth

        births_.at(*chain_f_it) = max_avail_b;  // give new birth
        available_birth.erase(max_avail_b_it);  // remove birth from availability
      } else {
        available_birth.erase(birth_it);
      }  // birth not available anymore, do not
    }    // modify *chain_f_it.

    //Following value can be erased, but it slowes the process down a bit, so I keep it as a remainder for now:
    // birth_ordering_.remove_birth(maxb);
    births_.erase(chain_fp);

    // Update persistence diagram with left interval [fil(b_max) ; fil(m))
    persistence_diagram_.emplace_back(cpx_.dimension(zzsh) - 1, maxb, num_arrow_);  //-1);//
  }

  /**
   * @brief Removes the given simplex by pushing up the matrix the corresponding column and erasing it.
   * 
   * @param zzsh Simplex handle of the simplex to remove.
   */
  void _process_backward_arrow(Simplex_handle zzsh) {
    Simplex_key simplexIndex = cpx_.key(zzsh); // cpx_.key(zzsh) is the key of the simplex we remove, not a new one
    // column whose key is the one of the removed simplex
    index curr_col = matrix_.get_column_with_pivot(simplexIndex);
    // Record all columns that get affected by the transpositions, i.e., have a coeff
    std::vector<index> modified_columns;
    const auto& row = matrix_.get_row(simplexIndex);
    modified_columns.reserve(row.size());
    std::transform(row.begin(), row.end(), std::back_inserter(modified_columns),
                   [](const auto& cell) { return cell.get_column_index(); });
    // Sort by left-to-right order in the matrix_ (no order maintained in rows)
    std::stable_sort(modified_columns.begin(), modified_columns.end(),
                     [this](index i1, index i2) { return matrix_.get_pivot(i1) < matrix_.get_pivot(i2); });

    // Modifies curr_col, not the other one.
    for (auto other_col_it = std::next(modified_columns.begin()); other_col_it != modified_columns.end();
         ++other_col_it) {
      curr_col = matrix_.vine_swap_with_z_eq_1_case(curr_col, *other_col_it);
    }

    // curr_col points to the column to remove by restriction of K to K-{\sigma}
    if (!matrix_.get_column(curr_col).is_paired()) {  // in F
      int dim_zzsh = cpx_.dimension(zzsh);
      auto it = births_.find(curr_col);
      if (dim_max_ == -1 || (dim_max_ != -1 && dim_zzsh < dim_max_)) {  // don't record intervals of max dim
        persistence_diagram_.emplace_back(dim_zzsh, it->second, num_arrow_);  // -1);
      }
      //Following value can be erased, but it slowes the process down a bit, so I keep it as a remainder for now:
      // birth_ordering_.remove_birth(it->second);
      births_.erase(it);
    } else {  // in H    -> paired with c_g, that now belongs to F now
      // maintain the <=b order
      birth_ordering_.add_birth_backward(num_arrow_);
      births_.try_emplace(matrix_.get_column(curr_col).get_paired_chain_index(), num_arrow_);
    }

    // cannot be in G as the removed simplex is maximal
    matrix_.remove_maximal_simplex(simplexIndex);
  }

  /**
   * @brief Returns the current persistence diagram ordered by length without infinit bars.
   *
   * @param shortest_interval Intervals shorter than the given value are ignored.
   * @return Vector of intervals.
   */
  std::vector<filtration_value_interval> _get_persistence_diagram(Filtration_value shortest_interval) {
    std::vector<filtration_value_interval> diag;
    diag.reserve(persistence_diagram_.size());

    std::stable_sort(filtration_values_.begin(), filtration_values_.end(),
                     [](std::pair<Simplex_key, Filtration_value> p1, std::pair<Simplex_key, Filtration_value> p2) {
                       return p1.first < p2.first;
                     });

    for (auto bar : persistence_diagram_) {
      Filtration_value birth, death;
      std::tie(birth, death) = map_index_to_filtration_value(bar.birth(), bar.death());
      if (birth > death) {
        std::swap(birth, death);
      }

      if (death - birth > shortest_interval) {
        diag.emplace_back(bar.dim(), birth, death);
      }
    }

    return diag;
  }

  /**
   * @brief Computes the births of the current essential cycles.
   *
   * @param diag Reference to vector where to store the infinit bars.
   */
  void _retrieve_infinit_bars(std::vector<filtration_value_interval>& diag) const {
    auto birth = [this](Simplex_key b_key) {
      auto it_b =  // lower_bound(x) returns leftmost y s.t. x <= y
          std::lower_bound(
              filtration_values_.begin(), filtration_values_.end(),
              std::pair<Simplex_key, Filtration_value>(b_key, std::numeric_limits<Filtration_value>::infinity()),
              [](std::pair<Simplex_key, Filtration_value> p1, std::pair<Simplex_key, Filtration_value> p2) {
                return p1.first < p2.first;
              });
      if (it_b == filtration_values_.end() || it_b->first > b_key) {
        --it_b;
      }
      return it_b->second;
    };

    for (auto& p : births_) {
      auto dim = matrix_.get_column(matrix_.get_column_with_pivot(p.first)).get_dimension();
      if (dim_max_ == -1 || (dim_max_ != -1 && dim < dim_max_)) 
        diag.emplace_back(dim, birth(p.second), std::numeric_limits<Filtration_value>::infinity());
    }
  }

 private:
  Complex cpx_;                                     /**< Complex in which the current simplices are stored. */
  int dim_max_;                                     /**< Maximal dimension of a bar to record. */
  matrix_type matrix_;                              /**< Matrix storing a base of the current chain complex. */
  std::unordered_map<index, int> births_;           /**< Map simplex index in F to corresponding birth. */
  birth_ordering birth_ordering_;                   /**< Maintains <b ordering of the births. */
  std::list<index_interval> persistence_diagram_;   /**< Stores current closed persistence intervals. */
  Simplex_key num_arrow_;                           /**< Current arrow number. */
  Filtration_value previous_filtration_value_;      /**< Filtration value of the previous arrow. */
  /**
   * @brief filtration_values_ stores consecutive pairs (i,f) , (j,f') with f != f',
   * meaning that all inserted simplices with key in [i;j-1] have filtration value f
   * i is the smallest simplex index whose simplex has filtration value f.
   */
  std::vector<std::pair<Simplex_key, Filtration_value>> filtration_values_;
};  // end class Zigzag_persistence

}  // namespace zigzag_persistence

}  // namespace Gudhi

#endif  // ZIGZAG_PERSISTENCE_H_
