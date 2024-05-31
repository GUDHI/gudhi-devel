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

#include <cmath>
#include <limits>
#include <list>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

#include <gudhi/Debug_utils.h>
#include <gudhi/matrix.h>

namespace Gudhi {
namespace zigzag_persistence {

/**
 * @brief Options for the internal matrix of @ref Zigzag_persistence.
 * 
 * @tparam column_type Column type of the matrix.
 */
template <Gudhi::persistence_matrix::Column_types column_type>
struct Zigzag_matrix_options : Gudhi::persistence_matrix::Default_options<column_type, true> {
  static const bool has_row_access = true;
  static const bool has_column_pairings = false;
  static const bool has_vine_update = true;
  static const bool is_of_boundary_type = false;
  static const bool has_map_column_container = true;
  static const bool has_removable_columns = true;
  static const bool has_removable_rows = true;
};

/**
 * @brief Default options for @ref Zigzag_persistence.
 */
struct Default_zigzag_options {
  using internal_key = int;         /**< Face ID used internaly, must be signed. */
  using face_key = int;             /**< Face ID used in the given boundaries. */
  using filtration_value = double;  /**< Filtration value type. */
  using dimension_type = int;       /**< Dimension value type. */
  /**
   * @brief Column type use by the internal matrix.
   */
  static const Gudhi::persistence_matrix::Column_types column_type =
      Gudhi::persistence_matrix::Column_types::INTRUSIVE_LIST;
};

//TODO: add the possibility of something else than Z2. Which means that the possibility of vineyards without Z2
//also needs to be implemented. The theory needs to be done first.
/** \class Zigzag_persistence Zigzag_persistence.h gudhi/Zigzag_persistence.h
 * \brief Class computating the zigzag persistent homology of a zigzag
 * filtration. Algorithm based on \cite zigzag.
 *
 * \ingroup zigzag_persistence
 *
 * \tparam ZigzagOptions Structure following the @ref ZigzagOptions concept. Default value: @ref Default_zigzag_options.
 */
template <class ZigzagOptions = Default_zigzag_options>
class Zigzag_persistence {
 public:
  using Options = ZigzagOptions;                                      /**< Zigzag options. */
  using Matrix_options = Zigzag_matrix_options<Options::column_type>; /**< Matrix options. */
  using internal_key = typename Options::internal_key;                /**< Key and index type, has to be signed. */
  using face_key = typename Options::face_key;                        /**< Face ID type from external inputs. */
  using filtration_value = typename Options::filtration_value;        /**< Type for filtration values. */
  using dimension_type = typename Options::dimension_type;            /**< Type for dimension values. */

  /** \brief Structure to store persistence intervals by their index values.
   *
   * \details By convention, interval [b;d] are
   * closed for finite indices b and d, and open for left-infinite and/or
   * right-infinite endpoints.
   */
  template <typename value_type>
  struct Interval {
    Interval() {}
    Interval(int dim, value_type b, value_type d) : dim_(dim), b_(b), d_(d) {}
    /** Returns the dimension of the homological feature corresponding to the interval. */
    int dim() const { return dim_; }
    /** Returns the birth index of the interval.*/
    value_type birth() const { return b_; }
    /** Returns the death index of the interval.*/
    value_type death() const { return d_; }

   protected:
    int dim_;       // homological dimension
    value_type b_;  // filtration value associated to birth index
    value_type d_;  // filtration value associated to death index
  };
  using Index_interval = Interval<internal_key>;

  /** \brief Structure to store persistence intervals by their filtration values.
   *
   * \details By convention, interval \f$[b;d]\f$ are
   * closed for finite indices b and d, and open for left-infinite and/or
   * right-infinite endpoints.
   */
  struct Filtration_value_interval : Interval<filtration_value> {
   private:
    using Base = Interval<filtration_value>;

   public:
    /**
     * @brief Default constructor
     */
    Filtration_value_interval() : Base() {}
    /**
     * @brief Construct a new interval with given parameters
     *
     * @param dim Dimension of the interval.
     * @param b Start value of the interval.
     * @param d End value of the interval.
     */
    Filtration_value_interval(int dim, filtration_value b, filtration_value d) : Base(dim, b, d) {}

    /**
     * @brief Returns the absolute length of the interval \f$|d-b|\f$.
     */
    filtration_value length() const {
      if (Base::b_ == Base::d_) {
        return 0;
      }  // otherwise inf - inf would return nan.
      return std::abs(Base::b_ - Base::d_);
    }
    /**
     * @brief Returns the absolute length of the log values of birth and death, i.e.  \f$|\log d - \log b|\f$.
     */
    filtration_value log_length() const {
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
   * 'birthToPos_[i] < birthToPos_[j]', with \f$0 <= i,j <= k\f$ indices in the quiver
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
  struct Birth_ordering {
    /**
     * @brief Default constructor
     */
    Birth_ordering() : birthToPos_(), maxBirthPos_(0), minBirthPos_(-1) {}

    /**
     * @brief Inserts arrow number in the ordering after an insertion.
     * When the arrow key-1 -> key is forward, key is larger than any other index
     * i < key in the birth ordering <b. We give key the largest value max_birth_pos_.
     *
     * @param arrow_number Forward arrow number.
     */
    void add_birth_forward(internal_key arrow_number) {  // amortized constant time
      birthToPos_.emplace_hint(birthToPos_.end(), arrow_number, maxBirthPos_);
      ++maxBirthPos_;
    }
    /**
     * @brief Inserts arrow number in the ordering after a removal.
     * When the arrow key-1 <- key is backward, key is smaller than any other index
     * i < key in the birth ordering <b. We give key the smallest value min_birth_pos_
     *
     * @param arrow_number Backward arrow number.
     */
    void add_birth_backward(internal_key arrow_number) {  // amortized constant time
      birthToPos_.emplace_hint(birthToPos_.end(), arrow_number, minBirthPos_);
      --minBirthPos_;
    }

    /**
     * @brief Removes the birth from the ordering.
     * When the row at index @p birth is removed from the homology matrix, we do not need
     * to maintain its position in <b anymore.
     *
     * @param birth Birth to remove.
     */
    void remove_birth(internal_key birth) { birthToPos_.erase(birth); }
    /**
     * @brief Increasing birth order <=b, true iff k1 <b k2.
     *
     * @param k1
     * @param k2
     * @return true if k1 <b k2, false otherwise.
     */
    bool birth_order(internal_key k1, internal_key k2) const { return birthToPos_.at(k1) < birthToPos_.at(k2); }
    /**
     * @brief Decreasing birth order <=b, true iff k1 >b k2.
     *
     * @param k1
     * @param k2
     * @return true if k1 >b k2, false otherwise.
     */
    bool reverse_birth_order(internal_key k1, internal_key k2) const {
      return birthToPos_.at(k1) > birthToPos_.at(k2);
    }

   private:
    std::unordered_map<internal_key, internal_key> birthToPos_; /**< birth_to_pos_[i] < birth_to_pos_[j] iff i <b j */
    internal_key maxBirthPos_;                                  /**< is strictly larger than any other birth so far */
    internal_key minBirthPos_;                                  /**< is strictly smaller than any other birth so far */
  };

  using Matrix_type = Gudhi::persistence_matrix::Matrix<Matrix_options>;
  using index = typename Matrix_type::index;

 public:
  /**
   * @brief Constructor of the Zigzag_persistence class.
   * @details After construction of the class, the zigzag filtration should be given in a streaming like way, i.e.,
   * call @ref insert_face or @ref remove_face for each step of the filtration in order of the filtration.
   * To retrieve the current persistence diagram at any moment of the filtration,
   * use @ref get_persistence_diagram or @ref get_index_persistence_diagram.
   *
   * @param minNumberOfFaces Minimum number of faces that will be inserted at some point in the filtration.
   * If the total number of faces is known in advance, the memory allocation can be better optimized.
   * Default value: 0.
   * @param ignoreCyclesAboveDim Ignores cycles in dimension larger or equal in the final diagram.
   * If -1, no cycles are ignored. Default value: -1.
   */
  Zigzag_persistence(unsigned int minNumberOfFaces = 0, int ignoreCyclesAboveDim = -1)
      : dimMax_(ignoreCyclesAboveDim),
        matrix_(
            minNumberOfFaces,
            [this](index columnIndex1, index columnIndex2) -> bool {
              if (matrix_.get_column(columnIndex1).is_paired()) {
                return matrix_.get_pivot(columnIndex1) < matrix_.get_pivot(columnIndex2);
              }
              return birthOrdering_.birth_order(births_.at(columnIndex1), births_.at(columnIndex2));
            },
            [this](index columnIndex1, index columnIndex2) -> bool { return false; }),
        numArrow_(-1),
        previousFiltrationValue_(std::numeric_limits<filtration_value>::infinity()) {}

  /**
   * @brief Updates the zigzag persistence diagram after the insertion of the given face.
   * 
   * @tparam BoundaryRange Range type needing begin and end members.
   * @param faceID ID representing the inserted face.
   * @param boundary Boundary of the inserted face. The range should be composed of the IDs of all faces contained in
   * the boundary (i.e. with non-zero coefficients), using the ID specified as `faceID` when the corresponding face
   * was previously inserted (recall that the faces should be inserted in order of filtration).
   * @param dimension Dimension of the inserted face.
   * @param filtrationValue Filtration value associated to the face.
   * Assumed to be larger or equal to previously used filtration values.
   */
  template <class BoundaryRange = std::initializer_list<face_key> >
  void insert_face(face_key faceID, 
                   const BoundaryRange& boundary, 
                   dimension_type dimension,
                   filtration_value filtrationValue) {
    if (dimMax_ != -1 && dimension > dimMax_) return;

    ++numArrow_;

    //TODO: to make it really stream like, we should stream out finished bars and remove unnecessary filtration values
    //from memory.
    if (filtrationValue != previousFiltrationValue_)  // check whether the filt value has changed
    {  // consecutive pairs (i,f), (j,f') mean faces of index k in [i,j-1] have
      previousFiltrationValue_ = filtrationValue;  // filtration value f
      filtrationValues_.emplace_back(numArrow_, previousFiltrationValue_);
    }

    [[maybe_unused]] auto res = handleToKey_.try_emplace(faceID, numArrow_);

    GUDHI_CHECK(res.second, "Zigzag_persistence::insert_face - face already in the complex");

    // Reduce the boundary of zzsh in the basis of cycles.
    // Compute the keys of the faces of the boundary of zzsh.
    std::set<internal_key> col_bsh;  // set maintains the natural order on indices
    for (auto b_sh : boundary) {
      col_bsh.insert(handleToKey_.at(b_sh));  // TODO: add possibilities of coefficients
    }

    _process_forward_arrow(col_bsh, dimension);
  }

  /**
   * @brief Updates the zigzag persistence diagram after the removal of the given face.
   * 
   * @param faceID ID representing the face to remove. Should be the same than the one used to insert it.
   * @param dimension Dimension of the face.
   * @param filtrationValue Filtration value associated to the removal.
   * Assumed to be larger or equal to previously used filtration values.
   */
  void remove_face(face_key faceID, dimension_type dimension, filtration_value filtrationValue) {
    if (dimMax_ != -1 && dimension > dimMax_) return;

    ++numArrow_;

    auto it = handleToKey_.find(faceID);
    GUDHI_CHECK(it != handleToKey_.end(), "Zigzag_persistence::remove_face - face not in the complex");

    if (filtrationValue != previousFiltrationValue_)  // check whether the filt value has changed
    {  // consecutive pairs (i,f), (j,f') mean faces of index k in [i,j-1] have
      previousFiltrationValue_ = filtrationValue;  // filtration value f
      filtrationValues_.emplace_back(numArrow_, previousFiltrationValue_);
    }

    _process_backward_arrow(it->second, dimension);
    handleToKey_.erase(it);
  }

  /**
   * @brief Returns the "index persistence diagram" of the current filtration, that is, the pairs of atomic arrow
   * numbers corresponding to a birth-death pair. Does not contain points at infinity, only the cycle classes which
   * already died are represented.
   *
   * @return Reference to the list of intervals.
   */
  const std::list<Index_interval>& get_index_persistence_diagram() const { return persistenceDiagram_; }

  /**
   * @brief Returns the filtration values \f$[f(b),f(d)]\f$ associated to the indices \f$[b,d]\f$ which are retrieved
   * by @ref get_index_persistence_diagram.
   *
   * @param birthKey Birth index
   * @param deathKey Death index
   * @return A pair of filtration values associated to the given indices.
   */
  std::pair<filtration_value, filtration_value> map_index_to_filtration_value(internal_key birthKey,
                                                                              internal_key deathKey) const {
    // filtration_values_ must be sorted by increasing keys.
    auto itBirth =  // lower_bound(x) returns leftmost y s.t. x <= y
        std::lower_bound(
            filtrationValues_.begin(), filtrationValues_.end(),
            std::pair<internal_key, filtration_value>(birthKey, std::numeric_limits<filtration_value>::infinity()),
            [](std::pair<internal_key, filtration_value> p1, std::pair<internal_key, filtration_value> p2) {
              return p1.first < p2.first;
            });
    if (itBirth == filtrationValues_.end() || itBirth->first > birthKey) {
      --itBirth;
    }
    // it points to the rightmost z such that z <= x

    auto itDeath =  //
        std::lower_bound(
            filtrationValues_.begin(), filtrationValues_.end(),
            std::pair<internal_key, filtration_value>(deathKey, std::numeric_limits<filtration_value>::infinity()),
            [](std::pair<internal_key, filtration_value> p1, std::pair<internal_key, filtration_value> p2) {
              return p1.first < p2.first;
            });
    if (itDeath == filtrationValues_.end() || itDeath->first > deathKey) {
      --itDeath;
    }

    return std::make_pair(itBirth->second, itDeath->second);
  }

  /**
   * @brief Returns the current persistence diagram ordered first by length, than by dimension,
   * than by birth value and finally by death value.
   *
   * @param shortestInterval Threshold. Every bar shorter than the given value will be ignored. Default value: 0.
   * @param includeInfinitBars If set to true, infinit bars are included in the diagram. Default value: false.
   * @return A vector of pairs of filtration values representing the persistence diagram.
   */
  std::vector<Filtration_value_interval> get_persistence_diagram(filtration_value shortestInterval = 0.,
                                                                 bool includeInfinitBars = false) {
    auto comp = [](Filtration_value_interval p, Filtration_value_interval q) {
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

    std::vector<Filtration_value_interval> diag = _get_persistence_diagram(shortestInterval);

    if (includeInfinitBars) {
      _retrieve_infinit_bars(diag);
    }

    std::stable_sort(diag.begin(), diag.end(), comp);

    return diag;
  }

 private:

  /**
   * @brief Express the boundary cycle of the new face as a sum of cycles in a matrix.
   * If some cycles are not boundary cycles, i.e., columns with F-index
   * in the matrix, it applies a surjective diamond to the zigzag module.
   * 
   * @param boundary Boundary of the inserted face.
   * @param dim Dimension of the inserted face.
   */
  void _process_forward_arrow(const std::set<internal_key>& boundary, dimension_type dim) {
    std::vector<index> chainsInF = matrix_.insert_boundary(numArrow_, boundary);

    if (!chainsInF.empty()) {
      _apply_surjective_reflection_diamond(dim, chainsInF);
    } else {
      birthOrdering_.add_birth_forward(numArrow_);
      births_.emplace_hint(births_.end(), matrix_.get_column_with_pivot(numArrow_), numArrow_);
    }
  }

  /**
   * @brief Applies the surjective reflection diamond principle to the current filtration.
   *
   * @details The vector chainsInF is sorted by decreasing lowest index values in the
   * columns corresponding to the chains, due to its computation in the reduction of
   * the boundary in _process_forward_arrow(...). It is equivalent to decreasing death index
   * order w.r.t. the <d ordering.
   * 
   * @param dim Dimension of the inserted face.
   * @param chainsInF Indices of the non paired columns in the matrix.
   */
  void _apply_surjective_reflection_diamond(dimension_type dim, const std::vector<index>& chainsInF) {
    // fp is the largest death index for <=d
    // Set col_fp: col_fp <- col_f1+...+col_fp (now in G); preserves lowest idx
    auto chainFp = chainsInF[0];  // col_fp, with largest death <d index.

    // chainsInF is ordered, from .begin() to end(), by decreasing lowest_idx_. The
    // lowest_idx_ is also the death of the chain in the right suffix of the
    // filtration (all backward arrows). Consequently, the chains in F are ordered by
    // decreasing death for <d.
    // Pair the col_fi, i = 1 ... p-1, according to the reflection diamond principle
    // Order the fi by reverse birth ordering <=_b
    auto cmp_birth = [this](internal_key k1, internal_key k2) -> bool {
      return birthOrdering_.reverse_birth_order(k1, k2);
    };  // true iff b(k1) >b b(k2)

    // available_birth: for all i by >d value of the d_i,
    // contains at step i all b_j, j > i, and maybe b_i if not stolen
    std::set<internal_key, decltype(cmp_birth)> availableBirth(cmp_birth);
    // for f1 to f_{p} (i by <=d), insertion in available_birth_to_fidx sorts by >=b
    for (auto& chainF : chainsInF) {
      availableBirth.insert(births_.at(chainF));
    }

    auto maxbIt = availableBirth.begin();  // max birth cycle
    auto maxb = *maxbIt;                   // max birth value, for persistence diagram
    availableBirth.erase(maxbIt);          // remove max birth cycle (stolen)

    auto lastModifiedChainIt = chainsInF.rbegin();

    // consider all death indices by increasing <d order i.e., increasing lowest_idx_
    for (auto chainFIt = chainsInF.rbegin();    // by increasing death order <d
         *chainFIt != chainFp; ++chainFIt)      // chain_fp=*begin() has max death
    {                                           // find which reduced col has this birth
      auto birthIt = availableBirth.find(births_.at(*chainFIt));
      if (birthIt == availableBirth.end())  // birth is not available. *chain_f_it
      {                                     // must become the sum of all chains in F with smaller death index.
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
        for (auto chainPassedIt = lastModifiedChainIt; chainPassedIt != chainFIt; ++chainPassedIt) {
          // all with smaller <d death
          matrix_.add_to(*chainPassedIt, *chainFIt);
        }
        lastModifiedChainIt = chainFIt;  // new cumulated c_i+...+c_1
        // remove the max available death
        auto maxAvailBirthIt = availableBirth.begin();  // max because order by deacr <b
        internal_key maxAvailBirth = *maxAvailBirthIt;      // max available birth

        births_.at(*chainFIt) = maxAvailBirth;  // give new birth
        availableBirth.erase(maxAvailBirthIt);  // remove birth from availability
      } else {
        availableBirth.erase(birthIt);
      }  // birth not available anymore, do not
    }    // modify *chain_f_it.

    // Following value can be erased, but it slowes the process down a bit, so I keep it as a remainder for now:
    //  birth_ordering_.remove_birth(maxb);
    births_.erase(chainFp);

    // Update persistence diagram with left interval [fil(b_max) ; fil(m))
    persistenceDiagram_.emplace_back(dim - 1, maxb, numArrow_);  //-1);//
  }

  /**
   * @brief Removes the given face by pushing up the matrix the corresponding column and erasing it.
   * 
   * @param faceID Internal ID of the face to remove. 
   * @param dim Dimension of the face to remove.
   */
  void _process_backward_arrow(internal_key faceID, dimension_type dim) {
    // column whose key is the one of the removed face
    index currCol = matrix_.get_column_with_pivot(faceID);

    // Record all columns that get affected by the transpositions, i.e., have a coeff
    std::vector<index> modifiedColumns;
    const auto& row = matrix_.get_row(faceID);
    modifiedColumns.reserve(row.size());
    std::transform(row.begin(), row.end(), std::back_inserter(modifiedColumns),
                   [](const auto& cell) { return cell.get_column_index(); });
    // Sort by left-to-right order in the matrix_ (no order maintained in rows)
    std::stable_sort(modifiedColumns.begin(), modifiedColumns.end(),
                     [this](index i1, index i2) { return matrix_.get_pivot(i1) < matrix_.get_pivot(i2); });

    // Modifies curr_col, not the other one.
    for (auto otherColIt = std::next(modifiedColumns.begin()); otherColIt != modifiedColumns.end(); ++otherColIt) {
      currCol = matrix_.vine_swap_with_z_eq_1_case(currCol, *otherColIt);
    }

    // curr_col points to the column to remove by restriction of K to K-{\sigma}
    if (!matrix_.get_column(currCol).is_paired()) {  // in F
      auto it = births_.find(currCol);
      if (dimMax_ == -1 || (dimMax_ != -1 && dim < dimMax_)) {        // don't record intervals over max dim
        persistenceDiagram_.emplace_back(dim, it->second, numArrow_);
      }
      // Following value can be erased, but it slowes the process down a bit, so I keep it as a remainder for now:
      //  birthOrdering_.remove_birth(it->second);
      births_.erase(it);
    } else {  // in H    -> paired with c_g, that now belongs to F now
      // maintain the <=b order
      birthOrdering_.add_birth_backward(numArrow_);
      births_.try_emplace(matrix_.get_column(currCol).get_paired_chain_index(), numArrow_);
    }

    // cannot be in G as the removed face is maximal
    matrix_.remove_maximal_face(faceID, {});
  }

  /**
   * @brief Returns the current persistence diagram ordered by length without infinit bars.
   *
   * @param shortestInterval Intervals shorter than the given value are ignored.
   * @return Vector of intervals.
   */
  std::vector<Filtration_value_interval> _get_persistence_diagram(filtration_value shortestInterval) {
    std::vector<Filtration_value_interval> diag;
    diag.reserve(persistenceDiagram_.size());

    std::stable_sort(filtrationValues_.begin(), filtrationValues_.end(),
                     [](std::pair<internal_key, filtration_value> p1, std::pair<internal_key, filtration_value> p2) {
                       return p1.first < p2.first;
                     });

    for (auto bar : persistenceDiagram_) {
      filtration_value birth, death;
      std::tie(birth, death) = map_index_to_filtration_value(bar.birth(), bar.death());
      if (birth > death) {
        std::swap(birth, death);
      }

      if (death - birth > shortestInterval) {
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
  void _retrieve_infinit_bars(std::vector<Filtration_value_interval>& diag) const {
    auto birth = [this](internal_key birthKey) {
      auto itBirth =  // lower_bound(x) returns leftmost y s.t. x <= y
          std::lower_bound(
              filtrationValues_.begin(), filtrationValues_.end(),
              std::pair<internal_key, filtration_value>(birthKey, std::numeric_limits<filtration_value>::infinity()),
              [](std::pair<internal_key, filtration_value> p1, std::pair<internal_key, filtration_value> p2) {
                return p1.first < p2.first;
              });
      if (itBirth == filtrationValues_.end() || itBirth->first > birthKey) {
        --itBirth;
      }
      return itBirth->second;
    };

    for (auto& p : births_) {
      auto dim = matrix_.get_column(matrix_.get_column_with_pivot(p.first)).get_dimension();
      if (dimMax_ == -1 || (dimMax_ != -1 && dim < dimMax_))
        diag.emplace_back(dim, birth(p.second), std::numeric_limits<filtration_value>::infinity());
    }
  }

 private:
  std::unordered_map<face_key,internal_key> handleToKey_; /**< Map from input keys to internal keys. */
  dimension_type dimMax_;                                 /**< Maximal dimension of a bar to record. */
  Matrix_type matrix_;                                    /**< Matrix storing a base of the current chain complex. */
  std::unordered_map<index, int> births_;                 /**< Map face index in F to corresponding birth. */
  Birth_ordering birthOrdering_;                          /**< Maintains <b ordering of the births. */
  std::list<Index_interval> persistenceDiagram_;          /**< Stores current closed persistence intervals. */
  internal_key numArrow_;                                 /**< Current arrow number. */
  filtration_value previousFiltrationValue_;              /**< Filtration value of the previous arrow. */
  /**
   * @brief filtrationValues_ stores consecutive pairs (i,f) , (j,f') with f != f',
   * meaning that all inserted faces with key in [i;j-1] have filtration value f,
   * i is the smallest face index whose face has filtration value f.
   */
  std::vector<std::pair<internal_key, filtration_value> > filtrationValues_;
};  // end class Zigzag_persistence

}  // namespace zigzag_persistence

}  // namespace Gudhi

#endif  // ZIGZAG_PERSISTENCE_H_
