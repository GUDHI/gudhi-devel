/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2021 Inria
 *
 *    Modification(s):
 *      - 2023/05 Hannah Schreiber: Rework of the interface, reorganization and debug
 *      - 2023/05 Hannah Schreiber: Addition of infinite bars
 *      - 2024/06 Hannah Schreiber: Separation of the zigzag algorithm from the filtration value management
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file zigzag_persistence.h
 * @author Clément Maria, Hannah Schreiber
 * @brief Contains the implementation of the @ref Gudhi::zigzag_persistence::Zigzag_persistence class.
 */

#ifndef ZIGZAG_PERSISTENCE_H_
#define ZIGZAG_PERSISTENCE_H_

#include <cmath>
#include <set>
#include <utility>
#include <vector>

#include <boost/iterator/indirect_iterator.hpp>
#if BOOST_VERSION >= 108100
#include <boost/unordered/unordered_flat_map.hpp>   //don't exist for lower versions of boost
// #include <boost/unordered/unordered_map.hpp>
#else
#include <unordered_map>
#endif

#include <gudhi/Matrix.h>

namespace Gudhi {
namespace zigzag_persistence {

/**
 * @ingroup zigzag_persistence
 *
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
 * @ingroup zigzag_persistence
 *
 * @brief Default options for @ref Zigzag_persistence.
 */
struct Default_zigzag_options {
  using internal_key = int;   /**< Face ID used internally, must be signed. */
  using dimension_type = int; /**< Dimension value type. */
  /**
   * @brief Column type use by the internal matrix.
   */
  static const Gudhi::persistence_matrix::Column_types column_type =
      Gudhi::persistence_matrix::Column_types::NAIVE_VECTOR;  //TODO: redo benchmark with oscillating rips
};

// TODO: add the possibility of something else than Z2. Which means that the possibility of vineyards without Z2
// also needs to be implemented. The theory needs to be done first.
// TODO: erase_birth_history will be moved to the options if it is proven to be useful. In the meantime
// it stays here undocumented to ease benchmarks.
/**
 * @class Zigzag_persistence zigzag_persistence.h gudhi/zigzag_persistence.h
 * @brief Class computing the zigzag persistent homology of a zigzag sequence. Algorithm based on \cite zigzag.
 * @details After construction of the class, the zigzag filtration should be given in a streaming like way, i.e.,
 * call @ref insert_face, @ref remove_face or @ref apply_identity for each step of the filtration in order of
 * the filtration. The pairs of birth and death indices are retrieved via the given callback method every time
 * a pair is closed. To retrieve the open pairs (corresponding to infinite bars),
 * use @ref get_current_infinite_intervals.
 *
 * ### Minimalistic example of usage
 *
 * #### Includes
 * ```
 * #include <gudhi/zigzag_persistence.h>
 * ```
 *
 * #### Useful aliases
 * ```
 * using Zigzag_persistence = Gudhi::zigzag_persistence::Zigzag_persistence<>;
 *
 * using dimension_type = Zigzag_persistence::dimension_type;
 * using index = Zigzag_persistence::index;
 * ```
 *
 * #### Construction with default values
 * ```
 * //Zigzag_persistence(callback) with for example callback method as a anonymous lambda
 * Zigzag_persistence zp([](dimension_type dim, index birth, index death) {
 *   std::cout << "[" << dim << "] " << birth << " - " << death << std::endl;
 * });
 * ```
 *
 * #### Input of the zigzag sequence/filtration
 * ```
 * // In all cases, it is important that the operations of insertions and removals are made **in the same order**
 * // as in the zigzag filtration ones wants to compute the barcode from.
 *
 * // A face has to be identified in the boundaries by the operation number the face was inserted with in the sequence.
 *
 * //inserts vertex 0 -> birth at 0 of 0-cycle
 * zp.insert_face({}, 0);
 * //inserts vertex 1 -> birth at 1 of 0-cycle
 * zp.insert_face({}, 0);
 * //inserts edge 2 = (0,1) -> death at 2 -> outputs (0, 1, 2)
 * zp.insert_face({0, 1}, 1);
 * //inserts vertex 3 -> birth at 3 of 0-cycle
 * zp.insert_face({}, 0);
 * //inserts edge 4 = (0,3) -> death at 4 -> outputs (0, 3, 4)
 * zp.insert_face({0, 3}, 1);
 * //inserts edge 5 = (1,3) -> birth at 5 of 1-cycle
 * zp.insert_face({1, 3}, 1);
 * //removes edge 4 -> death at 6 -> outputs (1, 5, 6)
 * zp.remove_face(4);
 * //removes edge 2 -> birth at 7 of 0-cycle
 * zp.remove_face(2);
 * ```
 *
 * #### Finalizations
 * ```
 * // Only the closed bars where output so far, so the open/infinite bars still need to be retrieved.
 *
 * //in this example, outputs (0, 0) and (0, 7)
 * zp.get_current_infinite_intervals([](dimension_type dim, index birth){
 *   std::cout << "[" << dim << "] " << birth << " - inf" << std::endl;
 * });
 * ```
 *
 * @ingroup zigzag_persistence
 *
 * @tparam ZigzagOptions Structure following the @ref ZigzagOptions concept. Default value: @ref Default_zigzag_options.
 */
template <class ZigzagOptions = Default_zigzag_options, bool erase_birth_history = true>
class Zigzag_persistence
{
 public:
  using Options = ZigzagOptions;                           /**< Zigzag options. */
  using index = typename Options::internal_key;            /**< Key and index type, has to be signed. */
  using dimension_type = typename Options::dimension_type; /**< Type for dimension values. */

 private:
#if BOOST_VERSION >= 108100
  using birth_dictionary = boost::unordered_flat_map<index, index>;      /**< Dictionary type. */
  // using birth_dictionary = boost::unordered_map<index, index>;           /**< Dictionary type. */
#else
  using birth_dictionary = std::unordered_map<index, index>;             /**< Dictionary type. */
#endif
  using Matrix_options = Zigzag_matrix_options<Options::column_type>;     /**< Matrix options. */
  using Matrix_type = Gudhi::persistence_matrix::Matrix<Matrix_options>;  /**< Matrix. */
  using matrix_index = typename Matrix_type::Index;                       /**< Matrix indexation type. */

  /** \brief Maintains the birth ordering \f$\leq_b\f$.
   *
   * \details Contains a map of size the number of
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
     * @param arrowNumber Forward arrow number.
     */
    void add_birth_forward(index arrowNumber) {  // amortized constant time
      birthToPos_.emplace_hint(birthToPos_.end(), arrowNumber, maxBirthPos_);
      ++maxBirthPos_;
    }
    /**
     * @brief Inserts arrow number in the ordering after a removal.
     * When the arrow key-1 <- key is backward, key is smaller than any other index
     * i < key in the birth ordering <b. We give key the smallest value min_birth_pos_
     *
     * @param arrowNumber Backward arrow number.
     */
    void add_birth_backward(index arrowNumber) {  // amortized constant time
      birthToPos_.emplace_hint(birthToPos_.end(), arrowNumber, minBirthPos_);
      --minBirthPos_;
    }

    /**
     * @brief Removes the birth from the ordering.
     * When the row at index @p birth is removed from the homology matrix, we do not need
     * to maintain its position in <b anymore.
     *
     * @param birth Birth to remove.
     */
    void remove_birth(index birth) { birthToPos_.erase(birth); }
    /**
     * @brief Increasing birth order <=b, true iff k1 <b k2.
     *
     * @param k1
     * @param k2
     * @return true if k1 <b k2, false otherwise.
     */
    bool birth_order(index k1, index k2) const { return birthToPos_.at(k1) < birthToPos_.at(k2); }
    /**
     * @brief Decreasing birth order <=b, true iff k1 >b k2.
     *
     * @param k1
     * @param k2
     * @return true if k1 >b k2, false otherwise.
     */
    bool reverse_birth_order(index k1, index k2) const { return birthToPos_.at(k1) > birthToPos_.at(k2); }

   private:
    birth_dictionary birthToPos_; /**< birth_to_pos_[i] < birth_to_pos_[j] iff i <b j */
    index maxBirthPos_;            /**< is strictly larger than any other birth so far */
    index minBirthPos_;            /**< is strictly smaller than any other birth so far */
  };

 public:
  /**
   * @brief Constructor of the Zigzag_persistence class.
   * @details After construction of the class, the zigzag filtration should be given in a streaming like way, i.e.,
   * call @ref insert_face, @ref remove_face or @ref apply_identity for each step of the filtration in order of
   * the filtration. The pairs of birth and death indices are retrieved via the given callback method every time
   * a pair is closed. To retrieve the open pairs (corresponding to infinite bars),
   * use @ref get_current_infinite_intervals.
   *
   * @param stream_interval Callback method to process the birth and death index pairs. Has to take three arguments
   * as input: first the dimension of the cycle, then the birth index of the cycle and third the death index of the
   * cycle. An index always corresponds to the arrow number the event occurred (one call to @ref insert_face,
   * @ref remove_face or @ref apply_identity is equal to one arrow and increases the arrow count by one).
   * @param preallocationSize Reserves space for @p preallocationSize faces in the internal data structure.
   * This is optional and just helps skip a few reallocations. The optimal value (no reallocation, no wasted space) is
   * the number of faces in the biggest complex of the filtration.
   * Default value: 0.
   */
  Zigzag_persistence(std::function<void(dimension_type, index, index)> stream_interval,
                     unsigned int preallocationSize = 0)
      : matrix_(
            preallocationSize,
            [this](matrix_index columnIndex1, matrix_index columnIndex2) -> bool {
              if (matrix_.get_column(columnIndex1).is_paired()) {
                return matrix_.get_pivot(columnIndex1) < matrix_.get_pivot(columnIndex2);
              }
              return birthOrdering_.birth_order(births_.at(columnIndex1), births_.at(columnIndex2));
            },
            [this](matrix_index columnIndex1, matrix_index columnIndex2) -> bool { return false; }),
        numArrow_(-1),
        stream_interval_(std::move(stream_interval)) {}
  /**
   * @brief Updates the zigzag persistence diagram after the insertion of the given face.
   *
   * @tparam BoundaryRange Range type needing size, begin and end members.
   * @param boundary Boundary of the inserted face. The boundary should be represented by all the faces with
   * non-zero coefficients generating it. A face should be represented by the arrow number when the face appeared for
   * the first time in the filtration (if a face was inserted and then removed and reinserted etc., only the last
   * insertion counts). The face range should be ordered by increasing arrow numbers.
   * @param dimension Dimension of the inserted face.
   * @return Number of the operation.
   */
  template <class BoundaryRange = std::initializer_list<index> >
  index insert_face(const BoundaryRange& boundary, dimension_type dimension) {
    ++numArrow_;
    _process_forward_arrow(boundary, dimension);
    return numArrow_;
  }

  /**
   * @brief Updates the zigzag persistence diagram after the removal of the given face.
   *
   * @param arrowNumber Arrow number of when the face to remove was inserted for the last time.
   * @return Number of the operation.
   */
  index remove_face(index arrowNumber) {
    ++numArrow_;
    _process_backward_arrow(arrowNumber);
    return numArrow_;
  }

  /**
   * @brief To use when a face is neither inserted nor removed, but the filtration moves along the identity operator
   * on homology level. Useful to keep the birth/death indices aligned when insertions/removals are purposely skipped
   * to avoid useless computation. Increases the arrow number by one.
   * @return Number of the operation.
   */
  index apply_identity() { return ++numArrow_; }

  /**
   * @brief Outputs through the given callback method all birth indices which are currently not paired with 
   * a death index.
   * 
   * @tparam F Type of the callback method. Takes two arguments: the dimension of the cycle and the birth index
   * of the cycle.
   * @param stream_infinite_interval Method processing the unpaired birth indices.
   */
  template <typename F>
  void get_current_infinite_intervals(F&& stream_infinite_interval) {
    for (auto& p : births_) {
      if constexpr (erase_birth_history) {
        auto& col = matrix_.get_column(p.first);
        stream_infinite_interval(col.get_dimension(), p.second);
      } else {
        try {
          auto& col = matrix_.get_column(p.first);
          if (!col.is_paired()) {
            stream_infinite_interval(col.get_dimension(), p.second);
          }
        } catch (const std::out_of_range&) {
          continue;
        }
      }
    }
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
  template <class BoundaryRange>
  void _process_forward_arrow(const BoundaryRange& boundary, dimension_type dim) {
    std::vector<matrix_index> chainsInF = matrix_.insert_boundary(numArrow_, boundary, dim);

    if (!chainsInF.empty()) {
      _apply_surjective_reflection_diamond(dim, chainsInF);
    } else {
      birthOrdering_.add_birth_forward(numArrow_);
      births_[matrix_.get_column_with_pivot(numArrow_)] = numArrow_;
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
  void _apply_surjective_reflection_diamond(dimension_type dim, const std::vector<matrix_index>& chainsInF) {
    // fp is the largest death index for <=d
    // Set col_fp: col_fp <- col_f1+...+col_fp (now in G); preserves lowest idx
    auto chainFp = chainsInF[0];  // col_fp, with largest death <d index.

    // chainsInF is ordered, from .begin() to end(), by decreasing lowest_idx_. The
    // lowest_idx_ is also the death of the chain in the right suffix of the
    // filtration (all backward arrows). Consequently, the chains in F are ordered by
    // decreasing death for <d.
    // Pair the col_fi, i = 1 ... p-1, according to the reflection diamond principle
    // Order the fi by reverse birth ordering <=_b
    auto cmp_birth = [this](index k1, index k2) -> bool {
      return birthOrdering_.reverse_birth_order(k1, k2);
    };  // true iff b(k1) >b b(k2)

    // availableBirth: for all i by >d value of the d_i,
    // contains at step i all b_j, j > i, and maybe b_i if not stolen
    std::set<index, decltype(cmp_birth)> availableBirth(cmp_birth);
    // for f1 to f_{p} (i by <=d), insertion in availableBirth sorts by >=b
    for (auto& chainF : chainsInF) {
      availableBirth.insert(births_.at(chainF));
    }

    auto maxbIt = availableBirth.begin();  // max birth cycle
    auto maxb = *maxbIt;                   // max birth value, for persistence diagram
    availableBirth.erase(maxbIt);          // remove max birth cycle (stolen)

    auto lastModifiedChainIt = chainsInF.rbegin();

    // consider all death indices by increasing <d order i.e., increasing lowest_idx_
    for (auto chainFIt = chainsInF.rbegin();  // by increasing death order <d
         *chainFIt != chainFp; ++chainFIt)    // chain_fp=*begin() has max death
    {                                         // find which reduced col has this birth
      auto birthIt = availableBirth.find(births_.at(*chainFIt));
      if (birthIt == availableBirth.end())  // birth is not available. *chain_f_it
      {                                     // must become the sum of all chains in F with smaller death index.
        // this gives as birth the maximal birth of all chains with strictly larger
        // death <=> the maximal available death.
        // Let c_1 ... c_f be the chains s.t. <[c_1+...+c_f]> is the kernel and
        //  death(c_i) >d death(c_i-1). If the birth of c_i is not available, we set
        // c_i <- c_i + c_i-1 + ... + c_1, which is [c_i + c_i-1 + ... + c_1] on
        // the right (of death the maximal<d death(c_i)), and is [c_i + c_i-1 + ... +
        // c_1] + kernel = [c_f + c_f-1 + ... + c_i+1] on the left (of birth the max<b
        // of the birth of the c_j, j>i  <=> the max<b available birth).
        // N.B. some of the c_k, k<i, have already been modified to be equal to
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
        auto maxAvailBirthIt = availableBirth.begin();  // max because order by decreasing <b
        index maxAvailBirth = *maxAvailBirthIt;         // max available birth

        births_.at(*chainFIt) = maxAvailBirth;  // give new birth
        availableBirth.erase(maxAvailBirthIt);  // remove birth from availability
      } else {
        availableBirth.erase(birthIt);
      }  // birth not available anymore, do not
    }  // modify *chain_f_it.

    if constexpr (erase_birth_history) {
      birthOrdering_.remove_birth(maxb);
      births_.erase(chainFp);
    }

    // Update persistence diagram with left interval [fil(b_max) ; fil(m))
    stream_interval_(dim - 1, maxb, numArrow_);
  }

  /**
   * @brief Removes the given face by pushing up the matrix the corresponding column and erasing it.
   *
   * @param faceID Internal ID of the face to remove.
   */
  void _process_backward_arrow(index faceID) {
    // column whose key is the one of the removed face
    matrix_index currCol = matrix_.get_column_with_pivot(faceID);

    // Record all columns that get affected by the transpositions, i.e., have a coeff
    std::vector<matrix_index> modifiedColumns;
    const auto& row = matrix_.get_row(faceID);
    modifiedColumns.reserve(row.size());
    std::transform(row.begin(), row.end(), std::back_inserter(modifiedColumns),
                   [](const auto& cell) { return cell.get_column_index(); });
    // Sort by left-to-right order in the matrix_ (no order maintained in rows)
    std::stable_sort(modifiedColumns.begin(), modifiedColumns.end(), [this](matrix_index i1, matrix_index i2) {
      return matrix_.get_pivot(i1) < matrix_.get_pivot(i2);
    });

    // Modifies curr_col, not the other one.
    for (auto otherColIt = std::next(modifiedColumns.begin()); otherColIt != modifiedColumns.end(); ++otherColIt) {
      currCol = matrix_.vine_swap_with_z_eq_1_case(currCol, *otherColIt);
    }

    // curr_col points to the column to remove by restriction of K to K-{\sigma}
    auto& col = matrix_.get_column(currCol);
    if (!col.is_paired()) {  // in F
      auto it = births_.find(currCol);
      stream_interval_(col.get_dimension(), it->second, numArrow_);
      if constexpr (erase_birth_history) {
        birthOrdering_.remove_birth(it->second);
        births_.erase(it);
      }
    } else {  // in H    -> paired with c_g, that now belongs to F now
      // maintain the <=b order
      birthOrdering_.add_birth_backward(numArrow_);
      births_[col.get_paired_chain_index()] = numArrow_;
    }

    // cannot be in G as the removed face is maximal
    matrix_.remove_maximal_face(faceID, {});  // also un-pairs c_g if in H
  }

 private:
  Matrix_type matrix_;           /**< Matrix storing a base of the current chain complex. */
  birth_dictionary births_;      /**< Map face index in F to corresponding birth. */
  Birth_ordering birthOrdering_; /**< Maintains <b ordering of the births. */
  index numArrow_;               /**< Current arrow number. */
  std::function<void(dimension_type, index, index)> stream_interval_; /**< Callback method for closed pairs. */
};  // end class Zigzag_persistence

}  // namespace zigzag_persistence
}  // namespace Gudhi

#endif  // ZIGZAG_PERSISTENCE_H_
