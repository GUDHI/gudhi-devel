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
 * @file filtered_zigzag_persistence.h
 * @author Clément Maria, Hannah Schreiber
 * @brief Contains the implementation of the @ref Gudhi::zigzag_persistence::Default_filtered_zigzag_options structure
 * and the @ref Gudhi::zigzag_persistence::Filtered_zigzag_persistence_with_storage and
 * @ref Gudhi::zigzag_persistence::Filtered_zigzag_persistence classes.
 */

#ifndef FILTERED_ZIGZAG_PERSISTENCE_H_
#define FILTERED_ZIGZAG_PERSISTENCE_H_

#include <cmath>
#include <limits>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

#include <gudhi/Debug_utils.h>
#include <gudhi/zigzag_persistence.h>
#include <gudhi/persistence_interval.h>

namespace Gudhi {
namespace zigzag_persistence {

/**
 * @ingroup zigzag_persistence
 *
 * @brief Default options for @ref Filtered_zigzag_persistence_with_storage and @ref Filtered_zigzag_persistence.
 */
struct Default_filtered_zigzag_options {
  using internal_key = int;        /**< Face ID used internally, must be signed. */
  using face_key = int;            /**< Face ID used in the given boundaries. */
  using filtration_value = double; /**< Filtration value type. */
  using dimension_type = int;      /**< Dimension value type. */
  /**
   * @brief Column type use by the internal matrix.
   */
  static const Gudhi::persistence_matrix::Column_types column_type =
      Gudhi::persistence_matrix::Column_types::NAIVE_VECTOR;
};

/**
 * @ingroup zigzag_persistence
 *
 * @brief Class computing the zigzag persistent homology of a zigzag filtration. Algorithm based on \cite zigzag.
 * Even though the insertions and removals are given in a "stream-like" way, the barcode and other values are
 * stored during the whole process and not removed. It is therefore suited for smaller filtrations where the clean
 * ups produce a higher overhead than the memory consumption.
 * @details After construction of the class, the zigzag filtration should be given in a streaming like way, i.e.,
 * call @ref insert_face, @ref remove_face or @ref apply_identity for each step of the filtration in order of
 * the filtration. To retrieve the current persistence diagram at any moment of the filtration,
 * use @ref get_persistence_diagram or @ref get_index_persistence_diagram.
 * 
 * @tparam FilteredZigzagOptions Structure following the @ref FilteredZigzagOptions concept.
 * Default value: @ref Default_filtered_zigzag_options.
 */
template <class FilteredZigzagOptions = Default_filtered_zigzag_options>
class Filtered_zigzag_persistence_with_storage
{
 public:
  using Options = FilteredZigzagOptions;                        /**< Zigzag options. */
  using internal_key = typename Options::internal_key;          /**< Key and index type, has to be signed. */
  using face_key = typename Options::face_key;                  /**< Face ID type from external inputs. */
  using filtration_value = typename Options::filtration_value;  /**< Type for filtration values. */
  using dimension_type = typename Options::dimension_type;      /**< Type for dimension values. */

  /**
   * @brief Persistence index interval type.
   */
  using Index_interval = Gudhi::persistence_matrix::Persistence_interval<dimension_type,internal_key>;
  /**
   * @brief Persistence filtration interval type.
   */
  using Filtration_value_interval = Gudhi::persistence_matrix::Persistence_interval<dimension_type,filtration_value>;

  /**
   * @brief Constructor.
   * @details After construction of the class, the zigzag filtration should be given in a streaming like way, i.e.,
   * call @ref insert_face, @ref remove_face or @ref apply_identity for each step of the filtration in order of
   * the filtration. To retrieve the current persistence diagram at any moment of the filtration,
   * use @ref get_persistence_diagram or @ref get_index_persistence_diagram.
   *
   * @param preallocationSize Space for @p preallocationSize faces are reserved in the underlying structure.
   * Theoretically, any values works therefore, but for better performances, it is better to be as close as possible
   * to the maximal value of the number of faces stored at the same time. At a same time are stored faces which were
   * inserted before that time but not removed until that time. Default value: 0.
   * @param ignoreCyclesAboveDim Ignores cycles in dimension larger or equal in the final diagram.
   * If -1, no cycles are ignored. Default value: -1.
   */
  Filtered_zigzag_persistence_with_storage(unsigned int preallocationSize = 0, int ignoreCyclesAboveDim = -1)
      : dimMax_(ignoreCyclesAboveDim),
        persistenceDiagram_(),
        numArrow_(-1),
        previousFiltrationValue_(std::numeric_limits<filtration_value>::infinity()),
        pers_(
            [&](dimension_type dim, internal_key birth, internal_key death) {
              if (dimMax_ == -1 || (dimMax_ != -1 && dim < dimMax_)) {  // don't record intervals over max dim
                persistenceDiagram_.emplace_back(birth, death, dim);
              }
            },
            preallocationSize) {}

  /**
   * @brief Updates the zigzag persistence diagram after the insertion of the given face.
   *
   * @tparam BoundaryRange Range type needing size, begin and end members.
   * @param faceID ID representing the inserted face.
   * @param boundary Boundary of the inserted face. The range should be composed of the IDs of all faces contained in
   * the boundary (i.e. with non-zero coefficients), using the ID specified as `faceID` when the corresponding face
   * was previously inserted (recall that the faces should be inserted in order of filtration).
   * @param dimension Dimension of the inserted face.
   * @param filtrationValue Filtration value associated to the face.
   * Assumed to be always larger or equal to previously used filtration values or always smaller or equal than previous
   * values, ie. the changes are monotonous.
   * @return Number of the operation.
   */
  template <class BoundaryRange = std::initializer_list<face_key> >
  internal_key insert_face(face_key faceID,
                           const BoundaryRange& boundary,
                           dimension_type dimension,
                           filtration_value filtrationValue)
  {
    ++numArrow_;

    if (dimMax_ != -1 && dimension > dimMax_) {
      pers_.apply_identity();
      return numArrow_;
    }

    _store_filtration_value(filtrationValue);

    [[maybe_unused]] auto res = handleToKey_.try_emplace(faceID, numArrow_);

    GUDHI_CHECK(res.second, "Zigzag_persistence::insert_face - face already in the complex");

    // Compute the keys of the faces of the boundary.
    std::set<internal_key> translatedBoundary;  // set maintains the natural order on indices
    for (auto b : boundary) {
      translatedBoundary.insert(handleToKey_.at(b));  // TODO: add possibilities of coefficients
    }

    pers_.insert_face(translatedBoundary, dimension);

    return numArrow_;
  }

  /**
   * @brief Updates the zigzag persistence diagram after the removal of the given face.
   *
   * @param faceID ID representing the face to remove. Should be the same than the one used to insert it.
   * @param dimension Dimension of the face.
   * @param filtrationValue Filtration value associated to the removal.
   * Assumed to be always larger or equal to previously used filtration values or always smaller or equal than previous
   * values, ie. the changes are monotonous.
   * @return Number of the operation.
   */
  internal_key remove_face(face_key faceID, dimension_type dimension, filtration_value filtrationValue) {
    ++numArrow_;

    if (dimMax_ != -1 && dimension > dimMax_) {
      pers_.apply_identity();
      return numArrow_;
    }

    auto it = handleToKey_.find(faceID);
    GUDHI_CHECK(it != handleToKey_.end(), "Zigzag_persistence::remove_face - face not in the complex");

    _store_filtration_value(filtrationValue);

    pers_.remove_face(it->second, dimension);
    handleToKey_.erase(it);

    return numArrow_;
  }

  /**
   * @brief To use when a face is neither inserted nor removed, but the filtration moves along the identity operator
   * on homology level. Useful to keep the birth/death indices aligned when insertions/removals are purposely skipped
   * to avoid useless computation.
   * @return Number of the operation.
   */
  internal_key apply_identity() {
    ++numArrow_;
    pers_.apply_identity();
    return numArrow_;
  }

  /**
   * @brief Returns the "index persistence diagram" of the current filtration, that is, the pairs of atomic arrow
   * numbers corresponding to a birth-death pair. Does not contain points at infinity, only the cycle classes which
   * already died are represented.
   *
   * @return Reference to the list of intervals.
   */
  const std::vector<Index_interval>& get_index_persistence_diagram() const { return persistenceDiagram_; }

  /**
   * @brief Returns the filtration value \f$f(idx)\f$ associated to the index \f$idx\f$ returned
   * by @ref get_index_persistence_diagram.
   * 
   * @param idx Birth or death index
   * @return filtration_value Filtration value associated to @p idx.
   */
  filtration_value get_filtration_value_from_index(internal_key idx) {
    // lower_bound(x) returns leftmost y s.t. x <= y
    auto itBirth =
        std::lower_bound(filtrationValues_.begin(),
                         filtrationValues_.end(),
                         idx,
                         [](std::pair<internal_key, filtration_value> p, internal_key k) { return p.first < k; });
    if (itBirth == filtrationValues_.end() || itBirth->first > idx) {
      --itBirth;
    }
    return itBirth->second;
  };

  /**
   * @brief Returns the current persistence diagram.
   *
   * @param shortestInterval Threshold. Every bar shorter than the given value will be ignored. Default value: 0.
   * @param includeInfiniteBars If set to true, infinite bars are included in the diagram. Default value: true.
   * @return A vector of pairs of filtration values representing the persistence diagram.
   */
  std::vector<Filtration_value_interval> get_persistence_diagram(filtration_value shortestInterval = 0.,
                                                                 bool includeInfiniteBars = true) {
    std::vector<Filtration_value_interval> diag = _get_persistence_diagram(shortestInterval);

    if (includeInfiniteBars) {
      _retrieve_infinite_bars(diag);
    }

    return diag;
  }

 private:
  std::unordered_map<face_key, internal_key> handleToKey_;  /**< Map from input keys to internal keys. */
  dimension_type dimMax_;                                   /**< Maximal dimension of a bar to record. */
  std::vector<Index_interval> persistenceDiagram_;          /**< Stores current closed persistence intervals. */
  internal_key numArrow_;                                   /**< Current arrow number. */
  filtration_value previousFiltrationValue_;                /**< Filtration value of the previous arrow. */
  /**
   * @brief filtrationValues_ stores consecutive pairs (i,f) , (j,f') with f != f',
   * meaning that all inserted faces with key in [i;j-1] have filtration value f,
   * i is the smallest face index whose face has filtration value f.
   */
  std::vector<std::pair<internal_key, filtration_value> > filtrationValues_;
  Zigzag_persistence<FilteredZigzagOptions, false> pers_; /**< Class computing the pairs. */

  /**
   * @brief Stores the filtration value if the value is new. Assumes that the given value is either greater (or equal)
   * than previous ones or smaller (or equal) than previous ones.
   * 
   * @param filtrationValue Filtration value to store.
   */
  void _store_filtration_value(filtration_value filtrationValue) {
    if (filtrationValue != previousFiltrationValue_)  // check whether the filt value has changed
    {
      // consecutive pairs (i,f), (j,f') mean faces of index k in [i,j-1] have filtration value f
      previousFiltrationValue_ = filtrationValue;
      filtrationValues_.emplace_back(numArrow_, previousFiltrationValue_);
    }
  }

  /**
   * @brief Returns the current persistence diagram without infinite bars.
   *
   * @param shortestInterval Intervals shorter than the given value are ignored.
   * @return Vector of intervals.
   */
  std::vector<Filtration_value_interval> _get_persistence_diagram(filtration_value shortestInterval) {
    std::vector<Filtration_value_interval> diag;
    diag.reserve(persistenceDiagram_.size());

    // std::stable_sort(filtrationValues_.begin(), filtrationValues_.end(),
    //                  [](std::pair<internal_key, filtration_value> p1, std::pair<internal_key, filtration_value> p2) {
    //                    return p1.first < p2.first;
    //                  });

    for (auto bar : persistenceDiagram_) {
      filtration_value birth = get_filtration_value_from_index(bar.birth);
      filtration_value death = get_filtration_value_from_index(bar.death);
      if (birth > death) {
        std::swap(birth, death);
      }

      if (death - birth > shortestInterval) {
        diag.emplace_back(birth, death, bar.dim);
      }
    }

    return diag;
  }

  /**
   * @brief Computes the births of the current essential cycles.
   *
   * @param diag Reference to vector where to store the infinite bars.
   */
  void _retrieve_infinite_bars(std::vector<Filtration_value_interval>& diag) {
    auto stream_infinite_interval = [&](dimension_type dim, internal_key birthIndex) {
      if (dimMax_ == -1 || (dimMax_ != -1 && dim < dimMax_))
        diag.emplace_back(get_filtration_value_from_index(birthIndex), Filtration_value_interval::inf, dim);
    };

    pers_.get_current_infinite_intervals(stream_infinite_interval);
  }
};  // end class Filtered_zigzag_persistence_with_storage

/**
 * @ingroup zigzag_persistence
 *
 * @brief Class computing the zigzag persistent homology of a zigzag filtration. Algorithm based on \cite zigzag.
 * @details After construction of the class, the zigzag filtration should be given in a streaming like way, i.e.,
 * call @ref insert_face, @ref remove_face or @ref apply_identity for each step of the filtration in order of
 * the filtration. The bars of the diagram are retrieved via the given callback method every time
 * a pair with non-zero length is closed. To retrieve the open/infinite bars, use @ref get_current_infinite_intervals.
 * 
 * @tparam FilteredZigzagOptions Structure following the @ref FilteredZigzagOptions concept.
 * Default value: @ref Default_filtered_zigzag_options.
 */
template <class FilteredZigzagOptions = Default_filtered_zigzag_options>
class Filtered_zigzag_persistence {
 public:
  using Options = FilteredZigzagOptions;                        /**< Zigzag options. */
  using internal_key = typename Options::internal_key;          /**< Key and index type, has to be signed. */
  using face_key = typename Options::face_key;                  /**< Face ID type from external inputs. */
  using filtration_value = typename Options::filtration_value;  /**< Type for filtration values. */
  using dimension_type = typename Options::dimension_type;      /**< Type for dimension values. */

  /**
   * @brief Constructor.
   * @details After construction of the class, the zigzag filtration should be given in a streaming like way, i.e.,
   * call @ref insert_face, @ref remove_face or @ref apply_identity for each step of the filtration in order of
   * the filtration. The bars of the diagram are retrieved via the given callback method every time
   * a pair with non-zero length is closed. To retrieve the open/infinite bars, use @ref get_current_infinite_intervals.
   *
   * @param stream_interval Callback method to process the birth and death values of a persistence bar.
   * Has to take three arguments as input: first the dimension of the cycle, then the birth value of the cycle
   * and third the death value of the cycle. The values corresponds to the filtration values which were given at
   * insertions or removals.
   * @param preallocationSize Space for @p preallocationSize faces are reserved in the underlying structure.
   * Theoretically, any values works therefore, but for better performances, it is better to be as close as possible
   * to the maximal value of the number of faces stored at the same time. At a same time are stored faces which were
   * inserted before that time but not removed until that time. Default value: 0.
   * @tparam F Type of callback method.
   */
  template <typename F>
  Filtered_zigzag_persistence(F&& stream_interval, unsigned int preallocationSize = 0)
      : handleToKey_(preallocationSize),
        numArrow_(-1),
        keyToFiltrationValue_(preallocationSize),
        pers_(
            [&, stream_interval = std::forward<F>(stream_interval)](dimension_type dim,
                                                                    internal_key birth,
                                                                    internal_key death) {
              auto itB = keyToFiltrationValue_.find(birth);
              auto itD = keyToFiltrationValue_.find(death);
              if (itB->second != itD->second) stream_interval(dim, itB->second, itD->second);
              keyToFiltrationValue_.erase(itB);
              keyToFiltrationValue_.erase(itD);
            },
            preallocationSize) {}

  /**
   * @brief Updates the zigzag persistence diagram after the insertion of the given face.
   *
   * @tparam BoundaryRange Range type needing size, begin and end members.
   * @param faceID ID representing the inserted face.
   * @param boundary Boundary of the inserted face. The range should be composed of the IDs of all faces contained in
   * the boundary (i.e. with non-zero coefficients), using the ID specified as `faceID` when the corresponding face
   * was previously inserted (recall that the faces should be inserted in order of filtration).
   * @param dimension Dimension of the inserted face.
   * @param filtrationValue Filtration value associated to the face.
   * Assumed to be always larger or equal to previously used filtration values or always smaller or equal than previous
   * values, ie. the changes are monotonous.
   */
  template <class BoundaryRange = std::initializer_list<face_key> >
  void insert_face(face_key faceID,
                   const BoundaryRange& boundary,
                   dimension_type dimension,
                   filtration_value filtrationValue) {
    ++numArrow_;

    [[maybe_unused]] auto res = handleToKey_.try_emplace(faceID, numArrow_);

    GUDHI_CHECK(res.second, "Zigzag_persistence::insert_face - face already in the complex");

    keyToFiltrationValue_.try_emplace(numArrow_, filtrationValue);

    // Compute the keys of the faces of the boundary.
    std::set<internal_key> translatedBoundary;  // set maintains the natural order on indices
    for (auto b : boundary) {
      translatedBoundary.insert(handleToKey_.at(b));  // TODO: add possibilities of coefficients
    }

    pers_.insert_face(translatedBoundary, dimension);
  }

  /**
   * @brief Updates the zigzag persistence diagram after the removal of the given face.
   *
   * @param faceID ID representing the face to remove. Should be the same than the one used to insert it.
   * @param dimension Dimension of the face.
   * @param filtrationValue Filtration value associated to the removal.
   * Assumed to be always larger or equal to previously used filtration values or always smaller or equal than previous
   * values, ie. the changes are monotonous.
   */
  void remove_face(face_key faceID, dimension_type dimension, filtration_value filtrationValue) {
    ++numArrow_;

    auto it = handleToKey_.find(faceID);
    GUDHI_CHECK(it != handleToKey_.end(), "Zigzag_persistence::remove_face - face not in the complex");

    keyToFiltrationValue_.try_emplace(numArrow_, filtrationValue);

    pers_.remove_face(it->second, dimension);
    handleToKey_.erase(it);
  }

  /**
   * @brief To use when a face is neither inserted nor removed, but the filtration moves along the identity operator
   * on homology level. Useful to keep the birth/death indices aligned when insertions/removals are purposely skipped
   * to avoid useless computation.
   */
  void apply_identity() {
    ++numArrow_;
    pers_.apply_identity();
  }

  /**
   * @brief Outputs through the given callback method all current infinite bars.
   * 
   * @tparam F Type of the callback method. Takes two arguments: the dimension of the cycle and the birth value
   * of the cycle.
   * @param stream_infinite_interval Method processing the unpaired birth values.
   */
  template <typename F>
  void get_current_infinite_intervals(F&& stream_infinite_interval) {
    pers_.get_current_infinite_intervals(
        [&](dimension_type dim, internal_key birth) { stream_infinite_interval(dim, keyToFiltrationValue_.at(birth)); });
  }

 private:
  template <typename key_type, typename value_type>
  using dictionary = std::unordered_map<key_type, value_type>;  // TODO: benchmark with other map types

  dictionary<face_key, internal_key> handleToKey_;                  /**< Map from input keys to internal keys. */
  internal_key numArrow_;                                           /**< Current arrow number. */
  dictionary<internal_key, filtration_value> keyToFiltrationValue_; /**< Face Key to filtration value map. */
  Zigzag_persistence<FilteredZigzagOptions, true> pers_;            /**< Class computing the pairs. */
};  // end class Filtered_zigzag_persistence

}  // namespace zigzag_persistence
}  // namespace Gudhi

#endif  // FILTERED_ZIGZAG_PERSISTENCE_H_
