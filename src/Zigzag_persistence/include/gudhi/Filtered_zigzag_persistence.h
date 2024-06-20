/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2021 Inria
 *
 *    Modification(s):
 *      - 2023/05 Hannah Schreiber: Rework of the interface, reorganization and debug
 *      - 2023/05 Hannah Schreiber: Addition of infinit bars
 *      - 2024/06 Hannah Schreiber: Separation of the zigzag algorithm from the filtration value management
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file Filtered_zigzag_persistence.h
 * @author Clément Maria, Hannah Schreiber
 * @brief Contains the implementation of the @ref Interval structure and the
 * @ref Gudhi::zigzag_persistence::Filtered_zigzag_persistence_with_storage and
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
#include <gudhi/Zigzag_persistence.h>

namespace Gudhi {
namespace zigzag_persistence {

/**
 * @ingroup zigzag_persistence
 *
 * @brief Structure to store persistence intervals by their birth and death values.
 * 
 * @tparam value_type Type for the birth and death indices.
 */
template <typename value_type>
struct Interval {
  /**
   * @brief Default constructor. All values are initialized with default values.
   */
  Interval() {}
  /**
   * @brief Constructor.
   * 
   * @param dim Dimension of the cycle.
   * @param b Birth index or value of the cycle.
   * @param d Death index or value of the cycle.
   */
  Interval(int dim, value_type b, value_type d) : dim_(dim), b_(b), d_(d) {}
  /**
   * @brief Returns the dimension of the homological feature corresponding to the interval.
   * 
   * @return Stored dimension.
   */
  int dim() const { return dim_; }
  /**
   * @brief Returns the birth value of the interval.
   * 
   * @return The stored birth.
   */
  value_type birth() const { return b_; }
  /**
   * @brief Returns the death value of the interval.
   * 
   * @return The stored death.
   */
  value_type death() const { return d_; }

 protected:
  int dim_;       /**< Homological dimension. */
  value_type b_;  /**< Value associated to the interval birth. */
  value_type d_;  /**< Value associated to the interval death. */
};

/**
 * @ingroup zigzag_persistence
 *
 * @brief Default options for @ref Filtered_zigzag_persistence_with_storage and @ref Filtered_zigzag_persistence.
 */
struct Default_filtered_zigzag_options {
  using internal_key = int;        /**< Face ID used internaly, must be signed. */
  using face_key = int;            /**< Face ID used in the given boundaries. */
  using filtration_value = double; /**< Filtration value type. */
  using dimension_type = int;      /**< Dimension value type. */
  /**
   * @brief Column type use by the internal matrix.
   */
  static const Gudhi::persistence_matrix::Column_types column_type =
      Gudhi::persistence_matrix::Column_types::INTRUSIVE_LIST;
};

/**
 * @ingroup zigzag_persistence
 *
 * @brief Class computating the zigzag persistent homology of a zigzag filtration. Algorithm based on \cite zigzag.
 * Eventhough the insertions and removals are given in a "stream-like" way, the barcode and other values are
 * stored during the whole process and not removed. It is therefore suited for smaller filtrations where the clean
 * ups produce a higher overhead than the memory consumption.
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
  using Index_interval = Interval<internal_key>;                /**< Persistence interval type. */

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
      return Base::d_ - Base::b_;
    }
    /**
     * @brief Returns the absolute length of the log values of birth and death, i.e.  \f$|\log d - \log b|\f$.
     */
    filtration_value log_length() const {
      if (Base::b_ == Base::d_) {
        return 0;
      }  // otherwise inf - inf would return nan.
      return std::log2(static_cast<double>(Base::d_)) - std::log2(static_cast<double>(Base::b_));
    }
  };

  /**
   * @brief Constructor.
   * @details After construction of the class, the zigzag filtration should be given in a streaming like way, i.e.,
   * call @ref insert_face, @ref remove_face or @ref apply_identity for each step of the filtration in order of
   * the filtration. To retrieve the current persistence diagram at any moment of the filtration,
   * use @ref get_persistence_diagram or @ref get_index_persistence_diagram.
   *
   * @param minNumberOfFaces Minimum number of faces that will be in a complex at some point in the filtration.
   * If the maximal number of faces is known in advance, the memory allocation can be better optimized.
   * Default value: 0.
   * @param ignoreCyclesAboveDim Ignores cycles in dimension larger or equal in the final diagram.
   * If -1, no cycles are ignored. Default value: -1.
   */
  Filtered_zigzag_persistence_with_storage(unsigned int minNumberOfFaces = 0, int ignoreCyclesAboveDim = -1)
      : dimMax_(ignoreCyclesAboveDim),
        persistenceDiagram_(),
        numArrow_(-1),
        previousFiltrationValue_(std::numeric_limits<filtration_value>::infinity()),
        pers_(
            [&](dimension_type dim, internal_key birth, internal_key death) {
              if (dimMax_ == -1 || (dimMax_ != -1 && dim < dimMax_)) {  // don't record intervals over max dim
                persistenceDiagram_.emplace_back(dim, birth, death);
              }
            },
            minNumberOfFaces) {}

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
    ++numArrow_;

    if (dimMax_ != -1 && dimension > dimMax_) {
      pers_.apply_identity();
      return;
    }

    if (filtrationValue != previousFiltrationValue_)  // check whether the filt value has changed
    {  // consecutive pairs (i,f), (j,f') mean faces of index k in [i,j-1] have
      previousFiltrationValue_ = filtrationValue;  // filtration value f
      filtrationValues_.emplace_back(numArrow_, previousFiltrationValue_);
    }

    [[maybe_unused]] auto res = handleToKey_.try_emplace(faceID, numArrow_);

    GUDHI_CHECK(res.second, "Zigzag_persistence::insert_face - face already in the complex");

    // Reduce the boundary of zzsh in the basis of cycles.
    // Compute the keys of the faces of the boundary of zzsh.
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
   * Assumed to be larger or equal to previously used filtration values.
   */
  void remove_face(face_key faceID, dimension_type dimension, filtration_value filtrationValue) {
    ++numArrow_;

    if (dimMax_ != -1 && dimension > dimMax_) {
      pers_.apply_identity();
      return;
    }

    auto it = handleToKey_.find(faceID);
    GUDHI_CHECK(it != handleToKey_.end(), "Zigzag_persistence::remove_face - face not in the complex");

    if (filtrationValue != previousFiltrationValue_)  // check whether the filt value has changed
    {  // consecutive pairs (i,f), (j,f') mean faces of index k in [i,j-1] have
      previousFiltrationValue_ = filtrationValue;  // filtration value f
      filtrationValues_.emplace_back(numArrow_, previousFiltrationValue_);
    }

    pers_.remove_face(it->second, dimension);
    handleToKey_.erase(it);
  }

  /**
   * @brief To use when a face is neither inserted nor removed, but the filtration moves along the identity operator
   * on homology level. Useful to keep the birth/death indices aligned when insertions/removals are purposly skipped
   * to avoid useless computation.
   */
  void apply_identity() {
    ++numArrow_;
    pers_.apply_identity();
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
   * @brief Returns the current persistence diagram.
   *
   * @param shortestInterval Threshold. Every bar shorter than the given value will be ignored. Default value: 0.
   * @param includeInfinitBars If set to true, infinit bars are included in the diagram. Default value: false.
   * @return A vector of pairs of filtration values representing the persistence diagram.
   */
  std::vector<Filtration_value_interval> get_persistence_diagram(filtration_value shortestInterval = 0.,
                                                                 bool includeInfinitBars = false) {
    std::vector<Filtration_value_interval> diag = _get_persistence_diagram(shortestInterval);

    if (includeInfinitBars) {
      _retrieve_infinit_bars(diag);
    }

    return diag;
  }

 private:
  std::unordered_map<face_key, internal_key> handleToKey_; /**< Map from input keys to internal keys. */
  dimension_type dimMax_;                                  /**< Maximal dimension of a bar to record. */
  std::vector<Index_interval> persistenceDiagram_;         /**< Stores current closed persistence intervals. */
  internal_key numArrow_;                                  /**< Current arrow number. */
  filtration_value previousFiltrationValue_;               /**< Filtration value of the previous arrow. */
  /**
   * @brief filtrationValues_ stores consecutive pairs (i,f) , (j,f') with f != f',
   * meaning that all inserted faces with key in [i;j-1] have filtration value f,
   * i is the smallest face index whose face has filtration value f.
   */
  std::vector<std::pair<internal_key, filtration_value> > filtrationValues_;
  Zigzag_persistence<FilteredZigzagOptions, false> pers_; /**< Class computing the pairs. */

  /**
   * @brief Returns the current persistence diagram without infinit bars.
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
  void _retrieve_infinit_bars(std::vector<Filtration_value_interval>& diag) {
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

    auto stream_infinit_interval = [&](dimension_type dim, internal_key birthIndex) {
      if (dimMax_ == -1 || (dimMax_ != -1 && dim < dimMax_))
        diag.emplace_back(dim, birth(birthIndex), std::numeric_limits<filtration_value>::infinity());
    };

    pers_.get_current_infinit_intervals(stream_infinit_interval);
  }
};  // end class Filtered_zigzag_persistence_with_storage

/**
 * @ingroup zigzag_persistence
 *
 * @brief Class computating the zigzag persistent homology of a zigzag filtration. Algorithm based on \cite zigzag.
 * 
 * @tparam FilteredZigzagOptions Structure following the @ref FilteredZigzagOptions concept.
 * Default value: @ref Default_filtered_zigzag_options.
 */
template <class FilteredZigzagOptions = Default_filtered_zigzag_options>
class Filtered_zigzag_persistence {
 public:
  using Options = FilteredZigzagOptions;                       /**< Zigzag options. */
  using internal_key = typename Options::internal_key;         /**< Key and index type, has to be signed. */
  using face_key = typename Options::face_key;                 /**< Face ID type from external inputs. */
  using filtration_value = typename Options::filtration_value; /**< Type for filtration values. */
  using dimension_type = typename Options::dimension_type;     /**< Type for dimension values. */

  /**
   * @brief Constructor.
   * @details After construction of the class, the zigzag filtration should be given in a streaming like way, i.e.,
   * call @ref insert_face, @ref remove_face or @ref apply_identity for each step of the filtration in order of
   * the filtration. The bars of the diagram are retrieved via the given callback method every time
   * a pair with non-zero length is closed. To retrieve the open/infinit bars, use @ref get_current_infinit_intervals.
   *
   * @param stream_interval Callback method to process the birth and death values of a persistence bar.
   * Has to take three arguments as input: first the dimension of the cycle, then the birth value of the cycle
   * and third the death value of the cycle. The values corresponds to the filtration values which were given at
   * insertions or removals.
   * @param minNumberOfFaces Minimum number of faces that will be in a complex at some point in the filtration.
   * If the maximal number of faces is known in advance, the memory allocation can be better optimized.
   * Default value: 0.
   */
  Filtered_zigzag_persistence(std::function<void(dimension_type, filtration_value, filtration_value)> stream_interval,
                              unsigned int minNumberOfFaces = 0)
      : handleToKey_(minNumberOfFaces),
        numArrow_(-1),
        keyToFiltrationValue_(minNumberOfFaces),
        stream_interval_(std::move(stream_interval)),
        pers_(
            [&](dimension_type dim, internal_key birth, internal_key death) {
              auto itB = keyToFiltrationValue_.find(birth);
              auto itD = keyToFiltrationValue_.find(death);
              if (itB->second != itD->second) stream_interval_(dim, itB->second, itD->second);
              keyToFiltrationValue_.erase(itB);
              keyToFiltrationValue_.erase(itD);
            },
            minNumberOfFaces) {}

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
    ++numArrow_;

    [[maybe_unused]] auto res = handleToKey_.try_emplace(faceID, numArrow_);

    GUDHI_CHECK(res.second, "Zigzag_persistence::insert_face - face already in the complex");

    keyToFiltrationValue_.try_emplace(numArrow_, filtrationValue);

    // Reduce the boundary of zzsh in the basis of cycles.
    // Compute the keys of the faces of the boundary of zzsh.
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
   * Assumed to be larger or equal to previously used filtration values.
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
   * on homology level. Useful to keep the birth/death indices aligned when insertions/removals are purposly skipped
   * to avoid useless computation.
   */
  void apply_identity() {
    ++numArrow_;
    pers_.apply_identity();
  }

  /**
   * @brief Outputs through the given callback method all current infinit bars.
   * 
   * @tparam F Type of the callback method. Takes two arguments: the dimension of the cycle and the birth value
   * of the cycle.
   * @param stream_infinit_interval Method processing the unpaired birth values.
   */
  template <typename F>
  void get_current_infinit_intervals(F&& stream_infinit_interval) {
    pers_.get_current_infinit_intervals(
        [&](dimension_type dim, internal_key birth) { stream_infinit_interval(dim, keyToFiltrationValue_.at(birth)); });
  }

 private:
  template <typename key_type, typename value_type>
  using dictionnary = std::unordered_map<key_type, value_type>;  // TODO: benchmark with other map types

  dictionnary<face_key, internal_key> handleToKey_;                   /**< Map from input keys to internal keys. */
  internal_key numArrow_;                                             /**< Current arrow number. */
  dictionnary<internal_key, filtration_value> keyToFiltrationValue_;  /**< Face Key to filtration value map. */
  std::function<void(int,filtration_value,filtration_value)> stream_interval_;  /**< Callback method for finite bars. */
  Zigzag_persistence<FilteredZigzagOptions, true> pers_;              /**< Class computing the pairs. */
};  // end class Filtered_zigzag_persistence

}  // namespace zigzag_persistence
}  // namespace Gudhi

#endif  // FILTERED_ZIGZAG_PERSISTENCE_H_
