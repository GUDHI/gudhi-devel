/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria and Hannah Schreiber
 *
 *    Copyright (C) 2023-25 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file oscillating_rips_edge_ranges.h
 * @author Clément Maria, Hannah Schreiber
 * @brief Contains the implementation of the @ref Gudhi::zigzag_persistence::Oscillating_rips_vertex_order_policy enum,
 * @ref Gudhi::zigzag_persistence::Oscillating_rips_edge_iterator_base class,
 * @ref Gudhi::zigzag_persistence::Oscillating_rips_edge_iterator_range class and
 * @ref Gudhi::zigzag_persistence::Oscillating_rips_edge_vector_range_constructor class.
 */

#ifndef ZIGZAG_OSCILLATING_RIPS_EDGE_RANGES_H_
#define ZIGZAG_OSCILLATING_RIPS_EDGE_RANGES_H_

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <stdexcept>

#include <boost/iterator/iterator_facade.hpp>

#ifdef GUDHI_USE_TBB
#include <tbb/tbb.h>
#endif

#include <gudhi/choose_n_farthest_points.h>
#include <gudhi/pick_n_random_points.h>
#include <gudhi/Debug_utils.h>
#include <gudhi/Zigzag_persistence/Zigzag_edge.h>
#include <gudhi/Zigzag_persistence/edge_modifiers.h>

namespace Gudhi {
namespace zigzag_persistence {

/**
 * @brief Order policy for the points.
 *
 * @ingroup zigzag_persistence
 */
enum Oscillating_rips_vertex_order_policy {
  ALREADY_ORDERED,         /**< The given range of points is considered ordered. */
  FARTHEST_POINT_ORDERING, /**< The points are reordered using @ref Gudhi::subsampling::choose_n_farthest_points.*/
  RANDOM_POINT_ORDERING    /**< The points are shuffled randomly. */
};

/**
 * @private
 * @ingroup zigzag_persistence
 * @brief Initialize epsilon values and distance matrix from the given parameters and point cloud.
 *
 * @tparam Filtration_value Filtration value type.
 */
template <typename Filtration_value>
class Oscillating_rips_initializer
{
 public:
  /**
   * @brief Initializes the distance matrix and epsilon values.
   *
   * @tparam PointRange Point range type.
   * @tparam DistanceFunction Type of the distance function.
   * @param epsilonValues Container for the epsilon values.
   * @param distanceMatrix Container for the distance matrix.
   * @param points Point cloud as a range.The format of a point has to correspond to the input format of the
   * distance function.
   * @param distance Distance function. Has to take two points as it from the range @p points as input parameters
   * and return the distance between those points.
   * @param orderPolicy Order policy for the points. Can be either
   * @ref Oscillating_rips_vertex_order_policy::FARTHEST_POINT_ORDERING,
   * @ref Oscillating_rips_vertex_order_policy::ALREADY_ORDERED or
   * @ref Oscillating_rips_vertex_order_policy::RANDOM_POINT_ORDERING.
   */
  template <typename PointRange, typename DistanceFunction>
  static void initialize(std::vector<Filtration_value>& epsilonValues,
                         std::vector<std::vector<std::pair<int, Filtration_value> > >& distanceMatrix,
                         const PointRange& points,
                         DistanceFunction&& distance,
                         Oscillating_rips_vertex_order_policy orderPolicy)
  {
    std::size_t n = points.size();  // number of points
    PointRange sortedPoints;
    sortedPoints.reserve(n);

    // compute epsilon values
    if (orderPolicy == Oscillating_rips_vertex_order_policy::ALREADY_ORDERED) {
      sortedPoints.assign(points.begin(), points.end());
      epsilonValues = _compute_epsilon_values(sortedPoints, distance);
    } else if (orderPolicy == Oscillating_rips_vertex_order_policy::FARTHEST_POINT_ORDERING) {
      epsilonValues.reserve(n);
      Gudhi::subsampling::choose_n_farthest_points(distance,
                                                   points,
                                                   n,  // final size
                                                   0,  // starting point
                                                   std::back_inserter(sortedPoints),
                                                   std::back_inserter(epsilonValues));
      // need to shift values output by subsampling:
      for (unsigned int i = 1; i < n; ++i) {
        epsilonValues[i - 1] = epsilonValues[i];
      }
      epsilonValues[n - 1] = 0;
    } else {
      Gudhi::subsampling::pick_n_random_points(points, n, std::back_inserter(sortedPoints));
      epsilonValues = _compute_epsilon_values(sortedPoints, distance);
    }

    // compute the distance matrix
    distanceMatrix = _compute_distance_matrix(sortedPoints, distance);
  }

  /**
   * @brief Initializes the distance matrix.
   * 
   * @tparam PointRange Point range type.
   * @tparam DistanceFunction Type of the distance function.
   * @param distanceMatrix Container for the distance matrix. The order will correspond to the order of the given
   * points.
   * @param sortedPoints Point cloud as a range.The format of a point has to correspond to the input format of the
   * distance function. The order of the points will not change.
   * @param distance Distance function. Has to take two points as it from the range @p points as input parameters
   * and return the distance between those points.
   */
  template <typename PointRange, typename DistanceFunction>
  static void initialize(std::vector<std::vector<std::pair<int, Filtration_value> > >& distanceMatrix,
                         const PointRange& sortedPoints,
                         DistanceFunction&& distance)
  {
    distanceMatrix = _compute_distance_matrix(sortedPoints, distance);
  }

  /**
   * @brief The two input types std::pair<int, Filtration_value> encode pairs
   * @f$(j, d(p_i,p_j))@f$ and @f$(k, d(p_i,p_k))@f$ for some fixed point @f$p_i@f$.
   * The operator() orders edges by length. By convention, if lengths are equal,
   * it orders pairs by taking the smaller vertex label between @f$j@f$ and @f$k@f$.
   */
  struct Point_distance_comp {
    bool operator()(const std::pair<int, Filtration_value>& p1, const std::pair<int, Filtration_value>& p2) const
    {
      {
        if (p1.second != p2.second) {
          return p1.second < p2.second;  // shorter first
        }
        return p1.first < p2.first;
      }
    }
  };

 private:
  /**
   * @brief Compute the epsilon values for an ordered set of points, measuring the
   * sparsity of the ordering.
   *
   * @details Let \f$P = \{p_0, \ldots, p_{n-1}\f$ be the ordered set of points. Then
   * the method sets <CODE>eps_range[i]<\CODE> with the value \$f\varepsilon_i\$f,
   * defined as \f$\varepsilon_i = d_H(P_i,P)\f$, the Hausdorff between the points
   * \f$P_i= \{p_0, \ldots, p_{i}\}\f$ and the entire point cloud
   * \f$P = \{p_0, \ldots, p_{n-1}\}\f$.
   * 
   * @tparam PointRange Point range type.
   * @tparam DistanceFunction Type of the distance function.
   * @param sortedPoints Point cloud as an ordered range. The format of a point has to correspond to the input
   * format of the distance function.
   * @param distance Distance function. Has to take two points as it from the range @p points as input parameters
   * and return the distance between those points.
   * @return Vector of decreasing epsilon values ending with 0.
   */
  template <typename PointRange, typename DistanceFunction>
  static std::vector<Filtration_value> _compute_epsilon_values(const PointRange& sortedPoints,
                                                               DistanceFunction&& distance)
  {
    std::size_t n = sortedPoints.size();
    std::vector<Filtration_value> eps_range(n, std::numeric_limits<double>::infinity());

    // compute all \f$\varepsilon_i\f$ values, such that eps_range[i] ==
    // eps_i==d_H(P_i,P), for i=0 ... n-1:
    for (std::size_t i = 0; i < n; ++i) {
      // entering step i, maintain eps_range[j] = eps_j for j<i, and
      // eps_range[k] = d(p_k, P_{i-1}) for k >= i.
#ifdef GUDHI_USE_TBB
      tbb::parallel_for(std::size_t(i + 1), n, [&](std::size_t k) {
        // set eps_range[k] <- d(p_k, P_i) ==
        //                            min{ d(p_k, P_{i-1}), d(p_k, p_i) }  for k >= i.
        double dist_i_k = distance(sortedPoints[i], sortedPoints[k]);
        if (dist_i_k < eps_range[k]) {
          eps_range[k] = dist_i_k;
        }
      });
#else
      for (std::size_t k = i + 1; k < n; ++k) {
        // set eps_range[k] <- d(p_k, P_i) ==
        //                            min{ d(p_k, P_{i-1}), d(p_k, p_i) }  for k >= i.
        double dist_i_k = distance(sortedPoints[i], sortedPoints[k]);
        if (dist_i_k < eps_range[k]) {
          eps_range[k] = dist_i_k;
        }
      }
#endif
      // we have now eps_range[k] = d(p_k, P_i) for k > i.
      // to do: implement parallel version by dividing the vector
      // set eps_range[i] <- eps_i = d_H(P_i,P) = max_{k>i} d(p_k, P_i)
      double eps_i = 0.;
      for (std::size_t k = i + 1; k < n; ++k) {
        if (eps_range[k] > eps_i) {
          eps_i = eps_range[k];
        }
      }
      eps_range[i] = eps_i;
    }

    return eps_range;
  }

  /**
   * @brief Returns the sparse distance matrix. A cell is represented as a pair of its row index and its value,
   * i.e., @f$(j, d(p_i,p_j))@f$ is stored in column @f$i@f$ to indicate that the distance of the @f$i^{th}@f$
   * point has distance @f$d(p_i,p_j)@f$ from the @f$j^{th}@f$ point.
   *
   * @tparam PointRange Point range type.
   * @tparam DistanceFunction Type of the distance function.
   * @param sortedPoints Point cloud as an ordered range. The format of a point has to correspond to the input
   * format of the distance function.
   * @param distance Distance function. Has to take two points as it from the range @p points as input parameters
   * and return the distance between those points.
   *
   * @return A vector of vectors of pairs of row indices and distances.
   */
  template <typename PointRange, typename DistanceFunction>
  static std::vector<std::vector<std::pair<int, Filtration_value> > > _compute_distance_matrix(
      const PointRange& sortedPoints,
      DistanceFunction&& distance)
  {
    std::vector<std::vector<std::pair<int, Filtration_value> > > distanceMatrix(sortedPoints.size());
#ifdef GUDHI_USE_TBB
    tbb::parallel_for(std::size_t(0), sortedPoints.size(), [&](std::size_t i) {
      distanceMatrix[i].resize(i);
      for (std::size_t j = 0; j < i; ++j) {
        distanceMatrix[i][j] = std::make_pair(j, distance(sortedPoints[i], sortedPoints[j]));
      }
      // distanceMatrix[i] is sorted by (j, d(p_i,p_j)) < (k, d(p_i,p_k)) iff
      // d(p_i,p_j) < d(p_i,p_k) or (j<k in case d(p_i,p_j) == d(p_i,p_k)).
      std::stable_sort(distanceMatrix[i].begin(), distanceMatrix[i].end(), Point_distance_comp());
    });
#else
    for (std::size_t i = 0; i < sortedPoints.size(); ++i) {  // for all vertices
      distanceMatrix[i].resize(i);
      for (std::size_t j = 0; j < i; ++j) {
        distanceMatrix[i][j] = std::make_pair(j, distance(sortedPoints[i], sortedPoints[j]));
      }
      std::stable_sort(distanceMatrix[i].begin(), distanceMatrix[i].end(), Point_distance_comp());
    }
#endif

    return distanceMatrix;
  }
};

/**
 * @class Oscillating_rips_edge_iterator_base oscillating_rips_edge_ranges.h \
 * gudhi/Zigzag_persistence/oscillating_rips_edge_ranges.h
 * @ingroup zigzag_persistence
 *
 * @brief Heavy base for a custom iterator over the edges of an oscillating rips filtration.
 * 
 * @tparam Filtration_value Filtration value type. Should be compatible with the edge modifier.
 * @tparam EdgeFiltrationTransformer Modifier for the edge filtration values. If no modifications are wanted,
 * use @ref Identity_edge_modifier. Default value: @ref Identity_edge_modifier.
 */
template <typename Filtration_value, class EdgeFiltrationTransformer = Identity_edge_modifier>
class Oscillating_rips_edge_iterator_base
{
 public:
  /**
   * @brief Construct.
   * 
   * @param nu Lower multiplier.
   * @param mu Upper multiplier.
   * @param epsilonValues Pointer to the epsilon values.
   * @param distanceMatrix Pointer to the distance matrix. The distance matrix is lower triangular such that
   * `distanceMatrix[i][j]` is only valid if `i < j`.
   */
  Oscillating_rips_edge_iterator_base(
      Filtration_value nu,
      Filtration_value mu,
      const std::vector<Filtration_value>* epsilonValues,
      const std::vector<std::vector<std::pair<int, Filtration_value> > >* distanceMatrix)
      : epsilonValues_(epsilonValues),
        distanceMatrix_(distanceMatrix),
        nu_(EdgeFiltrationTransformer::apply_inverse_modifier(nu)),
        mu_(EdgeFiltrationTransformer::apply_inverse_modifier(mu)),
        currentEdge_(0, 0, std::numeric_limits<Filtration_value>::infinity(), true),
        epsilonIndex_(0),
        rowIndex_(1),
        inPositiveDirection_(true),
        insertVertex_(true)
  {
    GUDHI_CHECK(epsilonValues->size() == distanceMatrix->size(),
                "Epsilon values and distance matrix are not compatible.");

    if (epsilonValues->empty()) {
      _set_end();
      return;
    }

    const auto& row = (*distanceMatrix_)[1];
    auto it = std::upper_bound(
        row.begin(),
        row.end(),
        std::pair<int, Filtration_value>(distanceMatrix_->size(), mu_ * (*epsilonValues_)[epsilonIndex_]),
        typename Oscillating_rips_initializer<Filtration_value>::Point_distance_comp());
    columnIndex_ = it - row.begin();
  }

  /**
   * @brief Default constructor. Corresponds to the end iterator.
   */
  Oscillating_rips_edge_iterator_base()
      : epsilonValues_(nullptr),
        distanceMatrix_(nullptr),
        nu_(0),
        mu_(0),
        currentEdge_(0, 0, 0, true),
        epsilonIndex_(0),
        rowIndex_(0),
        columnIndex_(0),
        inPositiveDirection_(true),
        insertVertex_(true)
  {}

  /**
   * @brief Indicates if two iterators are equal.
   */
  bool equal(Oscillating_rips_edge_iterator_base const& other) const
  {
    return rowIndex_ == other.rowIndex_ && currentEdge_ == other.currentEdge_;
  }

  /**
   * @brief Returns the value of the dereferenced iterator.
   *
   * @return Current @ref Zigzag_edge edge.
   */
  const Zigzag_edge<Filtration_value>& dereference() const { return currentEdge_; }

  /**
   * @brief Increments the iterator.
   */
  void increment()
  {
    auto size = distanceMatrix_->size();
    if (epsilonIndex_ < size - 1) {
      if (insertVertex_) {
        _update_edge_as_positive_vertex();
        insertVertex_ = false;
        return;
      }

      if (inPositiveDirection_ && _positive_col_index_is_not_valid()) {
        while (_positive_col_index_is_not_valid() && rowIndex_ <= epsilonIndex_) {
          ++rowIndex_;
          _initialize_positive_col_index();
        }
        if (_positive_col_index_is_not_valid()) {
          _initialize_negative_col_index();
          inPositiveDirection_ = false;
        }
      }

      if (!inPositiveDirection_ && _negative_col_index_is_not_valid()) {
        while (_negative_col_index_is_not_valid() && rowIndex_ > 1) {
          --rowIndex_;
          _initialize_negative_col_index();
        }
        if (_negative_col_index_is_not_valid()) {
          _initialize_positive_col_index();
          ++epsilonIndex_;
          inPositiveDirection_ = true;
          if (epsilonIndex_ == size - 1) {
            //   _set_end();
            rowIndex_ = size;
            _update_edge_as_negative_vertex();
            return;
          }
          _update_edge_as_positive_vertex();
          return;
        }
      }

      if (inPositiveDirection_) {
        --columnIndex_;
        _update_edge(epsilonIndex_, true);
        while (_positive_col_index_is_not_valid() && rowIndex_ <= epsilonIndex_) {
          ++rowIndex_;
          _initialize_positive_col_index();
        }
        if (_positive_col_index_is_not_valid()) {
          ++rowIndex_;
        }
        if (rowIndex_ == epsilonIndex_ + 2) {
          inPositiveDirection_ = false;
          --rowIndex_;
          _initialize_negative_col_index();
        }
        return;
      }

      _update_edge(epsilonIndex_, false);
      ++columnIndex_;
      while (_negative_col_index_is_not_valid() && rowIndex_ > 1) {
        --rowIndex_;
        _initialize_negative_col_index();
      }
      if (_negative_col_index_is_not_valid()) {
        --rowIndex_;
      }
      if (rowIndex_ == 0) {
        ++epsilonIndex_;
        if (epsilonIndex_ == size - 1) {
          rowIndex_ = size + 1;
          return;
        }
        insertVertex_ = true;
        inPositiveDirection_ = true;
        ++rowIndex_;
        _initialize_positive_col_index();
      }
      return;
    }
    if (rowIndex_ > 1) {
      --rowIndex_;
      _update_edge_as_negative_vertex();
      return;
    }

    _set_end();
  }

 private:
  const std::vector<Filtration_value>* epsilonValues_;                                 /**< Epsilon values. */
  /**
   * @brief Lower triangular distance matrix (`distanceMatrix[i][j]` is only valid if `i < j`).
   */
  const std::vector<std::vector<std::pair<int, Filtration_value> > >* distanceMatrix_;
  const Filtration_value nu_;                                                          /**< Lower multiplier. */
  const Filtration_value mu_;                                                          /**< Upper multiplier. */
  Zigzag_edge<Filtration_value> currentEdge_;         /**< Stores the current edge in the range. */
  std::size_t epsilonIndex_, rowIndex_, columnIndex_; /**< Indices indicating the next position in the range. */
  bool inPositiveDirection_, insertVertex_;           /**< Next direction and indicates if next "edge" is a vertex. */

  /**
   * @brief Set the iterator as the end iterator.
   */
  void _set_end()
  {
    rowIndex_ = 0;
    currentEdge_ = Zigzag_edge<Filtration_value>();
  }

  /**
   * @brief Initialize the column index for a new sequence of positive edges.
   */
  void _initialize_positive_col_index()
  {
    const auto& row = (*distanceMatrix_)[rowIndex_];
    auto it = std::upper_bound(
        row.begin(),
        row.end(),
        std::pair<int, Filtration_value>(distanceMatrix_->size(), mu_ * (*epsilonValues_)[epsilonIndex_]),
        typename Oscillating_rips_initializer<Filtration_value>::Point_distance_comp());
    columnIndex_ = it - row.begin();
  }

  /**
   * @brief Initialize the column index for a new sequence of negative edges.
   */
  void _initialize_negative_col_index()
  {
    const auto& row = (*distanceMatrix_)[rowIndex_];
    const auto eps = (*epsilonValues_)[epsilonIndex_ + 1];
    auto it = std::lower_bound(row.begin(),
                               row.end(),
                               std::pair<int, Filtration_value>(0, nu_ * eps),
                               typename Oscillating_rips_initializer<Filtration_value>::Point_distance_comp());
    while (it != row.end() && it->second == nu_ * eps) ++it;
    columnIndex_ = it - row.begin();
  }

  /**
   * @brief Indicates if the column index needs to be initialized or not.
   */
  bool _positive_col_index_is_not_valid()
  {
    return columnIndex_ == 0 ||
           (rowIndex_ != (epsilonIndex_ + 1) &&
            (*distanceMatrix_)[rowIndex_][columnIndex_ - 1].second <= nu_ * (*epsilonValues_)[epsilonIndex_]);
  }

  /**
   * @brief Indicates if the column index needs to be initialized or not.
   */
  bool _negative_col_index_is_not_valid()
  {
    const auto& row = (*distanceMatrix_)[rowIndex_];
    return columnIndex_ == row.size() || row[columnIndex_].second > mu_ * (*epsilonValues_)[epsilonIndex_];
  }

  /**
   * @brief Updates the current edge.
   *
   * @param i Epsilon value index.
   * @param direction Direction.
   */
  void _update_edge(std::size_t i, bool direction)
  {
    currentEdge_.set((*distanceMatrix_)[rowIndex_][columnIndex_].first,
                     rowIndex_,
                     EdgeFiltrationTransformer::apply_modifier((*epsilonValues_)[i]),
                     direction);
  }

  /**
   * @brief Update the current edge as a vertex to insert.
   */
  void _update_edge_as_positive_vertex()
  {
    currentEdge_.set(epsilonIndex_ + 1,
                     epsilonIndex_ + 1,
                     EdgeFiltrationTransformer::apply_modifier((*epsilonValues_)[epsilonIndex_]),
                     true);
  }

  /**
   * @brief Update the current edge as a vertex to remove.
   */
  void _update_edge_as_negative_vertex()
  {
    currentEdge_.set(rowIndex_ - 1, rowIndex_ - 1, -std::numeric_limits<Filtration_value>::infinity(), false);
  }
};

/**
 * @class Oscillating_rips_edge_iterator_range oscillating_rips_edge_ranges.h \
 * gudhi/Zigzag_persistence/oscillating_rips_edge_ranges.h
 * @ingroup zigzag_persistence
 *
 * @brief Ordered range of @ref Zigzag_edge in the oscillating Rips filtration generated by the given parameters.
 * Even though it is called "edge", the simplex can also be a vertex if the two ends in @ref Zigzag_edge have the
 * same value. That means, that the range corresponds in fact to the filtration with maximal dimension 1.
 *
 * @details The range is forward traversal only, as the value of the edge is computed on the fly at each increment
 * and information of other edges are not stored. This is useful when the filtration is very long. Note that the
 * distance matrix and the filtration values are stored however. Note also that each copy of the same iterator will
 * increment simultaneously. That is, for example:
 * ```
 * Oscillating_rips_edge_iterator_range r(nu, mu, ...);
 * auto it1 = r.begin();
 * auto it2 = r.begin();
 * auto it3 = it1;
 * ++it1;
 * ++it2; ++it2;
 * // it1 and it2 are independent, so both have different values now: it1 points to the second edge, where it2 points
 * // to the third edge.
 * // but it3 is the copy of it1 and therefore, even if it3 was not explicitly incremented, it will still point to
 * // the second edge and not the first anymore.
 * ```
 * If a more flexible range is needed, use @ref Oscillating_rips_edge_vector_range_constructor::make_range instead.
 * It will construct a std::vector of @ref Zigzag_edge "".
 * 
 * @tparam Filtration_value Filtration value type. Should be compatible with the edge modifier.
 * @tparam EdgeFiltrationTransformer Modifier for the edge filtration values. If no modifications are wanted,
 * use @ref Identity_edge_modifier. Default value: @ref Identity_edge_modifier.
 */
template <typename Filtration_value, class EdgeFiltrationTransformer = Identity_edge_modifier>
class Oscillating_rips_edge_iterator_range
{
 public:
  /**
   * @class Oscillating_rips_edge_iterator oscillating_rips_edge_ranges.h \
   * gudhi/Zigzag_persistence/oscillating_rips_edge_ranges.h
   * @brief Custom iterator over the edges of an oscillating rips filtration.
   *
   * Category: LegacyInputIterator.
   *
   * @warning Each **copy** of the same iterator is pointing to the same base and will therefore update
   * **simultaneously**. This is to make the iterators copyable in the first place. If each copy would have its own
   * base, a copy would be too heavy to build without caution. Note that the `begin()` method of
   * @ref Oscillating_rips_edge_iterator_range does **not** return copies of the same iterator.
   */
  class Oscillating_rips_edge_iterator : public boost::iterator_facade<Oscillating_rips_edge_iterator,
                                                                       const Zigzag_edge<Filtration_value>&,
                                                                       boost::forward_traversal_tag>
  {
   public:
    /**
     * @brief Constructor.
     * 
     * @param base Pointer to an @ref Oscillating_rips_edge_iterator_base instantiation.
     */
    Oscillating_rips_edge_iterator(
        Oscillating_rips_edge_iterator_base<Filtration_value, EdgeFiltrationTransformer>* base)
        : base_iterator_(base)
    {}

    /**
     * @brief Default constructor. Equivalent to the end iterator, but has a high chance to seg fault if incremented
     * or dereferenced.
     */
    Oscillating_rips_edge_iterator() : base_iterator_(nullptr) {}

   private:
    // mandatory for the boost::iterator_facade inheritance.
    friend class boost::iterator_core_access;

    /**
     * @brief Pointer to heavy iterator, to avoid copying it.
     */
    std::shared_ptr<Oscillating_rips_edge_iterator_base<Filtration_value, EdgeFiltrationTransformer> > base_iterator_;

    /**
     * @brief Mandatory for the boost::iterator_facade inheritance. Indicates if to iterators are equal.
     *
     * @param other Iterator to compare.
     * @return True, the iterators are pointing to the same position.
     * @return False, otherwise.
     */
    bool equal(Oscillating_rips_edge_iterator const& other) const
    {
      if (base_iterator_ == nullptr) {
        if (other.base_iterator_ == nullptr) return true;
        return other.base_iterator_->equal(
            Oscillating_rips_edge_iterator_base<Filtration_value, EdgeFiltrationTransformer>());
      }
      return base_iterator_->equal(*other.base_iterator_);
    }

    /**
     * @brief Mandatory for the boost::iterator_facade inheritance. Returns the value of the dereferenced iterator.
     *
     * @return Current edge.
     */
    const Zigzag_edge<Filtration_value>& dereference() const { return base_iterator_->dereference(); }

    /**
     * @brief Mandatory for the boost::iterator_facade inheritance. Increments the iterator.
     */
    void increment() { base_iterator_->increment(); }
  };

  /**
   * @brief Default constructor. The range will be empty.
   */
  Oscillating_rips_edge_iterator_range() : nu_(0), mu_(0) {}

  /**
   * @brief Constructor. Initializes all necessary data to deduce the vertices and edges of the filtration.
   * See the @ref zigzagrips "introduction page" for more details about the arguments.
   * 
   * @tparam PointRange Point range type.
   * @tparam DistanceFunction Type of the distance function.
   * @param nu Lower multiplier.
   * @param mu Upper multiplier.
   * @param points Point cloud as a range.The format of a point has to correspond to the input format of the
   * distance function.
   * @param distance Distance function. Has to take two points as it from the range @p points as input parameters
   * and return the distance between those points.
   * @param orderPolicy Order policy for the points. Can be either
   * @ref Oscillating_rips_vertex_order_policy::FARTHEST_POINT_ORDERING,
   * @ref Oscillating_rips_vertex_order_policy::ALREADY_ORDERED or
   * @ref Oscillating_rips_vertex_order_policy::RANDOM_POINT_ORDERING.
   */
  template <typename PointRange, typename DistanceFunction>
  Oscillating_rips_edge_iterator_range(
      Filtration_value nu,
      Filtration_value mu,
      const PointRange& points,
      DistanceFunction&& distance,
      Oscillating_rips_vertex_order_policy orderPolicy = Oscillating_rips_vertex_order_policy::FARTHEST_POINT_ORDERING)
      : nu_(nu), mu_(mu)
  {
    GUDHI_CHECK((nu <= mu) && (nu >= 0), "Invalid parameters mu and nu");
    Oscillating_rips_initializer<Filtration_value>::initialize(
        epsilonValues_, distanceMatrix_, points, distance, orderPolicy);
  }

  /**
   * @brief Constructor. Initializes all necessary data to deduce the vertices and edges of the filtration.
   * See the @ref zigzagrips "introduction page" for more details about the arguments.
   * 
   * @tparam PointRange Point range type.
   * @tparam DistanceFunction Type of the distance function.
   * @param nu Lower multiplier.
   * @param mu Upper multiplier.
   * @param orderedPoints Point cloud already ordered in filtration order. The format of a point has to correspond to
   * the input format of the distance function.
   * @param distance Distance function. Has to take two points as it from the range @p points as input parameters
   * and return the distance between those points.
   * @param epsilonValues Wanted epsilon values. Should be decreasing and end with 0.
   */
  template <typename PointRange, typename DistanceFunction>
  Oscillating_rips_edge_iterator_range(Filtration_value nu,
                                       Filtration_value mu,
                                       const PointRange& orderedPoints,
                                       DistanceFunction&& distance,
                                       const std::vector<Filtration_value>& epsilonValues)
      : epsilonValues_(epsilonValues), nu_(nu), mu_(mu)
  {
    GUDHI_CHECK((nu <= mu) && (nu >= 0), "Invalid parameters mu and nu");
    Oscillating_rips_initializer<Filtration_value>::initialize(distanceMatrix_, orderedPoints, distance);
  }

  /**
   * @brief Returns the begin iterator of a the range.
   *
   * @warning Forward traversal only. And each copy of an iterator will increment simultaneously. That is, for example:
   * ```
   * Oscillating_rips_edge_iterator_range r(nu, mu, ...);
   * auto it1 = r.begin();
   * auto it2 = r.begin();
   * auto it3 = it1;
   * ++it1;
   * ++it2; ++it2;
   * // it1 and it2 are independent, so both have different values now: it1 points to the second edge, where it2 points
   * // to the third edge.
   * // but it3 is the copy of it1 and therefore, even if it3 was not explicitly incremented, it will still point to
   * // the second edge and not the first anymore.
   * ```
   */
  Oscillating_rips_edge_iterator begin()
  {
    // shared pointer on the other side will take ownership
    // enables begin() to be called several times without invalidating precedent iterators
    // still have the inconvenience that all copies of a same iterator (sharing the same base) increment simultaneously
    return Oscillating_rips_edge_iterator(
        new Oscillating_rips_edge_iterator_base<Filtration_value, EdgeFiltrationTransformer>(
            nu_, mu_, &epsilonValues_, &distanceMatrix_));
  }

  /**
   * @brief Returns the end iterator of a the range.
   */
  Oscillating_rips_edge_iterator end() { return endIt_; }

  /**
   * @brief Sets the multipliers of the range. Old iterators are not invalidated and will still point to the old
   * filtration. Only new iterators will take into account the changes. Epsilon values and the distance matrix will
   * remain unchanged.
   */
  void set_multipliers(Filtration_value nu, Filtration_value mu)
  {
    nu_ = nu;
    mu_ = mu;
  }

 private:
  std::vector<Filtration_value> epsilonValues_;                                 /**< Epsilon values. */
  std::vector<std::vector<std::pair<int, Filtration_value> > > distanceMatrix_; /**< Distance matrix. */
  Filtration_value nu_;                                                         /**< Lower multiplier. */
  Filtration_value mu_;                                                         /**< Upper multiplier. */
  /**
   * @brief End iterator. Does not depend on any parameter and can therefore be shared.
   */
  inline static const Oscillating_rips_edge_iterator endIt_ = Oscillating_rips_edge_iterator(
      new Oscillating_rips_edge_iterator_base<Filtration_value, EdgeFiltrationTransformer>());
};

/**
 * @class Oscillating_rips_edge_vector_range_constructor oscillating_rips_edge_ranges.h \
 * gudhi/Zigzag_persistence/oscillating_rips_edge_ranges.h
 * @ingroup zigzag_persistence
 *
 * @brief Constructor class for standard vectors of @ref Zigzag_edge corresponding to the oscillating Rips filtration
 * generated by the given parameters. Even though it is called "edge", the simplex can also be a vertex if the two ends
 * in @ref Zigzag_edge have the same value. That means, that the range corresponds in fact to the filtration with
 * maximal dimension 1.
 * 
 * @tparam Filtration_value Filtration value type. Should be compatible with the edge modifier.
 * @tparam EdgeFiltrationTransformer Modifier for the edge filtration values. If no modifications are wanted,
 * use @ref Identity_edge_modifier. Default value: @ref Identity_edge_modifier.
 */
template <typename Filtration_value, class EdgeFiltrationTransformer = Identity_edge_modifier>
class Oscillating_rips_edge_vector_range_constructor
{
 public:
  /**
   * @brief Builds the oscillating Rips filtration generated by the given parameters with maximal dimension 1 as a
   * standard vector of @ref Zigzag_edge "". A vertex is represented by an edge whose two ends have the same value.
   * See the @ref zigzagrips "introduction page" for more details about the arguments.
   * 
   * @tparam PointRange Point range type.
   * @tparam DistanceFunction Type of the distance function.
   * @param nu Lower multiplier.
   * @param mu Upper multiplier.
   * @param points Point cloud as a range.The format of a point has to correspond to the input format of the
   * distance function.
   * @param distance Distance function. Has to take two points as it from the range @p points as input parameters
   * and return the distance between those points.
   * @param orderPolicy Order policy for the points. Can be either
   * @ref Oscillating_rips_vertex_order_policy::FARTHEST_POINT_ORDERING,
   * @ref Oscillating_rips_vertex_order_policy::ALREADY_ORDERED or
   * @ref Oscillating_rips_vertex_order_policy::RANDOM_POINT_ORDERING.
   */
  template <typename PointRange, typename DistanceFunction>
  static std::vector<Zigzag_edge<Filtration_value> > make_range(
      Filtration_value nu,
      Filtration_value mu,
      const PointRange& points,
      DistanceFunction&& distance,
      Oscillating_rips_vertex_order_policy orderPolicy = Oscillating_rips_vertex_order_policy::FARTHEST_POINT_ORDERING)
  {
    GUDHI_CHECK((nu <= mu) && (nu >= 0), "Invalid parameters mu and nu");

    std::vector<Filtration_value> epsilonValues;
    std::vector<std::vector<std::pair<int, Filtration_value> > > distanceMatrix;

    Oscillating_rips_initializer<Filtration_value>::initialize(
        epsilonValues, distanceMatrix, points, distance, orderPolicy);

    return _compute_vector_range(EdgeFiltrationTransformer::apply_inverse_modifier(nu),
                                 EdgeFiltrationTransformer::apply_inverse_modifier(mu),
                                 distanceMatrix,
                                 epsilonValues);
  }

  /**
   * @brief Builds the oscillating Rips filtration generated by the given parameters with maximal dimension 1 as a
   * standard vector of @ref Zigzag_edge "". A vertex is represented by an edge whose two ends have the same value.
   * See the @ref zigzagrips "introduction page" for more details about the arguments.
   * 
   * @tparam PointRange Point range type.
   * @tparam DistanceFunction Type of the distance function.
   * @param nu Lower multiplier.
   * @param mu Upper multiplier.
   * @param orderedPoints Point cloud already ordered in filtration order. The format of a point has to correspond to
   * the input format of the distance function.
   * @param distance Distance function. Has to take two points as it from the range @p points as input parameters
   * and return the distance between those points.
   * @param epsilonValues Wanted epsilon values. Should be decreasing and end with 0.
   */
  template <typename PointRange, typename DistanceFunction>
  static std::vector<Zigzag_edge<Filtration_value> > make_range(Filtration_value nu,
                                                                Filtration_value mu,
                                                                const PointRange& orderedPoints,
                                                                DistanceFunction&& distance,
                                                                const std::vector<Filtration_value>& epsilonValues)
  {
    GUDHI_CHECK((nu <= mu) && (nu >= 0), "Invalid parameters mu and nu");
    GUDHI_CHECK(
        orderedPoints.size() == epsilonValues.size(),
        std::invalid_argument("Epsilon values should be initialized and have the same size than the point container."));

    std::vector<std::vector<std::pair<int, Filtration_value> > > distanceMatrix;

    Oscillating_rips_initializer<Filtration_value>::initialize(distanceMatrix, orderedPoints, distance);

    return _compute_vector_range(EdgeFiltrationTransformer::apply_inverse_modifier(nu),
                                 EdgeFiltrationTransformer::apply_inverse_modifier(mu),
                                 distanceMatrix,
                                 epsilonValues);
  }

 private:
  /**
   * @brief Computes and return a vector with all edges in order of the oscillating Rips filtration.
   */
  static std::vector<Zigzag_edge<Filtration_value> > _compute_vector_range(
      Filtration_value nu,
      Filtration_value mu,
      const std::vector<std::vector<std::pair<int, Filtration_value> > >& distanceMatrix,
      const std::vector<Filtration_value>& epsilonValues)
  {
    std::vector<Zigzag_edge<Filtration_value> > edgeFiltration;
    auto n = epsilonValues.size();

    // edgesAdded[i] (resp. edgesRemoved[i]) == list of edges (i,j), with j<i, added (resp. removed) at eps_i
    // we also put there (later) vertices that are added. Note that vertices are removed
    // only at the very last step of the oscillating Rips filtration.
    std::vector<std::vector<Zigzag_edge<Filtration_value> > > edgesAdded, edgesRemoved;

    std::size_t number_of_arrows = _compute_edges(nu, mu, epsilonValues, distanceMatrix, edgesAdded, edgesRemoved);

    // Now, sort edges according to lengths, and put everything in edgeFiltration
    edgeFiltration.clear();
    edgeFiltration.reserve(number_of_arrows + n);  // count edges + vertices additions

    // initialise R({p_0}, +infinity)
    edgeFiltration.emplace_back(0,
                                0,  // add vertex p_0,+infty
                                std::numeric_limits<Filtration_value>::infinity(),
                                true);
    // epsilonValues[0], true);

    for (std::size_t i = 0; i < n - 1; ++i) {  // all ascending arrows eps_i
      // add p_{i+1},eps_i
      edgeFiltration.emplace_back(i + 1, i + 1, EdgeFiltrationTransformer::apply_modifier(epsilonValues[i]), true);
      for (auto edgeIt = edgesAdded[i].begin(); edgeIt != edgesAdded[i].end(); ++edgeIt) {
        edgeFiltration.push_back(*edgeIt);
      }
      for (auto edgeIt = edgesRemoved[i].rbegin();  // longest first
           edgeIt != edgesRemoved[i].rend();
           ++edgeIt) {
        edgeFiltration.push_back(*edgeIt);
      }
    }
    for (int i = n - 1; i >= 0; --i) {
      edgeFiltration.emplace_back(i, i, -std::numeric_limits<Filtration_value>::infinity(), false);
    }

    _canonically_sort_edges(edgeFiltration);

    return edgeFiltration;
  }

  /**
   * @brief Computes the edges that are added and the edges that are removed and stores them in two separate containers.
   * The edge container @f$e@f$ stores the @f$i^{th}@f$ edge induced by vertex @f$j@f$ at @f$e[j][i]@f$.
   *
   * @param[in] nu Lower multiplier.
   * @param[in] mu Upper multiplier.
   * @param[in] epsilonValues Epsilon values in decreasing order finishing with 0.
   * @param[in] distanceMatrix Spare distance matrix.
   * @param[out] edgesAdded Container for positive edges.
   * @param[out] edgesRemoved Container for negative edges.
   *
   * @return Total number of edges.
   */
  static std::size_t _compute_edges(Filtration_value nu,
                                    Filtration_value mu,
                                    const std::vector<Filtration_value>& epsilonValues,
                                    const std::vector<std::vector<std::pair<int, Filtration_value> > >& distanceMatrix,
                                    std::vector<std::vector<Zigzag_edge<Filtration_value> > >& edgesAdded,
                                    std::vector<std::vector<Zigzag_edge<Filtration_value> > >& edgesRemoved)
  {
    std::size_t number_of_arrows = 0;
    auto n = epsilonValues.size();
    edgesAdded.resize(n);
    edgesRemoved.resize(n);

    // edgesAdded[i] must contain all new edges (k,j), 0 <= k < j <= i+1,
    // inserted in inclusion:
    // R({p_0, ... , p_i}, nu * eps_i) -> R({p_0, ... , p_i, p_i+1 }, mu * eps_i)
    // and
    // edgesRemoved[i] must contain all edges (k,j), 0 <= k < j <= i+1,
    // removed in backward inclusion:
    // R({p_0, ... , p_{i+1}}, mu * eps_i) <- R({p_0, ... , p_{i+1}}, nu * eps_{i+1})
#ifdef GUDHI_USE_TBB
    // no need to consider the case i=n-1 in an oscillating Rips filtration
    tbb::parallel_for(std::size_t(0), n - 1, [&](std::size_t i) {
      typename std::vector<std::pair<int, Filtration_value> >::const_iterator it;
      //----edgesAdded[i]:
      // consider first all edges added in inclusion:
      // R({p_0, ... , p_i}, nu * eps_i) -> R({p_0, ... , p_i}, mu * eps_i),
      // i.e., all (p_j,p_k) with 0 <= k < j <= i with
      //                                           nu eps_i < d(p_j,p_k) <= mu eps_i
      // these new edges get filtration value epsilonValues[i]
      for (std::size_t j = 1; j <= i; ++j) {
        // get very first edge (k,j), over all k<j, strictly longer than  mu * eps_i
        // distanceMatrix[j][k] = d(p_j,p_k) with k<j
        it = std::upper_bound(distanceMatrix[j].begin(),
                              distanceMatrix[j].end(),
                              std::pair<int, Filtration_value>(n, mu * epsilonValues[i]),
                              typename Oscillating_rips_initializer<Filtration_value>::Point_distance_comp());

        while (it != distanceMatrix[j].begin()) {
          --it;
          // if edge already in R({p_0, ... , p_i}, nu * eps_i), stop
          if (it->second <= nu * epsilonValues[i]) {
            break;
          }
          edgesAdded[i].emplace_back(it->first, j, EdgeFiltrationTransformer::apply_modifier(epsilonValues[i]), true);
          ++number_of_arrows;
        }
      }
      // now consider all edges added in inclusion:
      // R({p_0, ... , p_i}, mu * eps_i) -> R({p_0, ... , p_i, p_i+1}, mu * eps_i)
      // i.e., all (p_j,p_i+1) with 0 <= j <= i with       d(p_j,p_i+1) <= mu eps_i
      // these new edges get filtration value epsilonValues[i]
      // first strictly longer edge
      it = std::upper_bound(distanceMatrix[i + 1].begin(),
                            distanceMatrix[i + 1].end(),
                            std::pair<int, Filtration_value>(n, mu * epsilonValues[i]),
                            typename Oscillating_rips_initializer<Filtration_value>::Point_distance_comp());

      while (it != distanceMatrix[i + 1].begin()) {
        --it;
        edgesAdded[i].emplace_back(it->first, i + 1, EdgeFiltrationTransformer::apply_modifier(epsilonValues[i]), true);
        ++number_of_arrows;
      }

      //----edgesRemoved[i]:
      // consider all edges removed in
      // R({p_0, ... , p_{i+1}}, mu * eps_i) <- R({p_0, ... , p_{i+1}}, nu * eps_{i+1})
      // i.e., all edges (p_k,p_j), 0<=k<j<=i+1, such that
      // nu eps_{i+1} < d(p_k,p_j) <= mu eps_i
      // these new edges get filtration value epsilonValues[i]
      for (std::size_t j = 1; j <= i + 1; ++j) {
        // get very first edge (k,j), over all k<j, strictly longer than  mu * eps_i
        // distanceMatrix[j][k] = d(p_j,p_k) with k<j
        it = std::upper_bound(distanceMatrix[j].begin(),
                              distanceMatrix[j].end(),
                              std::pair<int, Filtration_value>(n, mu * epsilonValues[i]),
                              typename Oscillating_rips_initializer<Filtration_value>::Point_distance_comp());

        while (it != distanceMatrix[j].begin()) {
          --it;
          // when reading an edge in R({p_0, ... , p_{i+1}}, nu * eps_{i+1}), stop
          if (it->second <= nu * epsilonValues[i + 1]) {
            break;
          }
          edgesRemoved[i].emplace_back(
              it->first, j, EdgeFiltrationTransformer::apply_modifier(epsilonValues[i]), false);
          ++number_of_arrows;
        }
      }
    });
#else  // GUDHI_USE_TBB not defined

    typename std::vector<std::pair<int, Filtration_value> >::const_iterator it;

    for (std::size_t i = 0; i < n - 1; ++i) {
      //----edgesAdded[i]:
      // consider first all edges added in inclusion:
      // R({p_0, ... , p_i}, nu * eps_i) -> R({p_0, ... , p_i}, mu * eps_i),
      // i.e., all (p_j,p_k) with 0 <= k < j <= i with
      //                                           nu eps_i < d(p_j,p_k) <= mu eps_i
      // these new edges get filtration value epsilonValues[i]
      for (std::size_t j = 1; j <= i; ++j) {
        // get very first edge (k,j), over all k<j, strictly longer than  mu * eps_i
        // distanceMatrix[j][k] = d(p_j,p_k) with k<j
        it = std::upper_bound(distanceMatrix[j].begin(),
                              distanceMatrix[j].end(),
                              std::pair<int, Filtration_value>(n, mu * epsilonValues[i]),
                              typename Oscillating_rips_initializer<Filtration_value>::Point_distance_comp());

        while (it != distanceMatrix[j].begin()) {
          --it;
          // if edge already in R({p_0, ... , p_i}, nu * eps_i), stop
          if (it->second <= nu * epsilonValues[i]) {
            break;
          }
          edgesAdded[i].emplace_back(it->first, j, EdgeFiltrationTransformer::apply_modifier(epsilonValues[i]), true);
          ++number_of_arrows;
        }
      }
      // now consider all edges added in inclusion:
      // R({p_0, ... , p_i}, mu * eps_i) -> R({p_0, ... , p_i, p_i+1}, mu * eps_i)
      // i.e., all (p_j,p_i+1) with 0 <= j <= i with       d(p_j,p_i+1) <= mu eps_i
      // these new edges get filtration value epsilonValues[i]
      // first strictly longer edge
      it = std::upper_bound(distanceMatrix[i + 1].begin(),
                            distanceMatrix[i + 1].end(),
                            std::pair<int, Filtration_value>(n, mu * epsilonValues[i]),
                            typename Oscillating_rips_initializer<Filtration_value>::Point_distance_comp());
      while (it != distanceMatrix[i + 1].begin()) {
        --it;
        edgesAdded[i].emplace_back(it->first, i + 1, EdgeFiltrationTransformer::apply_modifier(epsilonValues[i]), true);
        ++number_of_arrows;
      }

      //----edgesRemoved[i]:
      // consider all edges removed in
      // R({p_0, ... , p_{i+1}}, mu * eps_i) <- R({p_0, ... , p_{i+1}}, nu * eps_{i+1})
      // i.e., all edges (p_k,p_j), 0<=k<j<=i+1, such that
      // nu eps_{i+1} < d(p_k,p_j) <= mu eps_i
      // these new edges get filtration value epsilonValues[i]
      for (std::size_t j = 1; j <= i + 1; ++j) {
        // get very first edge (k,j), over all k<j, strictly longer than  mu * eps_i
        // distanceMatrix[j][k] = d(p_j,p_k) with k<j
        it = std::upper_bound(distanceMatrix[j].begin(),
                              distanceMatrix[j].end(),
                              std::pair<int, Filtration_value>(n, mu * epsilonValues[i]),
                              typename Oscillating_rips_initializer<Filtration_value>::Point_distance_comp());

        while (it != distanceMatrix[j].begin()) {
          --it;
          // when reading an edge in R({p_0, ... , p_{i+1}}, nu * eps_{i+1}), stop
          if (it->second <= nu * epsilonValues[i + 1]) {
            break;
          }
          edgesRemoved[i].emplace_back(
              it->first, j, EdgeFiltrationTransformer::apply_modifier(epsilonValues[i]), false);
          ++number_of_arrows;
        }
      }
    }
#endif

    return number_of_arrows;
  }

  /**
   * @brief Sorts canonically the edges: as much as possible, edges should be removed in
   * the reverse order of their insertion. We decide to insert shortest edges first,
   * with increasing lexicographical order, and remove larger edges first, with
   * decreasing lexicographic order.
   *
   * @param edges Edges to sort.
   */
  static void _canonically_sort_edges(std::vector<Zigzag_edge<Filtration_value> >& edges)
  {
    // filtration then dimension, then lex order for insertion
    auto edge_cmp = [](const Zigzag_edge<Filtration_value>& e1, const Zigzag_edge<Filtration_value>& e2) {
      if (e1.get_filtration_value() != e2.get_filtration_value()) {
        return e1.get_filtration_value() < e2.get_filtration_value();  // lower fil first
      }

      if (e1.get_smallest_vertex() == e1.get_biggest_vertex()) {  // e1 is a vertex, -> put vertices first
        if (e2.get_smallest_vertex() == e2.get_biggest_vertex()) {
          return e1.get_smallest_vertex() < e2.get_smallest_vertex();  //-> vertex of lower label
        } else {
          return true;  //-> always vertices before edges
        }
      }
      // e1 is an edge
      if (e2.get_smallest_vertex() == e2.get_biggest_vertex()) {
        return false;  // e2 vertex, -> put it first
      }
      // both are edges, lexicographic compare
      if (e1.get_smallest_vertex() != e2.get_smallest_vertex()) {
        return e1.get_smallest_vertex() < e2.get_smallest_vertex();  // lex order
      }
      if (e1.get_biggest_vertex() != e2.get_biggest_vertex()) {
        return e1.get_biggest_vertex() < e2.get_biggest_vertex();
      }
      return false;  // equality
    };

    // the inverse ordering for deletions
    auto inv_edge_cmp = [&](const Zigzag_edge<Filtration_value>& e1, const Zigzag_edge<Filtration_value>& e2) {
      if (e1.get_smallest_vertex() == e2.get_smallest_vertex() && e1.get_biggest_vertex() == e2.get_biggest_vertex()) {
        return false;  // equality => false
      }
      return !(edge_cmp(e1, e2));  // reverse order
    };

    // sort sequences of inclusions of same filtration with edge_cmp
    // sort sequences of removals of same filtration with inv_edge_cmp
    auto beg = edges.begin();
    auto end = edges.begin();
    auto curr_fil = beg->get_filtration_value();
    auto curr_type = beg->get_direction();
    while (beg != edges.end()) {
      while (end != edges.end() && end->get_filtration_value() == curr_fil && end->get_direction() == curr_type) {
        ++end;
      }
      if (curr_type) {  // sequence of insertions
#ifdef GUDHI_USE_TBB
        tbb::parallel_sort(beg, end, edge_cmp);
#else
        std::sort(beg, end, edge_cmp);
#endif
      } else {  // sequence of removals
#ifdef GUDHI_USE_TBB
        tbb::parallel_sort(beg, end, inv_edge_cmp);
#else
        std::sort(beg, end, inv_edge_cmp);
#endif
      }
      beg = end;
      curr_fil = beg->get_filtration_value();
      curr_type = beg->get_direction();
    }
  }
};

}  // namespace zigzag_persistence
}  // namespace Gudhi

#endif  // ZIGZAG_OSCILLATING_RIPS_EDGE_RANGES_H_
