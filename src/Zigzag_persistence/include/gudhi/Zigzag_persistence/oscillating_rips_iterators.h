/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria and Hannah Schreiber
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

 /**
 * @file oscillating_rips_iterators.h
 * @author Clément Maria, Hannah Schreiber
 * @brief Contains the implementation of the @ref Gudhi::zigzag_persistence::Zigzag_edge class, 
 * @ref Gudhi::zigzag_persistence::Identity_edge_modifier class, 
 * @ref Gudhi::zigzag_persistence::Square_root_edge_modifier class, 
 * @ref Gudhi::zigzag_persistence::Oscillating_rips_edge_range class and 
 * @ref Gudhi::zigzag_persistence::Oscillating_rips_simplex_range class.
 */

#ifndef ZIGZAG_OSCILLATING_RIPS_ITERATORS_H_
#define ZIGZAG_OSCILLATING_RIPS_ITERATORS_H_

#include <cmath>
#include <cstddef>
#include <vector>
#include <algorithm>
#include <tuple>

#include <boost/iterator/iterator_facade.hpp>

#ifdef GUDHI_USE_TBB
#include <tbb/tbb.h>
#endif

#include <gudhi/choose_n_farthest_points.h>
#include <gudhi/pick_n_random_points.h>
#include <gudhi/Debug_utils.h>

namespace Gudhi {
namespace zigzag_persistence {

/**
 * @class Zigzag_edge oscillating_rips_iterators.h gudhi/Zigzag_persistence/oscillating_rips_iterators.h
 * @brief Edge structure for the oscilliating Rips filtration.
 * 
 * @ingroup zigzag_persistence
 *
 * @details The edges of the filtration are computed first and the remaining simplices are deduced from them. 
 * Is also used to represent the vertices for technical reasons, by giving both vertices the same value.
 * 
 * @tparam Filtration_value Type of the filtration value.
 */
template <typename Filtration_value>
class Zigzag_edge {
 public:
  /**
   * @brief Constructor.
   * 
   * @param u First boundary vertex ID
   * @param v Second boundary vertex ID. If same than @p u, the edge is considered a vertex.
   * @param fil Filtration value
   * @param direction If true, forwards. If false, backwards.
   */
  Zigzag_edge(int u, int v, Filtration_value fil, bool direction)
      : u_(u), v_(v), fil_(fil), direction_(direction) {
        if (u > v) std::swap(u_, v_);
      }

  /**
   * @brief Default constructor. Initialize everything to zero and the direction to true.
   */
  Zigzag_edge() : u_(0), v_(0), fil_(0), direction_(true) {}

  /**
   * @brief Returns vertex with smaller ID.
   * 
   * @return Smallest ID among the boundary vertices.
   */
  int get_smallest_vertex() const { return u_; }
  /**
   * @brief Returns vertex with bigger ID.
   * 
   * @return Biggest ID among the boundary vertices.
   */
  int get_biggest_vertex() const { return v_; }
  /**
   * @brief Returns the filtration value of the edge.
   * 
   * @return Filtration value of the edge.
   */
  Filtration_value get_filtration_value() const { return fil_; }
  /**
   * @brief Gives the direction of the arrow corresponding to the edge.
   * 
   * @return True, if forward, i.e., an insertion.
   * @return False, if backward, i.e., a removal.
   */
  bool get_direction() const { return direction_; }

  void set(int u, int v, Filtration_value fil, bool direction){
    u_ = u;
    v_ = v;
    fil_ = fil;
    direction_ = direction;
  }

  /**
   * @brief Equality test
   * 
   * @param e Edge to compare.
   * @return True, if both edges are equal.
   * @return False, if both edges are not equal.
   */
  bool operator==(const Zigzag_edge& e) const {
    return ((e.u_ == u_) && (e.v_ == v_) && (e.fil_ == fil_) && (e.direction_ == direction_));
  }

 private:
  int u_;                   /**< Smaller vertex. */
  int v_;                   /**< Bigger vertex. */
  Filtration_value fil_;    /**< Filtration value. */
  bool direction_;          /**< Direction. True = forward, false = backward. */
};

/**
 * @class Identity_edge_modifier oscillating_rips_iterators.h gudhi/Zigzag_persistence/oscillating_rips_iterators.h
 * @brief Identity modifier, i.e., does nothing.
 *
 * @ingroup zigzag_persistence
 */
class Identity_edge_modifier {
 public:
  static constexpr bool isActive_ = false;  /**< Indicates that the modifier should be ignored. */

 private:
  /**
   * @brief Default constructor. Should not be called and therefore private.
   */
  Identity_edge_modifier() {}
};

/**
 * @class Square_root_edge_modifier oscillating_rips_iterators.h gudhi/Zigzag_persistence/oscillating_rips_iterators.h
 * @brief Modifier that square roots the filtration value of an edge.
 * 
 * @ingroup zigzag_persistence
 *
 * @details Useful in particular when geometric computations (edge length, etc) are 
 * run with squared Euclidean distance for performance.
 *
 * @tparam Filtration_value Filtration value type.
 */
template <typename Filtration_value>
class Square_root_edge_modifier {
 public:
  /**
   * @brief Returns the square root of the given value. The parameter it-self is not modified.
   * 
   * @param f Value to modify.
   * @return The modified value of @p f.
   */
  static Filtration_value apply_modifier(Filtration_value f) { return std::sqrt(f); }
  /**
   * @brief Returns the square of the given value. The parameter it-self is not modified.
   * 
   * @param f Value to modify.
   * @return The modified value of @p f.
   */
  static Filtration_value apply_inverse_modifier(Filtration_value f) { return f * f; }

  static constexpr bool isActive_ = true;   /**< Indicates that the modifier should not be ignored. */

 private:
  /**
   * @brief Default constructor. Should not be called and therefore private.
   */
  Square_root_edge_modifier() {}
};

/**
 * @class Oscillating_rips_edge_range oscillating_rips_iterators.h gudhi/Zigzag_persistence/oscillating_rips_iterators.h
 * @brief Gives access and computes different edge range types for an oscilliating Rips filtration.
 *
 * @ingroup zigzag_persistence
 *
 * @details There are two different kind of ranges: a vector range containing all computed edges and
 * a Boost range made from a custom iterator computing an edge on the fly at each incrementation.
 * The custom iterator is therefore only a forward iterator and can only be incremented.
 *
 * @tparam Filtration_value Filtration value type
 * @tparam EdgeModifier Modifier for the edge filtration values. If no modifications are wanted, 
 * use @ref Identity_edge_modifier. Default value: @ref Identity_edge_modifier.
 *
 * @warning As the custom iterator used for the Boost ranges stores possibly large ranges, avoid copying it. 
 * Use @ref get_iterator_range, @ref begin and @ref end wisely. If one needs to do something else than to simply 
 * iterate over the edges in order, then the vector range is for sure the better option.
 */
template <typename Filtration_value, class EdgeModifier = Identity_edge_modifier>
class Oscillating_rips_edge_range {
 public:
  /**
   * @brief Order policy for the points.
   */
  enum Order_policy { 
    ALREADY_ORDERED,            /**< The given range of points is considered ordered. */
    FARTHEST_POINT_ORDERING,    /**< The points are reordered using @ref Gudhi::subsampling::choose_n_farthest_points.*/
    RANDOM_POINT_ORDERING       /**< The points are shuffled randomly. */
  };

  /**
   * @class Oscillating_rips_edge_iterator oscillating_rips_iterators.h gudhi/Zigzag_persistence/oscillating_rips_iterators.h
   * @brief Custom iterator over the edges of an oscillating rips filtration.
   * 
   * It inherits from boost::iterator_facade.
   * 
   * The iterator stores the distance matrix and epsilon values to compute the current edge on the fly
   * at each incrementation.
   *
   * @warning As the iterator stores possibly large ranges, avoid copying it.
   */
  class Oscillating_rips_edge_iterator
      : public boost::iterator_facade<Oscillating_rips_edge_iterator, 
                                      const Zigzag_edge<Filtration_value>&,
                                      boost::forward_traversal_tag> {
   public:
    /**
     * @brief Constructor. Computes the distance matrix and the epsilon values it-self from the given parameters.
     *
     * @tparam PointRange Point range type.
     * @tparam DistanceFunction Type of the distance function.
     * @param nu Lower multiplicator.
     * @param mu Upper multiplicator.
     * @param points Point cloud as a range.The format of a point has to correspond to the input format of the
     * distance function.
     * @param distance Distance function. Has to take two points as it from the range @p points as input parameters
     * and return the distance between those points.
     * @param orderPolicy Order policy for the points. Can be either @ref Order_policy::FARTHEST_POINT_ORDERING,
     * @ref Order_policy::ALREADY_ORDERED or @ref Order_policy::RANDOM_POINT_ORDERING.
     */
    template <typename PointRange, typename DistanceFunction>
    Oscillating_rips_edge_iterator(Filtration_value nu, 
                                   Filtration_value mu, 
                                   const PointRange& points,
                                   DistanceFunction&& distance, 
                                   Order_policy orderPolicy)
        : nu_(nu),
          mu_(mu),
          currentEdge_(0, 0, std::numeric_limits<Filtration_value>::infinity(), true),
          epsilonIndex_(0),
          rowIndex_(1),
          inPositiveDirection_(true),
          insertVertex_(true) 
    {
      _initialize(nu_, mu_, epsilonValues_, distanceMatrix_, points, distance, orderPolicy);
      auto it =
          std::upper_bound(distanceMatrix_[1].begin(), distanceMatrix_[1].end(),
                           std::pair<int, Filtration_value>(distanceMatrix_.size(), mu_ * epsilonValues_[epsilonIndex_]),
                           Point_distance_comp());
      columnIndex_ = it - distanceMatrix_[1].begin();
    }

    /**
     * @brief Constructor. Takes already computed epsilon values as input, but assumes that the points are 
     * already ordered accordingly. Assumes also that the last epsilon value is 0.
     * 
     * @tparam PointRange Point range type.
     * @tparam DistanceFunction Type of the distance function.
     * @param nu Lower multiplicator.
     * @param mu Upper multiplicator.
     * @param orderedPoints Point cloud as an ordered range. The format of a point has to correspond to the input 
     * format of the distance function. The order of the points should be in correspondence with the order of 
     * the epsilon values.
     * @param distance Distance function. Has to take two points as it from the range @p points as input parameters
     * and return the distance between those points.
     * @param epsilonValues Epsilon values for the oscillating rips filtration. Should be in decreasing order.
     * And the last value should be 0.
     */
    template <typename PointRange, typename DistanceFunction>
    Oscillating_rips_edge_iterator(Filtration_value nu, 
                                   Filtration_value mu, 
                                   const PointRange& orderedPoints,
                                   DistanceFunction&& distance, 
                                   const std::vector<Filtration_value>& epsilonValues)
        : epsilonValues_(epsilonValues),
          nu_(nu),
          mu_(mu),
          currentEdge_(0, 0, std::numeric_limits<Filtration_value>::infinity(), true),
          epsilonIndex_(0),
          rowIndex_(1),
          inPositiveDirection_(true),
          insertVertex_(true)
    {
      GUDHI_CHECK(orderedPoints.size() == epsilonValues.size(),
                  "The number of points and the number of epsilon values should match.");
      GUDHI_CHECK((nu <= mu) && (nu >= 0), "Invalid parameters mu and nu");

      if constexpr (EdgeModifier::isActive_) {
        nu_ = EdgeModifier::apply_inverse_modifier(nu);
        mu_ = EdgeModifier::apply_inverse_modifier(mu);
      }

      // compute the distance matrix
      distanceMatrix_ = _compute_distance_matrix(orderedPoints, distance);

      auto it =
          std::upper_bound(distanceMatrix_[1].begin(), distanceMatrix_[1].end(),
                           std::pair<int, Filtration_value>(distanceMatrix_.size(), mu_ * epsilonValues_[epsilonIndex_]),
                           Point_distance_comp());
      columnIndex_ = it - distanceMatrix_[1].begin();
    }

    /**
     * @brief Default constructor. Corresponds to the end iterator.
     */
    Oscillating_rips_edge_iterator()
        : nu_(0),
          mu_(0),
          currentEdge_(0, 0, 0, true),
          epsilonIndex_(0),
          rowIndex_(0),
          columnIndex_(0),
          inPositiveDirection_(true),
          insertVertex_(true) {}

   private:
    //mandatory for the boost::iterator_facade inheritance.
    friend class boost::iterator_core_access;

    std::vector<Filtration_value> epsilonValues_;   /**< Epsilon values. */
    std::vector<std::vector<std::pair<int, Filtration_value> > > distanceMatrix_;   /**< Distance matrix. */
    Filtration_value nu_;                           /**< Lower multiplicator. */
    Filtration_value mu_;                           /**< Upper multiplicator. */
    Zigzag_edge<Filtration_value> currentEdge_;     /**< Stores the current edge in the range. */
    size_t epsilonIndex_, rowIndex_, columnIndex_;  /**< Indices indicating the next position in the range. */
    bool inPositiveDirection_, insertVertex_;       /**< Next direction and indicates if next ``edge'' is a vertex. */

    /**
     * @brief Mandatory for the boost::iterator_facade inheritance. Indicates if to iterators are equal.
     * 
     * @param other Iterator to compare.
     * @return True, the iterators are pointing to the same position.
     * @return False, otherwise.
     */
    bool equal(Oscillating_rips_edge_iterator const& other) const {
      return rowIndex_ == other.rowIndex_ && currentEdge_ == other.currentEdge_;
    }

    /**
     * @brief Mandatory for the boost::iterator_facade inheritance. Returns the value of the dereferenced iterator.
     * 
     * @return Current edge.
     */
    const Zigzag_edge<Filtration_value>& dereference() const { return currentEdge_; }

    /**
     * @brief Mandatory for the boost::iterator_facade inheritance. Increments the iterator.
     * 
     */
    void increment() {
      if (epsilonIndex_ < distanceMatrix_.size() - 1) {
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
            if (epsilonIndex_ == distanceMatrix_.size() - 1) {
              //   _set_end();
              rowIndex_ = distanceMatrix_.size();
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
          if (epsilonIndex_ == distanceMatrix_.size() - 1) {
            rowIndex_ = distanceMatrix_.size() + 1;
            return;
          }
          insertVertex_ = true;
          inPositiveDirection_ = true;
          ++rowIndex_;
          _initialize_positive_col_index();
        }
        return;
      } if (rowIndex_ > 1){
        --rowIndex_;
        _update_edge_as_negative_vertex();
        return;
      }

      _set_end();
    }

    /**
     * @brief Set the iterator as the end iterator.
     */
    void _set_end(){
        rowIndex_ = 0;
        currentEdge_ = Zigzag_edge<Filtration_value>();
    }

    /**
     * @brief Initialize the column index for a new sequence of positive edges.
     */
    void _initialize_positive_col_index() {
      auto it =
          std::upper_bound(distanceMatrix_[rowIndex_].begin(), distanceMatrix_[rowIndex_].end(),
                           std::pair<int, Filtration_value>(distanceMatrix_.size(), mu_ * epsilonValues_[epsilonIndex_]),
                           Point_distance_comp());
      columnIndex_ = it - distanceMatrix_[rowIndex_].begin();
    }

    /**
     * @brief Initialize the column index for a new sequence of negative edges.
     */
    void _initialize_negative_col_index() {
      auto it =
          std::lower_bound(distanceMatrix_[rowIndex_].begin(), distanceMatrix_[rowIndex_].end(),
                           std::pair<int, Filtration_value>(0, nu_ * epsilonValues_[epsilonIndex_ + 1]),
                           Point_distance_comp());
      while (it != distanceMatrix_[rowIndex_].end() && it->second == nu_ * epsilonValues_[epsilonIndex_ + 1]) ++it;
      columnIndex_ = it - distanceMatrix_[rowIndex_].begin();
    }

    /**
     * @brief Indicates if the column index needs to be initialized or not.
     */
    bool _positive_col_index_is_not_valid() {
        return columnIndex_ == 0 ||
               (rowIndex_ != (epsilonIndex_ + 1) &&
                distanceMatrix_[rowIndex_][columnIndex_ - 1].second <= nu_ * epsilonValues_[epsilonIndex_]);
    }

    /**
     * @brief Indicates if the column index needs to be initialized or not.
     */
    bool _negative_col_index_is_not_valid() {
        return columnIndex_ == distanceMatrix_[rowIndex_].size() ||
               distanceMatrix_[rowIndex_][columnIndex_].second > mu_ * epsilonValues_[epsilonIndex_];
    }

    /**
     * @brief Updates the current edge.
     * 
     * @param i Epsilon value index.
     * @param direction Direction.
     */
    void _update_edge(size_t i, bool direction) {
        if constexpr (EdgeModifier::isActive_)
            currentEdge_.set(distanceMatrix_[rowIndex_][columnIndex_].first, rowIndex_,
                             EdgeModifier::apply_modifier(epsilonValues_[i]), direction);
        else
            currentEdge_.set(distanceMatrix_[rowIndex_][columnIndex_].first, rowIndex_, epsilonValues_[i], direction);
    }

    /**
     * @brief Update the current edge as a vertex to insert.
     */
    void _update_edge_as_positive_vertex() {
        if constexpr (EdgeModifier::isActive_)
            currentEdge_.set(epsilonIndex_ + 1, epsilonIndex_ + 1,
                             EdgeModifier::apply_modifier(epsilonValues_[epsilonIndex_]), true);
        else
            currentEdge_.set(epsilonIndex_ + 1, epsilonIndex_ + 1, epsilonValues_[epsilonIndex_], true);
    }

    /**
     * @brief Update the current edge as a vertex to remove.
     */
    void _update_edge_as_negative_vertex() {
        currentEdge_.set(rowIndex_ - 1, rowIndex_ - 1, -std::numeric_limits<Filtration_value>::infinity(), false);
    }
  };

  /**
   * @brief Computes and return a vector with all edges in order of the oscillating Rips filtration.
   *
   * @tparam PointRange Point range type.
   * @tparam DistanceFunction Type of the distance function.
   * @param nu Lower multiplicator.
   * @param mu Upper multiplicator.
   * @param points Point cloud as a range.The format of a point has to correspond to the input format of the
   * distance function.
   * @param distance Distance function. Has to take two points as it from the range @p points as input parameters
   * and return the distance between those points.
   * @param orderPolicy Order policy for the points. Can be either @ref Order_policy::FARTHEST_POINT_ORDERING,
   * @ref Order_policy::ALREADY_ORDERED or @ref Order_policy::RANDOM_POINT_ORDERING.
   *
   * @return Vector of @ref Zigzag_edge.
   */
  template <typename PointRange, typename DistanceFunction>
  static std::vector<Zigzag_edge<Filtration_value> > compute_vector_range(
      Filtration_value nu, 
      Filtration_value mu, 
      const PointRange& points, 
      DistanceFunction&& distance,
      Order_policy orderPolicy = Order_policy::FARTHEST_POINT_ORDERING) {
    std::vector<Zigzag_edge<Filtration_value> > edgeFiltration;
    std::vector<Filtration_value> epsilonValues;
    std::vector<std::vector<std::pair<int, Filtration_value> > > distanceMatrix;
    auto n = points.size();

    _initialize(nu, mu, epsilonValues, distanceMatrix, points, distance, orderPolicy);

    // edgesAdded[i] (resp. edgesRemoved[i]) == list of edges (i,j), with j<i, added (resp. removed) at eps_i
    // we also put there (later) vertices that are added. Note that vertices are removed
    // only at the very last step of the oRzz filtration.
    std::vector<std::vector<Zigzag_edge<Filtration_value> > > edgesAdded, edgesRemoved;

    // auto it = std::upper_bound(distanceMatrix[1].begin(), distanceMatrix[1].end(),
    //                              std::pair<int, Filtration_value>(distanceMatrix.size(), mu * epsilonValues[0]),
    //                              Point_distance_comp());
    // std::cout << "start vect colind: " << (it - distanceMatrix[1].begin()) << ", (" << it->first << ", " << it->second << "), (" << distanceMatrix[1].begin()->first << ", " << distanceMatrix[1].begin()->second << ")\n";
    size_t number_of_arrows = _compute_edges(nu, mu, epsilonValues, distanceMatrix, edgesAdded, edgesRemoved);

    // Now, sort edges according to lengths, and put everything in edgeFiltration
    edgeFiltration.clear();
    edgeFiltration.reserve(number_of_arrows + n);  // count edges + vertices additions

    // initialise R({p_0}, +infinity)
    edgeFiltration.emplace_back(0, 0,  // add vertex p_0,+infty
                                std::numeric_limits<Filtration_value>::infinity(), true);
    // epsilonValues[0], true);

    for (size_t i = 0; i < n - 1; ++i) {  // all ascending arrows eps_i
      if constexpr (EdgeModifier::isActive_) {
        edgeFiltration.emplace_back(i + 1, i + 1, 
                                    EdgeModifier::apply_modifier(epsilonValues[i]),
                                    true);  // add p_{i+1},eps_i
      } else {
        edgeFiltration.emplace_back(i + 1, i + 1, epsilonValues[i], true);  // add p_{i+1},eps_i
      }
      for (auto edg_it = edgesAdded[i].begin(); edg_it != edgesAdded[i].end(); ++edg_it) {
        edgeFiltration.push_back(*edg_it);
      }
      for (auto edg_it = edgesRemoved[i].rbegin();  // longest first
           edg_it != edgesRemoved[i].rend(); ++edg_it) {
        edgeFiltration.push_back(*edg_it);
      }
    }
    for (int i = n - 1; i >= 0; --i) {
      edgeFiltration.emplace_back(i, i, -std::numeric_limits<Filtration_value>::infinity(), false);
    }

    _canonically_sort_edges(edgeFiltration);

    return edgeFiltration;
  }

  /**
   * @brief Returns a boost::iterator_range from @ref Oscillating_rips_edge_iterator. The edges are computed 
   * on the fly at each incrementation. The iterator is a forward iterator only.
   * 
   * @tparam PointRange Point range type.
   * @tparam DistanceFunction Type of the distance function.
   * @param nu Lower multiplicator.
   * @param mu Upper multiplicator.
   * @param points Point cloud as a range.The format of a point has to correspond to the input format of the
   * distance function.
   * @param distance Distance function. Has to take two points as it from the range @p points as input parameters
   * and return the distance between those points.
   * @param orderPolicy Order policy for the points. Can be either @ref Order_policy::FARTHEST_POINT_ORDERING,
   * @ref Order_policy::ALREADY_ORDERED or @ref Order_policy::RANDOM_POINT_ORDERING.
   *
   * @return boost::iterator_range of @ref Oscillating_rips_edge_iterator.
   * 
   * @warning Avoid copying the iterators as they are heavier than usual iterators. If begin and end iterators 
   * are needed but not the structure of the range, use @ref begin and @ref end instead.
   */
  template <typename PointRange, typename DistanceFunction>
  static boost::iterator_range<Oscillating_rips_edge_iterator> get_iterator_range(
      Filtration_value nu, 
      Filtration_value mu, 
      const PointRange& points, 
      DistanceFunction&& distance,
      Order_policy orderPolicy = Order_policy::FARTHEST_POINT_ORDERING) 
  {
    return boost::iterator_range<Oscillating_rips_edge_iterator>(
        Oscillating_rips_edge_iterator(nu, mu, points, distance, orderPolicy), Oscillating_rips_edge_iterator());
  }

  /**
   * @brief Returns a boost::iterator_range from @ref Oscillating_rips_edge_iterator. The edges are computed
   * on the fly at each incrementation. The iterator is a forward iterator only.
   * Takes already computed epsilon values as input, but assumes that the points are already ordered accordingly. 
   * Assumes also that the last epsilon value is 0.
   *
   * @tparam PointRange Point range type.
   * @tparam DistanceFunction Type of the distance function.
   * @param nu Lower multiplicator.
   * @param mu Upper multiplicator.
   * @param orderedPoints Point cloud as an ordered range. The format of a point has to correspond to the input
   * format of the distance function. The order of the points should be in correspondence with the order of
   * the epsilon values.
   * @param distance Distance function. Has to take two points as it from the range @p points as input parameters
   * and return the distance between those points.
   * @param epsilonValues Epsilon values for the oscillating rips filtration. Should be in decreasing order.
   * And the last value should be 0.
   *
   * @return boost::iterator_range of @ref Oscillating_rips_edge_iterator.
   * 
   * @warning Avoid copying the iterators as they are heavier than usual iterators. If begin and end iterators 
   * are needed but not the structure of the range, use @ref begin and @ref end instead.
   */
  template <typename PointRange, typename DistanceFunction>
  static boost::iterator_range<Oscillating_rips_edge_iterator> get_iterator_range(
      Filtration_value nu, 
      Filtration_value mu, 
      const PointRange& orderedPoints,
      DistanceFunction&& distance, 
      const std::vector<Filtration_value>& epsilonValues) 
  {
    return boost::iterator_range<Oscillating_rips_edge_iterator>(
        Oscillating_rips_edge_iterator(nu, mu, orderedPoints, distance, epsilonValues), Oscillating_rips_edge_iterator());
  }

  /**
   * @brief Returns the begin iterator of a the range of edges based on @ref Oscillating_rips_edge_iterator.
   * 
   * @tparam PointRange Point range type.
   * @tparam DistanceFunction Type of the distance function.
   * @param nu Lower multiplicator.
   * @param mu Upper multiplicator.
   * @param points Point cloud as a range.The format of a point has to correspond to the input format of the
   * distance function.
   * @param distance Distance function. Has to take two points as it from the range @p points as input parameters
   * and return the distance between those points.
   * @param orderPolicy Order policy for the points. Can be either @ref Order_policy::FARTHEST_POINT_ORDERING,
   * @ref Order_policy::ALREADY_ORDERED or @ref Order_policy::RANDOM_POINT_ORDERING.
   *
   * @return Instianciation of @ref Oscillating_rips_edge_iterator.
   *
   * @warning Avoid copying the iterator as it is heavier than usual iterators.
   */
  template <typename PointRange, typename DistanceFunction>
  static Oscillating_rips_edge_iterator begin(
      Filtration_value nu, 
      Filtration_value mu, 
      const PointRange& points, 
      DistanceFunction&& distance,
      Order_policy orderPolicy = Order_policy::FARTHEST_POINT_ORDERING) 
  {
    return Oscillating_rips_edge_iterator(nu, mu, points, distance, orderPolicy);
  }

  /**
   * @brief Returns the begin iterator of a the range of edges based on @ref Oscillating_rips_edge_iterator.
   * Takes already computed epsilon values as input, but assumes that the points are already ordered accordingly. 
   * Assumes also that the last epsilon value is 0.
   *
   * @tparam PointRange Point range type.
   * @tparam DistanceFunction Type of the distance function.
   * @param nu Lower multiplicator.
   * @param mu Upper multiplicator.
   * @param orderedPoints Point cloud as an ordered range. The format of a point has to correspond to the input
   * format of the distance function. The order of the points should be in correspondence with the order of
   * the epsilon values.
   * @param distance Distance function. Has to take two points as it from the range @p points as input parameters
   * and return the distance between those points.
   * @param epsilonValues Epsilon values for the oscillating rips filtration. Should be in decreasing order.
   * And the last value should be 0.
   *
   * @return Instianciation of @ref Oscillating_rips_edge_iterator.
   *
   * @warning Avoid copying the iterator as it is heavier than usual iterators.
   */
  template <typename PointRange, typename DistanceFunction>
  static Oscillating_rips_edge_iterator begin(
      Filtration_value nu, 
      Filtration_value mu, 
      const PointRange& orderedPoints,
      DistanceFunction&& distance, 
      const std::vector<Filtration_value>& epsilonValues) 
  {
    return Oscillating_rips_edge_iterator(nu, mu, orderedPoints, distance, epsilonValues);
  }

  /**
   * @brief Returns the end iterator of a the range of edges based on @ref Oscillating_rips_edge_iterator.
   * 
   * @return Default instianciation of @ref Oscillating_rips_edge_iterator.
   */
  static Oscillating_rips_edge_iterator end() 
  {
    return Oscillating_rips_edge_iterator();
  }

 private:
  /**
   * @brief Default constructor. Should not be called and therfore private. Use as a ``static'' class only.
   */
  Oscillating_rips_edge_range(){};

  /**
  * @brief The two input types std::pair<int, Filtration_value> encode pairs
   * @f$(j, d(p_i,p_j))@f$ and @f$(k, d(p_i,p_k))@f$ for some fixed point @f$p_i@f$.
   * The operator() orders edges by length. By convention, if lengths are equal,
   * it orders pairs by taking the smaller vertex label between @f$j@f$ and @f$k@f$.
  */
  struct Point_distance_comp {
    bool operator()(const std::pair<int, Filtration_value>& p1, const std::pair<int, Filtration_value>& p2) const {
      {
        if (p1.second != p2.second) {
          return p1.second < p2.second;  // shorter first
        }
        return p1.first < p2.first;
      }
    }
  };

  /**
   * @brief Initialize the distance function and epsilon values. Updates also the multiplicators if
   * the edge modifier is active.
   * 
   * @tparam PointRange Point range type.
   * @tparam DistanceFunction Type of the distance function.
   * @param[out] nu Lower multiplicator.
   * @param[out] mu Upper multiplicator.
   * @param[out] epsilonValues Container for the epsilon values.
   * @param[out] distanceMatrix Container for the distance matrices.
   * @param[in] points Point cloud as a range.The format of a point has to correspond to the input format of the
   * distance function.
   * @param[in] distance Distance function. Has to take two points as it from the range @p points as input parameters
   * and return the distance between those points.
   * @param[in] orderPolicy Order policy for the points. Can be either @ref Order_policy::FARTHEST_POINT_ORDERING,
   * @ref Order_policy::ALREADY_ORDERED or @ref Order_policy::RANDOM_POINT_ORDERING.
   */
  template <typename PointRange, typename DistanceFunction>
  static void _initialize(Filtration_value& nu, Filtration_value& mu, std::vector<Filtration_value>& epsilonValues,
                          std::vector<std::vector<std::pair<int, Filtration_value> > >& distanceMatrix,
                          const PointRange& points, DistanceFunction&& distance, Order_policy orderPolicy) {
    GUDHI_CHECK((nu <= mu) && (nu >= 0), "Invalid parameters mu and nu");

    size_t n = points.size();  // number of points
    PointRange sortedPoints;
    sortedPoints.reserve(n);

    if constexpr (EdgeModifier::isActive_) {
      nu = EdgeModifier::apply_inverse_modifier(nu);
      mu = EdgeModifier::apply_inverse_modifier(mu);
    }

    // compute epsilon values
    if (orderPolicy == Order_policy::ALREADY_ORDERED) {
      sortedPoints.assign(points.begin(), points.end());
      epsilonValues = _compute_epsilon_values(sortedPoints, distance);
    } else if (orderPolicy == Order_policy::FARTHEST_POINT_ORDERING) {
      epsilonValues.reserve(n);
      Gudhi::subsampling::choose_n_farthest_points(distance, points,
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
   * @brief Compute the epsilon values for an ordered set of points, measuring the
   * sparsity of the ordering.
   *
   * @details Let \f$P = \{p_0, \ldots, p_{n-1}\f$ be the ordered set of points. Then
   * the method sets <CODE>eps_range[i]<\CODE> with the value \$f\varepsilon_i\$f,
   * defined as \f$\varpesilon_i = d_H(P_i,P)\f$, the Hausdorff between the points
   * \f$P_i= \{p_0, \ldots, p_{i}\}\f$ and the entire point cloud
   * \f$P = \{p_0, \ldots, p_{n-1}\}\f$.
   * 
   * @tparam PointRange Point range type.
   * @tparam DistanceFunction Type of the distance function.
   * @param sortedPoints Point cloud as an ordered range. The format of a point has to correspond to the input
   * format of the distance function.
   * @param distance Distance function. Has to take two points as it from the range @p points as input parameters
   * and return the distance between those points.
   *
   * @return Vector of decreasing epsilon values ending with 0.
   */
  template <typename PointRange, typename DistanceFunction>
  static std::vector<Filtration_value> _compute_epsilon_values(const PointRange& sortedPoints, DistanceFunction&& distance) {
    size_t n = sortedPoints.size();
    std::vector<Filtration_value> eps_range(n, std::numeric_limits<double>::infinity());

    // compute all \f$\varepsilon_i\f$ values, such that eps_range[i] ==
    // eps_i==d_H(P_i,P), for i=0 ... n-1:
    for (size_t i = 0; i < n; ++i) {
      // entering step i, maintain eps_range[j] = eps_j for j<i, and
      // eps_range[k] = d(p_k, P_{i-1}) for k >= i.
#ifdef GUDHI_USE_TBB
      tbb::parallel_for(size_t(i + 1), n, [&](size_t k) {
        // set eps_range[k] <- d(p_k, P_i) ==
        //                            min{ d(p_k, P_{i-1}), d(p_k, p_i) }  for k >= i.
        double dist_i_k = distance(sortedPoints[i], sortedPoints[k]);
        if (dist_i_k < eps_range[k]) {
          eps_range[k] = dist_i_k;
        }
      });
#else
      for (size_t k = i + 1; k < n; ++k) {
        // set eps_range[k] <- d(p_k, P_i) ==
        //                            min{ d(p_k, P_{i-1}), d(p_k, p_i) }  for k >= i.
        double dist_i_k = distance(points[i], points[k]);
        if (dist_i_k < eps_range[k]) {
          eps_range[k] = dist_i_k;
        }
      }
#endif
      // we have now eps_range[k] = d(p_k, P_i) for k > i.
      // to do: implement parallel version by dividing the vector
      // set eps_range[i] <- eps_i = d_H(P_i,P) = max_{k>i} d(p_k, P_i)
      double eps_i = 0.;
      for (size_t k = i + 1; k < n; ++k) {
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
      const PointRange& sortedPoints, DistanceFunction&& distance) {
    std::vector<std::vector<std::pair<int, Filtration_value> > > distanceMatrix(sortedPoints.size());
#ifdef GUDHI_USE_TBB
    tbb::parallel_for(size_t(0), sortedPoints.size(), [&](size_t i) {
      // distanceMatrix[i] = std::vector< std::pair<int, Filtration_value> >();
      distanceMatrix[i].resize(i);
      for (size_t j = 0; j < i; ++j) {
        distanceMatrix[i][j] = std::make_pair(j, distance(sortedPoints[i], sortedPoints[j]));
      }
      // distanceMatrix[i] is sorted by (j, d(p_i,p_j)) < (k, d(p_i,p_k)) iff
      // d(p_i,p_j) < d(p_i,p_k) or (j<k in case d(p_i,p_j) == d(p_i,p_k)).
      std::stable_sort(distanceMatrix[i].begin(), distanceMatrix[i].end(), Point_distance_comp());
    });
#else
    for (size_t i = 0; i < sortedPoints.size(); ++i) {  // for all vertices
      // distanceMatrix[i] = std::vector< std::pair<int, Filtration_value> >();
      distanceMatrix[i].resize(i);
      for (size_t j = 0; j < i; ++j) {
        distanceMatrix[i][j] = std::make_pair(j, distance(sortedPoints[i], sortedPoints[j]));
      }
      std::stable_sort(distanceMatrix[i].begin(), distanceMatrix[i].end(), Point_distance_comp());
    }
#endif

    return distanceMatrix;
  }

  /**
   * @brief Computes the edges that are added and the edges that are removed and stores them in two separate containers.
   * The edge container @f$e@f$ stores the @f$i^{th}@f$ edge induced by vertex @f$j@f$ at @f$e[j][i]@f$.
   * 
   * @param[in] nu Lower multiplicator.
   * @param[in] mu Upper multiplicator.
   * @param[in] epsilonValues Epsilon values in decreasing order finishing with 0.
   * @param[in] distanceMatrix Spare distance matrix.
   * @param[out] edgesAdded Container for positive edges.
   * @param[out] edgesRemoved Container for negative edges.
   *
   * @return Total number of edges.
   */
  static size_t _compute_edges(Filtration_value nu, Filtration_value mu,
                               const std::vector<Filtration_value>& epsilonValues,
                               std::vector<std::vector<std::pair<int, Filtration_value> > >& distanceMatrix,
                               std::vector<std::vector<Zigzag_edge<Filtration_value> > >& edgesAdded,
                               std::vector<std::vector<Zigzag_edge<Filtration_value> > >& edgesRemoved) 
  {
    size_t number_of_arrows = 0;
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
    // no need to consider the case i=n-1 in an oRzz filtration
    tbb::parallel_for(size_t(0), n - 1, [&](size_t i) {
      typename std::vector<std::pair<int, Filtration_value> >::iterator it;
      //----edgesAdded[i]:
      // consider first all edges added in inclusion:
      // R({p_0, ... , p_i}, nu * eps_i) -> R({p_0, ... , p_i}, mu * eps_i),
      // i.e., all (p_j,p_k) with 0 <= k < j <= i with
      //                                           nu eps_i < d(p_j,p_k) <= mu eps_i
      // these new edges get filtration value epsilonValues[i]
      for (size_t j = 1; j <= i; ++j) {
        // get very first edge (k,j), over all k<j, strictly longer than  mu * eps_i
        // distanceMatrix[j][k] = d(p_j,p_k) with k<j
        it = std::upper_bound(distanceMatrix[j].begin(), distanceMatrix[j].end(),
                              std::pair<int, Filtration_value>(n, mu * epsilonValues[i]), Point_distance_comp());

        while (it != distanceMatrix[j].begin()) {
          --it;
          // if edge already in R({p_0, ... , p_i}, nu * eps_i), stop
          if (it->second <= nu * epsilonValues[i]) {
            break;
          }
          if constexpr (EdgeModifier::isActive_){
            edgesAdded[i].emplace_back(it->first, j, EdgeModifier::apply_modifier(epsilonValues[i]), true);
          } else {
            edgesAdded[i].emplace_back(it->first, j, epsilonValues[i], true);
          }
        //   std::cout << j << ", " << (it - distanceMatrix[j].begin()) << ", " << it->first << "\n";
          ++number_of_arrows;
        }
      }
      // now consider all edges added in inclusion:
      // R({p_0, ... , p_i}, mu * eps_i) -> R({p_0, ... , p_i, p_i+1}, mu * eps_i)
      // i.e., all (p_j,p_i+1) with 0 <= j <= i with       d(p_j,p_i+1) <= mu eps_i
      // these new edges get filtration value epsilonValues[i]
      // first striclty longer edge
      it = std::upper_bound(distanceMatrix[i + 1].begin(), distanceMatrix[i + 1].end(),
                            std::pair<int, Filtration_value>(n, mu * epsilonValues[i]), Point_distance_comp());

      while (it != distanceMatrix[i + 1].begin()) {
        --it;
        if constexpr (EdgeModifier::isActive_) {
          edgesAdded[i].emplace_back(it->first, i + 1, EdgeModifier::apply_modifier(epsilonValues[i]), true);
        } else {
          edgesAdded[i].emplace_back(it->first, i + 1, epsilonValues[i], true);
        }
        // std::cout << (i + 1) << ", " << (it - distanceMatrix[i+1].begin()) << ", " << it->first << "\n";
        ++number_of_arrows;
      }

      //----edgesRemoved[i]:
      // consider all edges removed in
      // R({p_0, ... , p_{i+1}}, mu * eps_i) <- R({p_0, ... , p_{i+1}}, nu * eps_{i+1})
      // i.e., all edges (p_k,p_j), 0<=k<j<=i+1, such that
      // nu eps_{i+1} < d(p_k,p_j) <= mu eps_i
      // these new edges get filtration value epsilonValues[i]
      for (size_t j = 1; j <= i + 1; ++j) {
        // get very first edge (k,j), over all k<j, strictly longer than  mu * eps_i
        // distanceMatrix[j][k] = d(p_j,p_k) with k<j
        it = std::upper_bound(distanceMatrix[j].begin(), distanceMatrix[j].end(),
                              std::pair<int, Filtration_value>(n, mu * epsilonValues[i]), Point_distance_comp());

        while (it != distanceMatrix[j].begin()) {
          --it;
          // when reading an edge in R({p_0, ... , p_{i+1}}, nu * eps_{i+1}), stop
          if (it->second <= nu * epsilonValues[i + 1]) {
            break;
          }
          // edgesRemoved[i].emplace_back(it->first, j, epsilonValues[i+1], false);
          if constexpr (EdgeModifier::isActive_){
            edgesRemoved[i].emplace_back(it->first, j, EdgeModifier::apply_modifier(epsilonValues[i]), false);
          } else {
            edgesRemoved[i].emplace_back(it->first, j, epsilonValues[i], false);
          }
        //   std::cout << j << ", " << (it - distanceMatrix[j].begin()) << ", " << it->first << "\n";
          ++number_of_arrows;
        }
      }
    });
#else  // GUDHI_USE_TBB not defined

    typename std::vector<std::pair<int, Filtration_value> >::iterator it;

    for (size_t i = 0; i < n - 1; ++i) {
      //----edgesAdded[i]:
      // consider first all edges added in inclusion:
      // R({p_0, ... , p_i}, nu * eps_i) -> R({p_0, ... , p_i}, mu * eps_i),
      // i.e., all (p_j,p_k) with 0 <= k < j <= i with
      //                                           nu eps_i < d(p_j,p_k) <= mu eps_i
      // these new edges get filtration value epsilonValues[i]
      for (size_t j = 1; j <= i; ++j) {
        // get very first edge (k,j), over all k<j, strictly longer than  mu * eps_i
        // distanceMatrix[j][k] = d(p_j,p_k) with k<j
        it = std::upper_bound(distanceMatrix[j].begin(), distanceMatrix[j].end(),
                              std::pair<int, Filtration_value>(n, mu * epsilonValues[i]), Point_distance_comp());

        while (it != distanceMatrix[j].begin()) {
          --it;
          // if edge already in R({p_0, ... , p_i}, nu * eps_i), stop
          if (it->second <= nu * epsilonValues[i]) {
            break;
          }
          if constexpr (EdgeModifier::isActive_){
            edgesAdded[i].emplace_back(it->first, j, EdgeModifier::apply_modifier(epsilonValues[i]), true);
          } else {
            edgesAdded[i].emplace_back(it->first, j, epsilonValues[i], true);
          }
          ++number_of_arrows;
        }
      }
      // now consider all edges added in inclusion:
      // R({p_0, ... , p_i}, mu * eps_i) -> R({p_0, ... , p_i, p_i+1}, mu * eps_i)
      // i.e., all (p_j,p_i+1) with 0 <= j <= i with       d(p_j,p_i+1) <= mu eps_i
      // these new edges get filtration value epsilonValues[i]
      // first striclty longer edge
      it = std::upper_bound(distanceMatrix[i + 1].begin(), distanceMatrix[i + 1].end(),
                            std::pair<int, Filtration_value>(n, mu * epsilonValues[i]), Point_distance_comp());
// if (i == 0){
//     std::cout << "start vect colind2: " << (it - distanceMatrix[1].begin()) << ", (" << it->first << ", " << it->second << "), (" << distanceMatrix[1].begin()->first << ", " << distanceMatrix[1].begin()->second << ")\n";
// }
      while (it != distanceMatrix[i + 1].begin()) {
        --it;
        if constexpr (EdgeModifier::isActive_) {
          edgesAdded[i].emplace_back(it->first, i + 1, EdgeModifier::apply_modifier(epsilonValues[i]), true);
        } else {
          edgesAdded[i].emplace_back(it->first, i + 1, epsilonValues[i], true);
        //   if (i==0) std::cout << "added: " << it->first << ", " << (i + 1) << ", " << epsilonValues[i] << "\n";
        }
        ++number_of_arrows;
      }

      //----edgesRemoved[i]:
      // consider all edges removed in
      // R({p_0, ... , p_{i+1}}, mu * eps_i) <- R({p_0, ... , p_{i+1}}, nu * eps_{i+1})
      // i.e., all edges (p_k,p_j), 0<=k<j<=i+1, such that
      // nu eps_{i+1} < d(p_k,p_j) <= mu eps_i
      // these new edges get filtration value epsilonValues[i]
      for (size_t j = 1; j <= i + 1; ++j) {
        // get very first edge (k,j), over all k<j, strictly longer than  mu * eps_i
        // distanceMatrix[j][k] = d(p_j,p_k) with k<j
        it = std::upper_bound(distanceMatrix[j].begin(), distanceMatrix[j].end(),
                              std::pair<int, Filtration_value>(n, mu * epsilonValues[i]), Point_distance_comp());

        while (it != distanceMatrix[j].begin()) {
          --it;
          // when reading an edge in R({p_0, ... , p_{i+1}}, nu * eps_{i+1}), stop
          if (it->second <= nu * epsilonValues[i + 1]) {
            break;
          }
          // edgesRemoved[i].emplace_back(it->first, j, epsilonValues[i+1], false);
          if constexpr (EdgeModifier::isActive_){
            edgesRemoved[i].emplace_back(it->first, j, EdgeModifier::apply_modifier(epsilonValues[i]), false);
          } else {
            edgesRemoved[i].emplace_back(it->first, j, epsilonValues[i], false);
          }
          ++number_of_arrows;
        }
      }
    }
#endif

    return number_of_arrows;
  }

  /**
   * @brief Sorts canonically the edges: as much as possible, edges should be removed in
   * the reverse order of their insertion. We decide to insert shorted edges first,
   * with increasing lexicographical order, and remove larger edges first, with
   * decreasing lexicographic order.
   *
   * @param edges Edges to sort.
   */
  static void _canonically_sort_edges(std::vector<Zigzag_edge<Filtration_value> >& edges) {
    // filtration then dimension, then lex order for insertion
    auto edge_cmp = [](const Zigzag_edge<Filtration_value>& e1, const Zigzag_edge<Filtration_value>& e2) {
      if (e1.get_filtration_value() != e2.get_filtration_value()) {
        return e1.get_filtration_value() < e2.get_filtration_value();   // lower fil first
      }

      if (e1.get_smallest_vertex() == e1.get_biggest_vertex()) {    // e1 is a vertex, -> put vertices first
        if (e2.get_smallest_vertex() == e2.get_biggest_vertex()) {
          return e1.get_smallest_vertex() < e2.get_smallest_vertex();   //-> vertex of lower label
        } else {
          return true;  //-> always vertices before edges
        }
      }
      // e1 is an edge
      if (e2.get_smallest_vertex() == e2.get_biggest_vertex()) {
        return false;   // e2 vertex, -> put it first
      }
      // both are edges, lexigraphic compare
      if (e1.get_smallest_vertex() != e2.get_smallest_vertex()) {
        return e1.get_smallest_vertex() < e2.get_smallest_vertex(); // lex order
      }
      if (e1.get_biggest_vertex() != e2.get_biggest_vertex()) {
        return e1.get_biggest_vertex() < e2.get_biggest_vertex();
      }
      return false;     // equality
    };

    // the inverse ordering for deletions
    auto inv_edge_cmp = [&](const Zigzag_edge<Filtration_value>& e1, const Zigzag_edge<Filtration_value>& e2) {
      if (e1.get_smallest_vertex() == e2.get_smallest_vertex() && e1.get_biggest_vertex() == e2.get_biggest_vertex()) {
        return false;               //equality => false
      }
      return !(edge_cmp(e1, e2));   // reverse order
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
      } else {          // sequence of removals
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

/**
 * @class Oscillating_rips_simplex_range oscillating_rips_iterators.h gudhi/Zigzag_persistence/oscillating_rips_iterators.h
 * @brief Gives access to a range of the ordered simplices in an oscilliating Rips filtration.
 *
 * @ingroup zigzag_persistence
 *
 * @details The simplices are returned as tuples of three elements: 
 * the first is the simplex handle of the simplex in the given complex,
 * the second is the filtration value of the corresponding arrow,
 * the third is the direction of the arrow, i.e., indicates if the simplex is inserted or removed.
 * For more information, see @ref Oscillating_rips_iterator.
 * 
 * @tparam StableFilteredComplex Filtered complex structure that has stable simplex handles, 
 * that is, they do not invalidates after an insertion or a removal (except for the removed simplices).
 * @tparam EdgeRangeIterator Type of the edge range iterator. 
 *
 * @warning As the custom iterator stores a possibly large range, avoid copying it. 
 * Use @ref get_iterator_range, @ref begin and @ref end wisely.
 */
template <class StableFilteredComplex, typename EdgeRangeIterator>
class Oscillating_rips_simplex_range {
 public:
  /**
   * @class Oscillating_rips_iterator oscillating_rips_iterators.h gudhi/Zigzag_persistence/oscillating_rips_iterators.h
   * @brief Custom iterator over the simplices of an oscillating rips filtration. Is a forward iterator only.
   *
   * It inherits from boost::iterator_facade.
   *
   * For each positive edge given by the edge range iterator, the given complex computes and adds all possible
   * cofaces and they respective simplex handle are stored. For each negative edge, the given complex computes
   * the simplex handles of all cofaces and they are stored.
   * The simplex handles induced by an edge are stored only the time needed and are removed when reaching the next edge.
   * That is, we start by storing the possible cofaces of an edge, then at each incrementation,
   * the next element within the stored simplex handles is retrieved and once the end is reached,
   * the simplex handles are erased and replaced by the possible cofaces of the next edge and so on...
   * If the edge is negative, then a simplex is removed from the complex at the next incrementation.
   * So, the maximum size of the complex corresponds to the maximum size of a complex in the zigzag filtration.
   *
   * The simplices are returned by the iterator as tuples of three elements:
   * the first is the simplex handle of the simplex in the given complex,
   * the second is the filtration value of the corresponding arrow,
   * the third is the direction of the arrow, i.e., indicates if the simplex is inserted or removed.
   *
   * @warning As the iterator stores possibly a large range, avoid copying it.
   */
  class Oscillating_rips_iterator
      : public boost::iterator_facade<Oscillating_rips_iterator,
                                      const std::tuple<typename StableFilteredComplex::Simplex_handle,
                                                       typename StableFilteredComplex::Filtration_value, bool>&,
                                      boost::forward_traversal_tag> 
  {
   public:
    using Filtration_value = typename StableFilteredComplex::Filtration_value;  /**< Filtration value type. */
    using Simplex_handle = typename StableFilteredComplex::Simplex_handle;      /**< Simplex handle type. */
    using Simplex_key = typename StableFilteredComplex::Simplex_key;            /**< Key type. */

    /**
     * @brief Constructor. The edge iterators and the complex are not copied, so do not modifiy or invalidate 
     * them outside as long as this iterator is different from the end iterator.
     * 
     * @param edgeStartIterator Begin iterator of the oscilliating Rips edge range. 
     * Is moved, so the original iterator will probably be invalidated.
     * @param edgeEndIterator End iterator of the oscilliating Rips edge range.
     * Is moved, so the original iterator will probably be invalidated.
     * @param complex Empty complex which will be used to store current simplices. 
     * Only the address of the complex is stored, so the state of the complex can be 
     * consulted between each incrementation.
     * @param maxDimension Maximal dimension to which to expand the complex. If set to -1, there is no limit.
     * Default value: -1.
     */
    Oscillating_rips_iterator(EdgeRangeIterator& edgeStartIterator, EdgeRangeIterator& edgeEndIterator,
                              StableFilteredComplex& complex, int maxDimension = -1)
        : complex_(&complex),
          currentSimplexIndex_(0),
          currentEdgeIt_(std::move(edgeStartIterator)),     // TODO: move constructor for custom iterator?
          endEdgeIt_(std::move(edgeEndIterator)),
          currentDirection_(true),
          maxDimension_(maxDimension),
          currentArrowNumber_(0)
    {
      if (currentEdgeIt_ == endEdgeIt_) {
        _set_end();
        return;
      }

      // first simplex is the vertex (0,0) which is the only one with its filtration value.
      complex_->insert_edge_as_flag(currentEdgeIt_->get_smallest_vertex(), currentEdgeIt_->get_biggest_vertex(),
                                    currentEdgeIt_->get_filtration_value(), maxDimension_, currentSimplices_);
      ++currentEdgeIt_;

      std::get<0>(currentArrow_) = currentSimplices_[currentSimplexIndex_];
      std::get<1>(currentArrow_) = complex_->filtration(currentSimplices_[currentSimplexIndex_]);
      std::get<2>(currentArrow_) = currentDirection_;

      complex_->assign_key(currentSimplices_[currentSimplexIndex_], 0);
    }

    /**
     * @brief Default constructor. Corresponds to the end iterator.
     */
    Oscillating_rips_iterator()
        : complex_(nullptr),
          currentSimplexIndex_(0),
          currentDirection_(true),
          maxDimension_(0),
          currentArrowNumber_(0) {}

   private:
    //mandatory for the boost::iterator_facade inheritance.
    friend class boost::iterator_core_access;

    /**
     * @brief Reverse lexicographical order for the simplex handles.
     */
    struct reverse_lexicographic_order {
      explicit reverse_lexicographic_order(StableFilteredComplex* st) : st_(st) {}

      bool operator()(const Simplex_handle sh1, const Simplex_handle sh2) const {
        auto rg1 = st_->simplex_vertex_range(sh1);
        auto rg2 = st_->simplex_vertex_range(sh2);
        auto it1 = rg1.begin();
        auto it2 = rg2.begin();
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
      StableFilteredComplex* st_;
    };

    std::vector<Simplex_handle> currentSimplices_;  /**< Stores current simplex handles. */
    StableFilteredComplex* complex_;                /**< Pointer to the complex. */
    size_t currentSimplexIndex_;                    /**< Index to current position in currentSimplices_. */
    EdgeRangeIterator currentEdgeIt_;               /**< Iterator pointing to the next edge. */
    EdgeRangeIterator endEdgeIt_;                   /**< End edge iterator. */
    bool currentDirection_;                         /**< Current direction. */
    const int maxDimension_;                        /**< Maximal dimension of expansion. */
    std::tuple<Simplex_handle,Filtration_value,bool> currentArrow_;     /**< Current return element. */
    Simplex_key currentArrowNumber_;                /**< Number of incrementations. */

    /**
     * @brief Mandatory for the boost::iterator_facade inheritance. Indicates if two iterators are equal.
     * 
     * @param other Iterator to compare.
     * @return True, if the two iterators point to the same position.
     * @return False, otherwise.
     */
    bool equal(Oscillating_rips_iterator const& other) const {
      if (complex_ == nullptr) return other.complex_ == nullptr;

      return complex_ == other.complex_ && currentEdgeIt_ == other.currentEdgeIt_ &&
             currentSimplexIndex_ == other.currentSimplexIndex_;
    }

    /**
     * @brief Mandatory for the boost::iterator_facade inheritance. Dereference the iterator.
     * 
     * @return Value of the current return element.
     */
    const std::tuple<Simplex_handle, Filtration_value, bool>& dereference() const {
        return currentArrow_;
    }

    /**
     * @brief Mandatory for the boost::iterator_facade inheritance. Increments the iterator.
     */
    void increment() {
      if (!currentDirection_) complex_->remove_maximal_simplex(currentSimplices_[currentSimplexIndex_]);

      ++currentSimplexIndex_;
      ++currentArrowNumber_;

      if (currentSimplexIndex_ == currentSimplices_.size()) {
        if (currentEdgeIt_ == endEdgeIt_) {
          _set_end();
          return;
        }

        currentSimplices_.clear();

        auto fil = currentEdgeIt_->get_filtration_value();
        currentDirection_ = currentEdgeIt_->get_direction();

        std::get<1>(currentArrow_) = fil;
        std::get<2>(currentArrow_) = currentDirection_;

        if (currentDirection_) {
          _update_positive_current_simplices(fil);
        } else {
          _update_negative_current_simplices(fil);
        }
        currentSimplexIndex_ = 0;
      }

      std::get<0>(currentArrow_) = currentSimplices_[currentSimplexIndex_];
      if (currentDirection_) complex_->assign_key(currentSimplices_[currentSimplexIndex_], currentArrowNumber_);
    }

    /**
     * @brief Sets the iterator as the end iterator.
     */
    void _set_end() { complex_ = nullptr; }

    /**
     * @brief Updates @p currentSimplices_ for insertions.
     * 
     * @param fil Current filtration value.
     */
    void _update_positive_current_simplices(Filtration_value fil) {
      while (currentEdgeIt_ != endEdgeIt_ && currentEdgeIt_->get_direction() &&
             currentEdgeIt_->get_filtration_value() == fil) 
      {
        complex_->insert_edge_as_flag(currentEdgeIt_->get_smallest_vertex(), currentEdgeIt_->get_biggest_vertex(),
                                      currentEdgeIt_->get_filtration_value(), maxDimension_, currentSimplices_);
        ++currentEdgeIt_;
      }
#ifdef GUDHI_USE_TBB
      tbb::parallel_sort(currentSimplices_.begin(), currentSimplices_.end(),
                         reverse_lexicographic_order(complex_));
#else
      std::sort(currentSimplices_.begin(), currentSimplices_.end(),
                reverse_lexicographic_order(complex_));
#endif
    }

    /**
     * @brief Updates @p currentSimplices_ for removals.
     * 
     * @param fil Current filtration value.
     */
    void _update_negative_current_simplices(Filtration_value fil) {
      unsigned int count = 0;
      while (currentEdgeIt_ != endEdgeIt_ && !currentEdgeIt_->get_direction() &&
             currentEdgeIt_->get_filtration_value() == fil) 
      {
        Simplex_handle sh = complex_->find({currentEdgeIt_->get_smallest_vertex()});
        if (currentEdgeIt_->get_smallest_vertex() != currentEdgeIt_->get_biggest_vertex()) {
          sh = sh->second.children()->members().find(currentEdgeIt_->get_biggest_vertex());
        }
        auto toRemove = complex_->star_simplex_range(sh);
        currentSimplices_.insert(currentSimplices_.end(), toRemove.begin(), toRemove.end());
        ++currentEdgeIt_;
        ++count;
      }
#ifdef GUDHI_USE_TBB
      tbb::parallel_sort(
          currentSimplices_.begin(), currentSimplices_.end(),
          [&](Simplex_handle sh1, Simplex_handle sh2) -> bool { return complex_->key(sh1) > complex_->key(sh2); });
#else
      std::sort(
          currentSimplices_.begin(), currentSimplices_.end(),
          [&](Simplex_handle sh1, Simplex_handle sh2) -> bool { return complex_->key(sh1) > complex_->key(sh2); });
#endif
      if (count > 1) {  // more than 1 edge inserted: cofaces can be duplicated
        auto last = std::unique(currentSimplices_.begin(), currentSimplices_.end(),
                                [&](Simplex_handle sh1, Simplex_handle sh2) -> bool {
                                  return complex_->key(sh1) == complex_->key(sh2);
                                });                              // equal simplex handles means equal key
        currentSimplices_.erase(last, currentSimplices_.end());  // remove duplicated cofaces
      }
    }
  };

  /**
   * @brief Returns a boost::iterator_range from @ref Oscillating_rips_iterator.
   * The iterator is a forward iterator only.
   *
   * @param edgeStartIterator Begin iterator of the oscilliating Rips edge range.
   * Is moved, so the original iterator will probably be invalidated.
   * @param edgeEndIterator End iterator of the oscilliating Rips edge range.
   * Is moved, so the original iterator will probably be invalidated.
   * @param complex Empty complex which will be used to store current simplices.
   * Only the address of the complex is stored, so the state of the complex can be
   * consulted between each incrementation.
   * @param maxDimension Maximal dimension to which to expand the complex. If set to -1, there is no limit.
   * Default value: -1.
   *
   * @return boost::iterator_range of @ref Oscillating_rips_iterator.
   * 
   * @warning Avoid copying the iterators as they are heavier than usual iterators. If begin and end iterators 
   * are needed but not the structure of the range, use @ref begin and @ref end instead.
   */
  static boost::iterator_range<Oscillating_rips_iterator> get_iterator_range(EdgeRangeIterator& edgeStartIterator,
                                                                             EdgeRangeIterator& edgeEndIterator,
                                                                             StableFilteredComplex& complex,
                                                                             int maxDimension = -1) {
    return boost::iterator_range<Oscillating_rips_iterator>(
        Oscillating_rips_iterator(edgeStartIterator, edgeEndIterator, complex, maxDimension),
        Oscillating_rips_iterator());
  }

  /**
   * @brief Returns the begin iterator of a the range of simplices based on @ref Oscillating_rips_iterator.
   * 
   * @param edgeStartIterator Begin iterator of the oscilliating Rips edge range.
   * Is moved, so the original iterator will probably be invalidated.
   * @param edgeEndIterator End iterator of the oscilliating Rips edge range.
   * Is moved, so the original iterator will probably be invalidated.
   * @param complex Empty complex which will be used to store current simplices.
   * Only the address of the complex is stored, so the state of the complex can be
   * consulted between each incrementation.
   * @param maxDimension Maximal dimension to which to expand the complex. If set to -1, there is no limit.
   * Default value: -1.
   *
   * @return Instianciation of @ref Oscillating_rips_iterator.
   *
   * @warning Avoid copying the iterator as it is heavier than usual iterators.
   */
  static Oscillating_rips_iterator begin(EdgeRangeIterator& edgeStartIterator, EdgeRangeIterator& edgeEndIterator,
                                         StableFilteredComplex& complex, int maxDimension = -1) {
    return Oscillating_rips_iterator(edgeStartIterator, edgeEndIterator, complex, maxDimension);
  }

  /**
   * @brief Returns the end iterator of a the range of simplices based on @ref Oscillating_rips_iterator.
   * 
   * @return Default instianciation of @ref Oscillating_rips_iterator.
   */
  static Oscillating_rips_iterator end() { return Oscillating_rips_iterator(); }

 private:
  /**
   * @brief Default constructor. Should not be called and therfore private. Use as a ``static'' class only.
   */
  Oscillating_rips_simplex_range(){};
};

}  // namespace zigzag_persistence
}  // namespace Gudhi

#endif  // ZIGZAG_OSCILLATING_RIPS_ITERATORS_H_
