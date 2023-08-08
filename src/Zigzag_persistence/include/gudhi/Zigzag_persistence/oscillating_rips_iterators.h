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

template <typename Filtration_value>
class Zigzag_edge {
 public:
  Zigzag_edge(int u, int v, Filtration_value fil, bool direction)
      : u_(u), v_(v), fil_(fil), direction_(direction) {
        if (u > v) std::swap(u_, v_);
      }

  Zigzag_edge() : u_(0), v_(0), fil_(0), direction_(true) {}

  /* Returns vertex with smaller label. */
  int get_smallest_vertex() const { return u_; }
  /* Returns vertex with bigger label. */
  int get_biggest_vertex() const { return v_; }
  /* Returns the filtration value of the edge. */
  Filtration_value get_filtration_value() const { return fil_; }
  /* Returns true if insertion of the edge, false if removal. */
  bool get_direction() const { return direction_; }

  void set(int u, int v, Filtration_value fil, bool direction){
    u_ = u;
    v_ = v;
    fil_ = fil;
    direction_ = direction;
  }

  bool operator==(const Zigzag_edge& e) const {
    return ((e.u_ == u_) && (e.v_ == v_) && (e.fil_ == fil_) && (e.direction_ == direction_));
  }

//   bool operator<(const Zigzag_edge& e) const {
//     if (e.fil_ != fil_) return fil_ < e.fil_;
//     if (e.direction_ != direction_) return direction_;
//     if (e.u_ != u_) return u_ < e.u_;
//     return v_ < e.v_;
//   }

 private:
  int u_;
  int v_;
  Filtration_value fil_;
  bool direction_;
};

template <typename Filtration_value>
class Identity_edge_modifier {
 public:
  static constexpr bool isActive_ = false;

 private:
  Identity_edge_modifier() {}
};

template <typename Filtration_value>
class Square_root_edge_modifier {
 public:
  static Filtration_value apply_modifier(Filtration_value f) { return std::sqrt(f); }
  static Filtration_value apply_inverse_modifier(Filtration_value f) { return f * f; }

  static constexpr bool isActive_ = true;

 private:
  Square_root_edge_modifier() {}
};

//assumes that eps_n-1 == 0
template <typename Filtration_value, class EdgeModifier = Identity_edge_modifier<Filtration_value> >
class Oscillating_rips_edge_range {
 public:
  enum Order_policy { ALREADY_ORDERED, FARTHEST_POINT_ORDERING, RANDOM_POINT_ORDERING };

  class Oscillating_rips_edge_iterator
      : public boost::iterator_facade<Oscillating_rips_edge_iterator, 
                                      const Zigzag_edge<Filtration_value>&,
                                      boost::forward_traversal_tag> {
   public:
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
    friend class boost::iterator_core_access;

    std::vector<Filtration_value> epsilonValues_;
    std::vector<std::vector<std::pair<int, Filtration_value> > > distanceMatrix_;
    Filtration_value nu_;
    Filtration_value mu_;
    Zigzag_edge<Filtration_value> currentEdge_;
    size_t epsilonIndex_, rowIndex_, columnIndex_;
    bool inPositiveDirection_, insertVertex_;

    bool equal(Oscillating_rips_edge_iterator const& other) const {
      return rowIndex_ == other.rowIndex_ && currentEdge_ == other.currentEdge_;
    }

    const Zigzag_edge<Filtration_value>& dereference() const { return currentEdge_; }

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

    void _set_end(){
        rowIndex_ = 0;
        currentEdge_ = Zigzag_edge<Filtration_value>();
    }

    void _initialize_positive_col_index() {
      auto it =
          std::upper_bound(distanceMatrix_[rowIndex_].begin(), distanceMatrix_[rowIndex_].end(),
                           std::pair<int, Filtration_value>(distanceMatrix_.size(), mu_ * epsilonValues_[epsilonIndex_]),
                           Point_distance_comp());
      columnIndex_ = it - distanceMatrix_[rowIndex_].begin();
    }

    void _initialize_negative_col_index() {
      auto it =
          std::lower_bound(distanceMatrix_[rowIndex_].begin(), distanceMatrix_[rowIndex_].end(),
                           std::pair<int, Filtration_value>(0, nu_ * epsilonValues_[epsilonIndex_ + 1]),
                           Point_distance_comp());
      while (it != distanceMatrix_[rowIndex_].end() && it->second == nu_ * epsilonValues_[epsilonIndex_ + 1]) ++it;
      columnIndex_ = it - distanceMatrix_[rowIndex_].begin();
    }

    bool _positive_col_index_is_not_valid() {
        return columnIndex_ == 0 ||
               (rowIndex_ != (epsilonIndex_ + 1) &&
                distanceMatrix_[rowIndex_][columnIndex_ - 1].second <= nu_ * epsilonValues_[epsilonIndex_]);
    }

    bool _negative_col_index_is_not_valid() {
        return columnIndex_ == distanceMatrix_[rowIndex_].size() ||
               distanceMatrix_[rowIndex_][columnIndex_].second > mu_ * epsilonValues_[epsilonIndex_];
    }

    void _update_edge(size_t i, bool direction) {
        if constexpr (EdgeModifier::isActive_)
            currentEdge_.set(distanceMatrix_[rowIndex_][columnIndex_].first, rowIndex_,
                             EdgeModifier::apply_modifier(epsilonValues_[i]), direction);
        else
            currentEdge_.set(distanceMatrix_[rowIndex_][columnIndex_].first, rowIndex_, epsilonValues_[i], direction);
    }

    void _update_edge_as_positive_vertex() {
        if constexpr (EdgeModifier::isActive_)
            currentEdge_.set(epsilonIndex_ + 1, epsilonIndex_ + 1,
                             EdgeModifier::apply_modifier(epsilonValues_[epsilonIndex_]), true);
        else
            currentEdge_.set(epsilonIndex_ + 1, epsilonIndex_ + 1, epsilonValues_[epsilonIndex_], true);
    }

    void _update_edge_as_negative_vertex() {
        currentEdge_.set(rowIndex_ - 1, rowIndex_ - 1, -std::numeric_limits<Filtration_value>::infinity(), false);
    }
  };

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

  //as Oscillating_rips_edge_iterator is a heavy iterator to copy, it should not be used as a usual iterator
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

  static Oscillating_rips_edge_iterator end() 
  {
    return Oscillating_rips_edge_iterator();
  }

 private:
  Oscillating_rips_edge_range(){};

  /* The two input types std::pair<int, Filtration_value> encode pairs
   * (j, d(p_i,p_j)) and (k, d(p_i,p_k)) for some fixed point p_i.
   * The operator() orders edges by length. By convention, if lengths are equal,
   * it orders pairs by taking the smaller vertex label between j and k.
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

  /** \brief Compute the epsilon values for an ordered set of points, measuring the
   * sparsity of the ordering.
   *
   * \details Let \f$P = \{p_0, \ldots, p_{n-1}\f$ be the ordered set of points. Then
   * the method sets <CODE>eps_range[i]<\CODE> with the value \$f\varepsilon_i\$f,
   * defined as \f$\varpesilon_i = d_H(P_i,P)\f$, the Hausdorff between the points
   * \f$P_i= \{p_0, \ldots, p_{i}\}\f$ and the entire point cloud
   * \f$P = \{p_0, \ldots, p_{n-1}\}\f$.
   *
   * @param[in] points Range of points.
   * @param[in] distance Distance function that can be called on two points.
   * @param[in] eps_range   Vector in which the epsilon values are written,
   *                        <CODE>eps_range[i]<\CODE> is \$f\varepsilon_i\$f.
   *                        Satisfyies <CODE>epsilonValues[i] >=
   *                        epsilonValues[j] >= 0<\CODE> whenever <CODE>i>=j<\CODE>.
   *                        The range must be of same size as the number of points.
   */
  template <typename PointRange, typename DistanceFunction>
  static std::vector<Filtration_value> _compute_epsilon_values(const PointRange& points, DistanceFunction&& distance) {
    size_t n = points.size();
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
        double dist_i_k = distance(points[i], points[k]);
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

  static void _canonically_sort_edges(std::vector<Zigzag_edge<Filtration_value> >& edges) {
    // canonical sort of the edges: as much as possible, edges should be removed in
    // the reverse order of their insertion. We decide to insert shorted edges first,
    // with increasing lexicographical order, and remove larger edges first, with
    // decreasing lexicographic order.

    // filtration then dimension, then lex order for insertion
    auto edge_cmp = [](const Zigzag_edge<Filtration_value>& e1, const Zigzag_edge<Filtration_value>& e2) {
      if (e1.get_filtration_value() != e2.get_filtration_value()) {
        return e1.get_filtration_value() < e2.get_filtration_value();
      }  // lower fil first

      if (e1.get_smallest_vertex() == e1.get_biggest_vertex()) {  // e1 is a vertex, -> put vertices first
        if (e2.get_smallest_vertex() == e2.get_biggest_vertex()) {
          return e1.get_smallest_vertex() < e2.get_smallest_vertex();
        }  //-> vertex of lower label
        else {
          return true;
        }  //-> always vertices before edges
      }
      // e1 is an edge
      if (e2.get_smallest_vertex() == e2.get_biggest_vertex()) {
        return false;
      }       // e2 vertex, -> put it first
      // both are edges, lexigraphic compare
      if (e1.get_smallest_vertex() != e2.get_smallest_vertex()) {
        return e1.get_smallest_vertex() < e2.get_smallest_vertex();
      }  // lex order
      if (e1.get_biggest_vertex() != e2.get_biggest_vertex()) {
        return e1.get_biggest_vertex() < e2.get_biggest_vertex();
      }
      return false;  // equality
    };
    // the inverse ordering for deletions
    auto inv_edge_cmp = [&](const Zigzag_edge<Filtration_value>& e1, const Zigzag_edge<Filtration_value>& e2) {
      if (e1.get_smallest_vertex() == e2.get_smallest_vertex() && e1.get_biggest_vertex() == e2.get_biggest_vertex()) {
        return false;
      }                            //== => false
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
      if (curr_type) {
#ifdef GUDHI_USE_TBB
        tbb::parallel_sort(beg, end, edge_cmp);
#else
        std::sort(beg, end, edge_cmp);
#endif
      }  // sequence of insertions
      else {
#ifdef GUDHI_USE_TBB
        tbb::parallel_sort(beg, end, inv_edge_cmp);
#else
        std::sort(beg, end, inv_edge_cmp);
#endif
      }  // sequence of removals
      beg = end;
      curr_fil = beg->get_filtration_value();
      curr_type = beg->get_direction();
    }
  }
};

//needs Filtered_complex to have stable simplex handles
template <class Filtered_complex, typename EdgeRangeIterator>
class Oscillating_rips_simplex_range {
 public:
  class Oscillating_rips_iterator
      : public boost::iterator_facade<Oscillating_rips_iterator,
                                      const std::tuple<typename Filtered_complex::Simplex_handle,typename Filtered_complex::Filtration_value,bool>&,
                                      boost::forward_traversal_tag> {
   public:
    using Filtration_value = typename Filtered_complex::Filtration_value;
    using Simplex_handle = typename Filtered_complex::Simplex_handle;
    using Simplex_key = typename Filtered_complex::Simplex_key;

    // edges and complex are not copied, so do not modifiy outside as long as iterator != end.
    // TODO: move constructor for iterator?
    Oscillating_rips_iterator(EdgeRangeIterator& edgeStartIterator, EdgeRangeIterator& edgeEndIterator,
                              Filtered_complex& complex, int maxDimension = -1)
        : complex_(&complex),
          currentSimplexIndex_(0),
          currentEdgeIt_(std::move(edgeStartIterator)),
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

    Oscillating_rips_iterator()
        : complex_(nullptr),
          currentSimplexIndex_(0),
          currentDirection_(true),
          maxDimension_(0),
          currentArrowNumber_(0) {}

   private:
    friend class boost::iterator_core_access;

    struct reverse_lexicographic_order {
      explicit reverse_lexicographic_order(Filtered_complex* st) : st_(st) {}

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
      Filtered_complex* st_;
    };

    std::vector<Simplex_handle> currentSimplices_;
    Filtered_complex* complex_;
    size_t currentSimplexIndex_;
    EdgeRangeIterator currentEdgeIt_;
    EdgeRangeIterator endEdgeIt_;
    bool currentDirection_;
    const int maxDimension_;
    //   bool isEnd_;
    std::tuple<Simplex_handle, Filtration_value, bool> currentArrow_;
    Simplex_key currentArrowNumber_;

    bool equal(Oscillating_rips_iterator const& other) const {
      if (complex_ == nullptr) return other.complex_ == nullptr;

      return complex_ == other.complex_ && currentEdgeIt_ == other.currentEdgeIt_ &&
             currentSimplexIndex_ == other.currentSimplexIndex_;
    }

    const std::tuple<Simplex_handle, Filtration_value, bool>& dereference() const {
        return currentArrow_;
    }

    void increment() {
      if (!currentDirection_) complex_->remove_maximal_simplex(currentSimplices_[currentSimplexIndex_]);

      ++currentSimplexIndex_;
      ++currentArrowNumber_;

      if (currentSimplexIndex_ == currentSimplices_.size()) {
        if (currentEdgeIt_ == endEdgeIt_) {
          _set_end();
          return;
        }

        // if (!currentDirection_) {
        //   for (auto& sh : currentSimplices_) complex_->remove_maximal_simplex(sh);
        // }
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

    void _set_end() { complex_ = nullptr; }

    void _update_positive_current_simplices(Filtration_value fil) {
      while (currentEdgeIt_ != endEdgeIt_ && currentEdgeIt_->get_direction() &&
             currentEdgeIt_->get_filtration_value() == fil) 
      {
        complex_->insert_edge_as_flag(currentEdgeIt_->get_smallest_vertex(), currentEdgeIt_->get_biggest_vertex(),
                                      currentEdgeIt_->get_filtration_value(), maxDimension_, currentSimplices_);
        // std::cout << "add edge: " << currentEdgeIt_->get_smallest_vertex() << " " << currentEdgeIt_->get_biggest_vertex() << " - " << currentEdgeIt_->get_filtration_value() << "\n";
        // std::cout << "current:\n";
        // for (auto sh : currentSimplices_){
        //     for (auto v : complex_->simplex_vertex_range(sh)) std::cout << v << " ";
        //     std::cout << "\n";
        // }
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

  static boost::iterator_range<Oscillating_rips_iterator> get_iterator_range(EdgeRangeIterator& edgeStartIterator,
                                                                             EdgeRangeIterator& edgeEndIterator,
                                                                             Filtered_complex& complex,
                                                                             int maxDimension = -1) {
    return boost::iterator_range<Oscillating_rips_iterator>(
        Oscillating_rips_iterator(edgeStartIterator, edgeEndIterator, complex, maxDimension),
        Oscillating_rips_iterator());
  }

  static Oscillating_rips_iterator begin(EdgeRangeIterator& edgeStartIterator, EdgeRangeIterator& edgeEndIterator,
                                         Filtered_complex& complex, int maxDimension = -1) {
    return Oscillating_rips_iterator(edgeStartIterator, edgeEndIterator, complex, maxDimension);
  }

  static Oscillating_rips_iterator end() { return Oscillating_rips_iterator(); }

 private:
  Oscillating_rips_simplex_range(){};
};

}  // namespace zigzag_persistence
}  // namespace Gudhi

#endif  // ZIGZAG_OSCILLATING_RIPS_ITERATORS_H_
