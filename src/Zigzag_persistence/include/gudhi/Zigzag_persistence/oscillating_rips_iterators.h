/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria and Hannah Schreiber
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef ZIGZAG_OSCILLATING_RIPS_ITERATORS_H_
#define ZIGZAG_OSCILLATING_RIPS_ITERATORS_H_

#include <cmath>
#include <cstddef>
#include <vector>

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
  Zigzag_edge(size_t u, size_t v, Filtration_value fil, bool direction)
      : u_(u), v_(v), fil_(fil), direction_(direction) {
        if (u > v) std::swap(u_, v_);
      }

  Zigzag_edge() : u_(0), v_(0), fil_(0), direction_(true) {}

  /* Returns vertex with smaller label. */
  size_t get_smallest_vertex() const { return u_; }
  /* Returns vertex with bigger label. */
  size_t get_biggest_vertex() const { return v_; }
  /* Returns the filtration value of the edge. */
  Filtration_value get_filtration_value() const { return fil_; }
  /* Returns true if insertion of the edge, false if removal. */
  bool get_direction() const { return direction_; }

  void set(size_t u, size_t v, Filtration_value fil, bool direction){
    u_ = u;
    v_ = v;
    fil_ = fil;
    direction_ = direction;
  }

  bool operator==(const Zigzag_edge& e) const {
    return ((e.u_ == u_) && (e.v_ == v_) && (e.fil_ == fil_) && (e.direction_ == direction_));
  }

  void assign_filtration(Filtration_value fil) { fil_ = fil; }

 private:
  size_t u_;
  size_t v_;
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

template <typename Filtration_value, class EdgeModifier = Identity_edge_modifier<Filtration_value> >
class Oscillating_rips_edge_range {
 private:
  class Oscillating_rips_edge_iterator;

 public:
  enum Order_policy { ALREADY_ORDERED, FARTHEST_POINT_ORDERING, RANDOM_POINT_ORDERING };

  template <typename PointRange, typename DistanceFunction>
  static std::vector<Zigzag_edge<Filtration_value> > compute_oscillating_rips_edges(Filtration_value nu,
                                                                                    Filtration_value mu,
                                                                                    const PointRange& points,
                                                                                    DistanceFunction&& distance,
                                                                                    Order_policy orderPolicy) 
  {
    std::vector<Zigzag_edge<Filtration_value> > edgeFiltration;
    std::vector<Filtration_value> epsilonValues;
    std::vector<std::vector<std::pair<int, Filtration_value> > > distanceMatrix;
    auto n = points.size();

    initialize_(nu, mu, epsilonValues, distanceMatrix, points, distance, orderPolicy);

    // edgesAdded[i] (resp. edgesRemoved[i]) == list of edges (i,j), with j<i, added (resp. removed) at eps_i
    // we also put there (later) vertices that are added. Note that vertices are removed
    // only at the very last step of the oRzz filtration.
    std::vector<std::vector<Zigzag_edge<Filtration_value> > > edgesAdded, edgesRemoved;

    // auto it = std::upper_bound(distanceMatrix[1].begin(), distanceMatrix[1].end(),
    //                              std::pair<int, Filtration_value>(distanceMatrix.size(), mu * epsilonValues[0]),
    //                              Point_distance_comp());
    // std::cout << "start vect colind: " << (it - distanceMatrix[1].begin()) << ", (" << it->first << ", " << it->second << "), (" << distanceMatrix[1].begin()->first << ", " << distanceMatrix[1].begin()->second << ")\n";
    size_t number_of_arrows = compute_edges_(nu, mu, epsilonValues, distanceMatrix, edgesAdded, edgesRemoved);

    // Now, sort edges according to lengths, and put everything in edgeFiltration
    edgeFiltration.clear();
    edgeFiltration.reserve(number_of_arrows + n);  // count edges + vertices additions

    // initialise R({p_0}, +infinity)
    edgeFiltration.emplace_back(0, 0,  // add vertex p_0,+infty
                                std::numeric_limits<Filtration_value>::infinity(), true);
    // epsilonValues[0], true);

    if constexpr (EdgeModifier::isActive_) {
      for (size_t i = 0; i < n - 1; ++i) {                                  // all ascending arrows eps_i
        edgeFiltration.emplace_back(i + 1, i + 1, EdgeModifier::apply_modifier(epsilonValues[i]), true);  // add p_{i+1},eps_i
        for (auto edg_it = edgesAdded[i].begin(); edg_it != edgesAdded[i].end(); ++edg_it) {
          edgeFiltration.push_back(*edg_it);
        }
        for (auto edg_it = edgesRemoved[i].rbegin();  // longest first
             edg_it != edgesRemoved[i].rend(); ++edg_it) {
          edgeFiltration.push_back(*edg_it);
        }
      }
    } else {
      for (size_t i = 0; i < n - 1; ++i) {                                  // all ascending arrows eps_i
        edgeFiltration.emplace_back(i + 1, i + 1, epsilonValues[i], true);  // add p_{i+1},eps_i
        for (auto edg_it = edgesAdded[i].begin(); edg_it != edgesAdded[i].end(); ++edg_it) {
          edgeFiltration.push_back(*edg_it);
        }
        for (auto edg_it = edgesRemoved[i].rbegin();  // longest first
             edg_it != edgesRemoved[i].rend(); ++edg_it) {
          edgeFiltration.push_back(*edg_it);
        }
      }
    }

    // what remains is removed in the zigzag iterator with -infinity values. If eps_n-1
    //== 0, which is the usual case, the remaining simplices in the filtration are
    // the n vertices.
    // cannot inforce this here.

    return edgeFiltration;
  }

  template <typename PointRange, typename DistanceFunction>
  static boost::iterator_range<Oscillating_rips_edge_iterator> compute_oscillating_rips_edges_as_iterator(
      Filtration_value nu, 
      Filtration_value mu, 
      const PointRange& points, 
      DistanceFunction&& distance,
      Order_policy orderPolicy) 
  {
    auto start = Oscillating_rips_edge_iterator(nu, mu, points, distance, orderPolicy);
    auto end = Oscillating_rips_edge_iterator();
    return boost::iterator_range<Oscillating_rips_edge_iterator>(
        start, end);
  }

 private:
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
        : nu_(nu), mu_(mu), currentEdge_(0, 0, std::numeric_limits<Filtration_value>::infinity(), true), isEnd_(false), epsValIndex_(0), rowIndex_(1), inPosEd_(true), insertVertex_(true)
    {
      initialize_(nu_, mu_, epsilonValues_, distanceMatrix_, points, distance, orderPolicy);
      auto it = std::upper_bound(distanceMatrix_[1].begin(), distanceMatrix_[1].end(),
                                 std::pair<int, Filtration_value>(distanceMatrix_.size(), mu_ * epsilonValues_[epsValIndex_]),
                                 Point_distance_comp());
      columnIndex_ = it - distanceMatrix_[1].begin();
    //   std::cout << "start colind: " << columnIndex_ << ", (" << it->first << ", " << it->second << "), (" << distanceMatrix_[1].begin()->first << ", " << distanceMatrix_[1].begin()->second << ")\n";
    }

    template <typename PointRange, typename DistanceFunction>
    Oscillating_rips_edge_iterator(Filtration_value nu, Filtration_value mu, const PointRange& orderedPoints,
                                   DistanceFunction&& distance, const std::vector<Filtration_value>& epsilonValues)
        : epsilonValues_(epsilonValues),
          nu_(nu),
          mu_(mu),
          currentEdge_(0, 0, std::numeric_limits<Filtration_value>::infinity(), true),
          isEnd_(false), epsValIndex_(0), rowIndex_(1), inPosEd_(true), insertVertex_(true) {
      GUDHI_CHECK(orderedPoints.size() == epsilonValues.size(),
                  "The number of points and the number of epsilon values should match.");
      GUDHI_CHECK((nu <= mu) && (nu >= 0), "Invalid parameters mu and nu");

      if constexpr (EdgeModifier::isActive_) {
        nu_ = EdgeModifier::apply_inverse_modifier(nu);
        mu_ = EdgeModifier::apply_inverse_modifier(mu);
      }

      // compute the distance matrix
      distanceMatrix_ = compute_distance_matrix_(orderedPoints, distance);

      auto it = std::upper_bound(distanceMatrix_[1].begin(), distanceMatrix_[1].end(),
                                 std::pair<int, Filtration_value>(distanceMatrix_.size(), mu_ * epsilonValues_[epsValIndex_]),
                                 Point_distance_comp());
      columnIndex_ = it - distanceMatrix_[1].begin();
    }

    Oscillating_rips_edge_iterator() : nu_(0), mu_(0), currentEdge_(0, 0, 0, true), isEnd_(true),  epsValIndex_(0), rowIndex_(1), columnIndex_(0), inPosEd_(true), insertVertex_(true) {}

   private:
    friend class boost::iterator_core_access;

    std::vector<Filtration_value> epsilonValues_;
    std::vector<std::vector<std::pair<int, Filtration_value> > > distanceMatrix_;
    Filtration_value nu_;
    Filtration_value mu_;
    Zigzag_edge<Filtration_value> currentEdge_;
    bool isEnd_;    //to replace
    size_t epsValIndex_, rowIndex_, columnIndex_;
    bool inPosEd_, insertVertex_;   //to replace

    bool equal(Oscillating_rips_edge_iterator const& other) const { return isEnd_ == other.isEnd_ && currentEdge_ == other.currentEdge_; }

    const Zigzag_edge<Filtration_value>& dereference() const { return currentEdge_; }

    void increment() {
      if (epsValIndex_ < distanceMatrix_.size() - 1) {
        if (insertVertex_) {
          currentEdge_.set(epsValIndex_ + 1, epsValIndex_ + 1, epsilonValues_[epsValIndex_], true);
          insertVertex_ = false;
          return;
        }

        while ((epsValIndex_ < distanceMatrix_.size() - 1 && inPosEd_ && col_index_is_not_valid_pos()) ||
               (epsValIndex_ < distanceMatrix_.size() - 1 && !inPosEd_ && col_index_is_not_valid_neg())) {
          if (inPosEd_) {
            ini_col_index_pos();
          } else {
            ini_col_index_neg();
            if (epsValIndex_ == distanceMatrix_.size() - 1) {
                set_end_();
                return;
            }
            if (insertVertex_) {
              currentEdge_.set(epsValIndex_ + 1, epsValIndex_ + 1, epsilonValues_[epsValIndex_], true);
              insertVertex_ = false;
              return;
            }
          }
        }

        if (inPosEd_) {
          up_edge_pos(epsValIndex_);
          if (rowIndex_ == epsValIndex_ + 2) {
            inPosEd_ = false;
            --rowIndex_;
            ini_col_index_low();
          }
          return;
        }

        up_edge_neg(epsValIndex_);
        if (rowIndex_ == 0) {
          ++epsValIndex_;
          if (epsValIndex_ == distanceMatrix_.size() - 1) return;
          insertVertex_ = true;
          inPosEd_ = true;
          ++rowIndex_;
          ini_col_index_up();
        }
        return;
      }

      set_end_();
    }

    void set_end_(){
        isEnd_ = true;
        currentEdge_ = Zigzag_edge<Filtration_value>();
    }

    void ini_col_index_up() {
      auto it =
          std::upper_bound(distanceMatrix_[rowIndex_].begin(), distanceMatrix_[rowIndex_].end(),
                           std::pair<int, Filtration_value>(distanceMatrix_.size(), mu_ * epsilonValues_[epsValIndex_]),
                           Point_distance_comp());
      columnIndex_ = it - distanceMatrix_[rowIndex_].begin();
    }

    void ini_col_index_low() {
      auto it =
          std::lower_bound(distanceMatrix_[rowIndex_].begin(), distanceMatrix_[rowIndex_].end(),
                           std::pair<int, Filtration_value>(0, nu_ * epsilonValues_[epsValIndex_ + 1]),
                           Point_distance_comp());
      while (it != distanceMatrix_[rowIndex_].end() && it->second == nu_ * epsilonValues_[epsValIndex_ + 1]) ++it;
      columnIndex_ = it - distanceMatrix_[rowIndex_].begin();
    }

    void ini_col_index_pos(){
        while (col_index_is_not_valid_pos() && rowIndex_ <= epsValIndex_) {
            ++rowIndex_;
            ini_col_index_up();
        }
        if (col_index_is_not_valid_pos()){
            ini_col_index_low();
            inPosEd_ = false;
        }
    }

    void ini_col_index_neg(){
        while (col_index_is_not_valid_neg() && rowIndex_ > 1) {
            --rowIndex_;
            ini_col_index_low();
        }
        if (col_index_is_not_valid_neg()){
            ini_col_index_up();
            ++epsValIndex_;
            insertVertex_ = true;
            inPosEd_ = true;
        }
    }

    bool col_index_is_not_valid_pos() {
        return columnIndex_ == 0 ||
               (rowIndex_ != (epsValIndex_ + 1) &&
                distanceMatrix_[rowIndex_][columnIndex_ - 1].second <= nu_ * epsilonValues_[epsValIndex_]);
    }

    bool col_index_is_not_valid_neg() {
        return columnIndex_ == distanceMatrix_[rowIndex_].size() ||
               distanceMatrix_[rowIndex_][columnIndex_].second > mu_ * epsilonValues_[epsValIndex_];
    }

    void up_edge_pos(size_t i) {
        --columnIndex_;
        if constexpr (EdgeModifier::isActive_)
            currentEdge_.set(distanceMatrix_[rowIndex_][columnIndex_].first, rowIndex_, EdgeModifier::apply_modifier(epsilonValues_[i]), true);
        else 
            currentEdge_.set(distanceMatrix_[rowIndex_][columnIndex_].first, rowIndex_, epsilonValues_[i], true);
        // std::cout << "pos: " << rowIndex_ << ", " << columnIndex_ << ", " << distanceMatrix_[rowIndex_][columnIndex_].first << "\n";
        while (col_index_is_not_valid_pos() && rowIndex_ <= i) {
            ++rowIndex_;
            ini_col_index_up();
        }
        if (col_index_is_not_valid_pos()) {
            ++rowIndex_;
        }
    }

    void up_edge_neg(size_t i) {
        if constexpr (EdgeModifier::isActive_)
            currentEdge_.set(distanceMatrix_[rowIndex_][columnIndex_].first, rowIndex_, EdgeModifier::apply_modifier(epsilonValues_[i]), false);
        else 
            currentEdge_.set(distanceMatrix_[rowIndex_][columnIndex_].first, rowIndex_, epsilonValues_[i], false);
        // std::cout << "neg: " << rowIndex_ << ", " << columnIndex_ << ", " << distanceMatrix_[rowIndex_][columnIndex_].first << "\n";
        ++columnIndex_;
        while (col_index_is_not_valid_neg() && rowIndex_ > 1) {
            --rowIndex_;
            ini_col_index_low();
        }
        if (col_index_is_not_valid_neg()) {
            --rowIndex_;
        }
    }
  };

  template <typename PointRange, typename DistanceFunction>
  static void initialize_(Filtration_value& nu, Filtration_value& mu, std::vector<Filtration_value>& epsilonValues,
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
      epsilonValues = compute_epsilon_values_(sortedPoints, distance);
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
      epsilonValues = compute_epsilon_values_(sortedPoints, distance);
    }

    // compute the distance matrix
    distanceMatrix = compute_distance_matrix_(sortedPoints, distance);
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
  static std::vector<Filtration_value> compute_epsilon_values_(const PointRange& points, DistanceFunction&& distance) {
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
  static std::vector<std::vector<std::pair<int, Filtration_value> > > compute_distance_matrix_(
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

  static size_t compute_edges_(Filtration_value nu, Filtration_value mu,
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
};

}  // namespace zigzag_persistence
}  // namespace Gudhi

#endif  // ZIGZAG_OSCILLATING_RIPS_ITERATORS_H_
