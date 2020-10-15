/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - YYYY/MM author: description of the modification.
 *         
 */

#ifndef SIMPLEX_TREE_ZIGZAG_ITERATORS_H_
#define SIMPLEX_TREE_ZIGZAG_ITERATORS_H_

#include <iostream>
#include <fstream>
#include <gudhi/choose_n_farthest_points.h>
#include <gudhi/pick_n_random_points.h>
#include <gudhi/Discrete_morse_theory.h>

#ifdef GUDHI_USE_TBB
#include <tbb/tbb.h>
#endif

/** Represents an edge for a zigzag filtrations of flag complexes. 
  * The edge must have two endpoints, encoded by Vertex_handles, a filtration 
  * value and a type (insertion or deletion) represented by a bool.
  *
  * A sequence of such edges represents a full flag zigzag filtration.
  *
  * An edge with u_ == v_ represents a vertex labeled u_.
  */
template< class FilteredComplex >
class Zigzag_edge {
public:
  Zigzag_edge( typename FilteredComplex::Vertex_handle u
             , typename FilteredComplex::Vertex_handle v
             , typename FilteredComplex::Filtration_value fil
             , bool type)
  : u_(u), v_(v), fil_(fil), type_(type) {}

  Zigzag_edge() : u_(-1), v_(-1), fil_(-1) {}

/* Returns vertex with smaller label. */
  typename FilteredComplex::Vertex_handle    u()    const { return u_; }
/* Returns vertex with bigger label. */
  typename FilteredComplex::Vertex_handle    v() const { return v_; }
/* Returns the filtration value of the edge. */
  typename FilteredComplex::Filtration_value fil()  { return fil_; }
/* Returns true if insertion of the edge, false if removal. */
  bool                                       type() { return type_; }

  bool operator==(const Zigzag_edge &e) const {
    return ( (e.u_ == u_) && (e.v_ == v_) && 
             (e.fil_ == fil_) && (e.type_ == type_) );
  }

  void assign_filtration(typename FilteredComplex::Filtration_value fil) { fil_ = fil; }

private:
  typename FilteredComplex::Vertex_handle    u_;
  typename FilteredComplex::Vertex_handle    v_;
  typename FilteredComplex::Filtration_value fil_;
  bool                                       type_;
};

/** \brief Computes the 1-skeleton-filtration allowing to generate the oscillating 
  * Rips zigzag filtration of a discrete metric space (such as a set of points).
  *
  * Let \f$P = \{p_0 ... p_{n-1}\}\f$ be an ordered, finite, set of points in a 
  * metric space with distance \f$d(p_i,p_j)\f$. 
  * Let \f$(\varepsilon_0, \ldots , \varepsilon_{n-1})\f$ be a decreasing sequence 
  * of positive scalars, 
  * i.e., such that \f$ \varepsilon_i \geq \varepsilon_j \geq 0 \f$ for any 
  * \f$n>i>j\geq 0\f$.
  * Let \f$\mu\f$ and \f$\nu\f$ be two scalars such that \f$0 < \nu \leq \mu\f$.
  *
  * Specifically, denote by \f$P_i\f$ the subset of points \f$\{p_0 ... p_i\}\f$ 
  * for \f$0\leq i \leq n-1\f$, and \f$\mathcal{R}^{\rho}(P_i)\f$ the Rips complex 
  * of threshold \f$\rho \geq 0\f$ built on points \f$P_i\f$. 
  *
  * The oscillating Rips zigzag filtration associated to this data is the zigzag 
  * filtration is the filtration:
  *
  * \f$
  *            \mathcal{R}^{+\infty}(P_0) \ \ \ \ \ 
  *   \supseteq \mathcal{R}^{\nu \cdot \varepsilon_0}(P_0) 
  *    \subseteq \mathcal{R}^{\mu \cdot \varepsilon_0}(P_1) 
  *     \supseteq \mathcal{R}^{\nu \cdot \varepsilon_1}(P_1) 
  *      \subseteq \mathcal{R}^{\mu \cdot \varepsilon_1}(P_2) 
  *       \supseteq \mathcal{R}^{\nu \cdot \varepsilon_2}(P_2) 
  *        \subseteq \cdots 
  *         \supseteq \mathcal{R}^{\nu \cdot \varepsilon_{n-1}}(P_{n-1})
  *          \ \ \ \ \ \supseteq \mathcal{R}^{-\infty}(\emptyset) \f$
  *
  * By convention, we start with a Rips complex on a single vertex p_0, 
  * with \f$\varpesilon_{-1}=+\infty\f$ threshold, and finish with an empty Rips 
  * complex with \f$\varpesilon_{n}=-\infty\f$ threshold. These extra epsilon 
  * values at indices -1 and n are in particular used in zigzag persistent homology.
  * More specifically:
  *
  * A simplex inserted in the inclusion 
  * \f$\mathcal{R}^{\nu \cdot \varepsilon_i}(P_i) 
  *     \subseteq \mathcal{R}^{\mu \cdot \varepsilon_i}(P_{i+1})\f$ 
  * is given filtration value \f$\varepsilon_i\f$. The vertex $p_0$ inserted as 
  * first simplex is given filtration value \f$\varepsilon_{-1} = +\infty\f$. 
  *
  * A simplex removed in the inclusion 
  * \f$\mathcal{R}^{\mu \cdot \varepsilon_i}(P_{i+1}) 
  *     \supseteq \mathcal{R}^{\nu \cdot \varepsilon_{i+1}}(P_{i+1})\f$ 
  * is given filtration value \f$\varepsilon_{i+1}\f$. The vertices removed when 
  * emptying the complex in the rightmost arrow are given filtration value 
  * \f$\varpesilon_{n}=-\infty\f$.
  *
  * The method computes the sequence of insertions and removals 
  * of the oscillating Rips zigzag filtration restricted to vertices and edges, 
  * i.e., the 1-skeleton. 
  * Because all complexes of the filtration are flag complexes, the 
  * 1-skeleton-filtration encodes fully the entire filtration, in a much more 
  * compact way.
  *
  * Edge_t must be of type Zigzag_edge.
  * 
  * @param[in] dist_matrix The distance matrix representing the discrete metric 
  *   space with ordered points \f$P = \{p_0 ... p_{n-1}\}\f$. dist_mat[i][j] must 
  *   return the value \f$d(p_i,p_j)\f$ whenever \f$i>j\f$. 
  * @param[in] eps_values Range of epsilon values, such that eps_values[i] >= 
  *   eps_values[j] >= 0 whenever i>j. The range must be of same size as the number 
  *   of points in the discrete metric space.
  * @param[in] mu, nu The two parameters \f$\mu,\nu\f$ such that
  *   \f$\mu\geq\nu\geq0\f$
  * @param[in] edge_filtration Vector in which the 1-skeleton filtration is 
  *   written. The vector is cleared beforehand. 
  */
template< typename DistanceMatrix,
          typename FiltrationValueRange,
          typename FiltrationValue,
          typename Edge_t >
void zigzag_filtration_one_skeleton(DistanceMatrix            & dist_mat,
                                    FiltrationValueRange      & eps_values,
                                    FiltrationValue const       nu,
                                    FiltrationValue const       mu,
                                    std::vector<Edge_t>       & edge_filtration )
{ 
  GUDHI_CHECK(dist_mat.size() == eps_values.size(), "invalid number of points and epsilon values");
  GUDHI_CHECK((nu <= mu) && (nu >= 0), "invalid parameters mu and nu");
  bool decreasing_eps = true;
  FiltrationValue prev_eps = std::numeric_limits<FiltrationValue>::infinity();
  for(auto eps : eps_values) {
    if(eps > prev_eps) { decreasing_eps = false; break; }
    else { prev_eps = eps; }
  }
  GUDHI_CHECK(decreasing_eps, "non-decreasing sequence of epsilon values");

  std::cout << "D\n";

  //initialize a distance matrix where dist_matrix[i][j] containing the pair 
  //(j, d(p_i,p_j)) for j < i. We sort edges according to length later.
  size_t n = dist_mat.size();//number of points
  std::vector< std::vector< std::pair<int, FiltrationValue> > > dist_matrix;
  dist_matrix.resize(n);
  //total number of insertion and removal of vertices and edges in the zigzag 
  //filtration
  size_t number_of_arrows = 0;

 /* The two input types std::pair<int, FiltrationValue> encode pairs 
  * (j, d(p_i,p_j)) and (k, d(p_i,p_k)) for some fixed point p_i. 
  * The operator() orders edges by length. By convention, if lengths are equal, 
  * it orders pairs by taking the smaller vertex label between j and k.
  */
  struct point_distance_cmp {
    bool operator()( std::pair<int, FiltrationValue> p1
                   , std::pair<int, FiltrationValue> p2 ) {
      { 
        if(p1.second != p2.second) {return p1.second < p2.second;} //shorter first
        return p1.first < p2.first; 
      }
    }
  };
  point_distance_cmp cmp;

  std::cout << "E\n";

#ifdef GUDHI_USE_TBB
  tbb::parallel_for(size_t(0), n, [&](size_t i) {
    // dist_matrix[i] = std::vector< std::pair<int, FiltrationValue> >();
    dist_matrix[i].resize(i);
    for(size_t j=0; j<i; ++j) {
      dist_matrix[i][j] = std::make_pair(j, dist_mat[i][j]);
                 // = std::pair<int, FiltrationValue>(j, dist_mat[i][j]);
    }
  //dist_matrix[i] is sorted by (j, d(p_i,p_j)) < (k, d(p_i,p_k)) iff 
  //d(p_i,p_j) < d(p_i,p_k) or (j<k in case d(p_i,p_j) == d(p_i,p_k)).
    std::stable_sort(dist_matrix[i].begin(), dist_matrix[i].end(), cmp);
  } );
#else
  for(size_t i=0; i<n; ++i) {//for all vertices
    // dist_matrix[i] = std::vector< std::pair<int, FiltrationValue> >();
    dist_matrix[i].resize(i);
    for(size_t j=0; j<i; ++j) {
      dist_matrix[i][j] = std::make_pair(j, dist_mat[i][j]);
                // = std::pair<int, FiltrationValue>(j, dist_mat[i][j]);
    }
    std::stable_sort(dist_matrix[i].begin(), dist_matrix[i].end(), cmp);
  }
#endif

//edges_added[i] (resp. edges_removed[i]) == list of edges (i,j), with j<i, added (resp. removed) at eps_i
//we also put there (later) vertices that are added. Note that vertices are removed 
//only at the very last step of the oRzz filtration.
  std::vector< std::vector< Edge_t > > edges_added, edges_removed; 
  edges_added.resize(n);  edges_removed.resize(n);

  //edges_added[i] must contain all new edges (k,j), 0 <= k < j <= i+1, 
  //inserted in inclusion: 
  //R({p_0, ... , p_i}, nu * eps_i) -> R({p_0, ... , p_i, p_i+1 }, mu * eps_i)
  //and
  //edges_removed[i] must contain all edges (k,j), 0 <= k < j <= i+1, 
  //removed in backward inclusion:
  //R({p_0, ... , p_{i+1}}, mu * eps_i) <- R({p_0, ... , p_{i+1}}, nu * eps_{i+1})
#ifdef GUDHI_USE_TBB
  //no need to consider the case i=n-1 in an oRzz filtration
  tbb::parallel_for(size_t(0), n-1, [&](size_t i) {
    typename std::vector< std::pair<int, FiltrationValue> >::iterator it;
  //----edges_added[i]:
    //consider first all edges added in inclusion:
    //R({p_0, ... , p_i}, nu * eps_i) -> R({p_0, ... , p_i}, mu * eps_i),
    //i.e., all (p_j,p_k) with 0 <= k < j <= i with
    //                                           nu eps_i < d(p_j,p_k) <= mu eps_i  
    for(size_t j = 1; j <= i ; ++j) { 
      //get very first edge (k,j), over all k<j, strictly longer than  mu * eps_i
      //dist_matrix[j][k] = d(p_j,p_k) with k<j
      it = std::upper_bound(dist_matrix[j].begin(), dist_matrix[j].end()
                   , std::pair<int, FiltrationValue>(n, mu * eps_values[i])
                   , cmp);

      while(it != dist_matrix[j].begin()) {
        --it;
        //if edge already in R({p_0, ... , p_i}, nu * eps_i), stop
        if(it->second <= nu * eps_values[i]) { break; }
        edges_added[i].emplace_back(it->first, j, eps_values[i], true);
        ++number_of_arrows;
      }
    }
    //now consider all edges added in inclusion:
    //R({p_0, ... , p_i}, mu * eps_i) -> R({p_0, ... , p_i, p_i+1}, mu * eps_i)
    //i.e., all (p_j,p_i+1) with 0 <= j <= i with       d(p_j,p_i+1) <= mu eps_i  
    //first striclty longer edge
    it = std::upper_bound(dist_matrix[i+1].begin(), dist_matrix[i+1].end(), 
            std::pair<int, FiltrationValue>(n, mu * eps_values[i]), cmp); 

    while(it != dist_matrix[i+1].begin()) {
      --it;
      edges_added[i].emplace_back(it->first, i+1, eps_values[i], true);
      ++number_of_arrows;
    }

  //----edges_removed[i]:
    //consider all edges removed in
    //R({p_0, ... , p_{i+1}}, mu * eps_i) <- R({p_0, ... , p_{i+1}}, nu * eps_{i+1})
    //i.e., all edges (p_k,p_j), 0<=k<j<=i+1, such that 
    // nu eps_{i+1} < d(p_k,p_j) <= mu eps_i
    for(size_t j = 1; j <= i+1; ++j) {
      //get very first edge (k,j), over all k<j, strictly longer than  mu * eps_i
      //dist_matrix[j][k] = d(p_j,p_k) with k<j
      it = std::upper_bound(dist_matrix[j].begin(), dist_matrix[j].end(), 
             std::pair<int, FiltrationValue>(n, mu * eps_values[i]), cmp ); 

      while(it != dist_matrix[j].begin()) {
        --it;
        //when reading an edge in R({p_0, ... , p_{i+1}}, nu * eps_{i+1}), stop
        if(it->second <= nu * eps_values[i+1]) { break; }
        edges_removed[i].emplace_back(it->first, j, eps_values[i+1], false);
        ++number_of_arrows;
      }
    }
  } );
#else //GUDHI_USE_TBB not defined

  typename std::vector< std::pair<int, FiltrationValue> >::iterator it;

  for(size_t i=0; i<n-1; ++i) {
  //----edges_added[i]:
    //consider first all edges added in inclusion:
    //R({p_0, ... , p_i}, nu * eps_i) -> R({p_0, ... , p_i}, mu * eps_i),
    //i.e., all (p_j,p_k) with 0 <= k < j <= i with
    //                                           nu eps_i < d(p_j,p_k) <= mu eps_i  
    for(size_t j = 1; j <= i ; ++j) { 
      //get very first edge (k,j), over all k<j, strictly longer than  mu * eps_i
      //dist_matrix[j][k] = d(p_j,p_k) with k<j
      it = std::upper_bound(dist_matrix[j].begin(), dist_matrix[j].end()
                   , std::pair<int, FiltrationValue>(n, mu * eps_values[i])
                   , cmp);

      while(it != dist_matrix[j].begin()) {
        --it;
        //if edge already in R({p_0, ... , p_i}, nu * eps_i), stop
        if(it->second <= nu * eps_values[i]) { break; }
        edges_added[i].emplace_back(it->first, j, eps_values[i], true);
        ++number_of_arrows;
      }
    }
    //now consider all edges added in inclusion:
    //R({p_0, ... , p_i}, mu * eps_i) -> R({p_0, ... , p_i, p_i+1}, mu * eps_i)
    //i.e., all (p_j,p_i+1) with 0 <= j <= i with       d(p_j,p_i+1) <= mu eps_i  
    //first striclty longer edge
    it = std::upper_bound(dist_matrix[i+1].begin(), dist_matrix[i+1].end(), 
            std::pair<int, FiltrationValue>(n, mu * eps_values[i]), cmp); 

    while(it != dist_matrix[i+1].begin()) {
      --it;
      edges_added[i].emplace_back(it->first, i+1, eps_values[i], true);
      ++number_of_arrows;
    }

  //----edges_removed[i]:
    //consider all edges removed in
    //R({p_0, ... , p_{i+1}}, mu * eps_i) <- R({p_0, ... , p_{i+1}}, nu * eps_{i+1})
    //i.e., all edges (p_k,p_j), 0<=k<j<=i+1, such that 
    // nu eps_{i+1} < d(p_k,p_j) <= mu eps_i
    for(size_t j = 1; j <= i+1; ++j) {
      //get very first edge (k,j), over all k<j, strictly longer than  mu * eps_i
      //dist_matrix[j][k] = d(p_j,p_k) with k<j
      it = std::upper_bound(dist_matrix[j].begin(), dist_matrix[j].end(), 
             std::pair<int, FiltrationValue>(n, mu * eps_values[i]), cmp ); 

      while(it != dist_matrix[j].begin()) {
        --it;
        //when reading an edge in R({p_0, ... , p_{i+1}}, nu * eps_{i+1}), stop
        if(it->second <= nu * eps_values[i+1]) { break; }
        edges_removed[i].emplace_back(it->first, j, eps_values[i+1], false);
        ++number_of_arrows;
      }
    }
  }  
#endif

  //Now, sort edges according to lengths, and put everything in edge_filtration
  edge_filtration.clear();
  edge_filtration.reserve(number_of_arrows + n); //count edges + vertices additions

  // Compare edges by distance first (shorter is smaller), and lexicographic 
  // ordering 
  // otherwise. This agrees with the usual Rips filtration in standard persistence.
  struct edge_cmp {
    edge_cmp(DistanceMatrix *dist_mat) 
    : dist_ptr_(dist_mat) {}

    bool operator()(Edge_t e1, Edge_t e2) 
    { //lengths of edges e1 and e2
      int u1 = e1.u(); int v1 = e1.v(); if(u1>v1) { std::swap(u1,v1); }
      int u2 = e2.u(); int v2 = e2.v(); if(u2>v2) { std::swap(u2,v2); }
      FiltrationValue dist1 = (*dist_ptr_)[v1][u1];
      FiltrationValue dist2 = (*dist_ptr_)[v2][u2];;
      if(dist1  != dist2)  {return dist1 < dist2;}
      if(u1 != u2) {return u1 < u2;}
      return v1 < v2; 
    }

    DistanceMatrix * dist_ptr_;
  };

//sort insertions and deletion by edge length, then lexicographic order
  edge_cmp e_cmp(&dist_mat);

#ifdef GUDHI_USE_TBB
  tbb::parallel_for(size_t(0), n-1, [&](size_t i) {
    //add shortest edges first
    std::stable_sort(edges_added[i].begin(), edges_added[i].end(), e_cmp);
    //remove longest edges first (read from right to left), see below
    std::stable_sort(edges_removed[i].begin(), edges_removed[i].end(), e_cmp);
  } );
#else
  for(size_t i = 0; i < n-1; ++i) {
    //add shortest edges first
    std::stable_sort(edges_added[i].begin(), edges_added[i].end(), e_cmp);
    //remove longest edges first (read from right to left), see below
    std::stable_sort(edges_removed[i].begin(), edges_removed[i].end(), e_cmp);
  }
#endif
 
  std::cout << "H\n";


  //initialise R({p_0}, +infinity)
  edge_filtration.emplace_back(0, 0, //add vertex p_0,+infty
                          std::numeric_limits<FiltrationValue>::infinity(), true);
  for(size_t i = 0; i < n-1; ++i) {//all ascending arrows eps_i
    edge_filtration.emplace_back(i+1, i+1, eps_values[i], true);//add p_{i+1},eps_i
    for(auto edg_it = edges_added[i].begin(); 
             edg_it != edges_added[i].end(); ++edg_it) {
      edge_filtration.push_back(*edg_it);
    }
    for(auto edg_it = edges_removed[i].rbegin(); //longest first
             edg_it != edges_removed[i].rend(); ++edg_it) {
      edge_filtration.push_back(*edg_it);
    }
  }

  std::cout << "I\n";


  //what remains is removed in the zigzag iterator with -infinity values. If eps_n-1
  //== 0, which is the usual case, the remaining simplices in the filtration are 
  //the n vertices.
}



/** Compute the epsilon values for an ordered set of points 
  * \f$P = \{p_0, \ldots, p_{n-1}\f$, with  
  * eps_range[i] <- eps_i = d_H(P_i,P), where d_H is the Hausdorff distance and 
  * P_i = \{p_0, \ldots, p_{n-1}\}.
  */
template< typename PointRange,
          typename Distance, //furnish()
          typename EpsilonRange >
void compute_epsilon_values(PointRange const & points,//ordered set of points
                            Distance           distance,
                            EpsilonRange     & eps_range) 
{
  size_t n = points.size();
  eps_range.resize(n,std::numeric_limits<double>::infinity());

  //compute all \f$\varepsilon_i\f$ values, such that eps_range[i] == 
  //eps_i==d_H(P_i,P), for i=0 ... n-1:
  for(size_t i=0; i<n; ++i) {
    //entering step i, maintain eps_range[j] = eps_j for j<i, and
    //eps_range[k] = d(p_k, P_{i-1}) for k >= i.
#ifdef GUDHI_USE_TBB
    tbb::parallel_for(size_t(i+1), n, [&](size_t k) {
    //set eps_range[k] <- d(p_k, P_i) == 
    //                           min{ d(p_k, P_{i-1}), d(p_k, p_i) }  for k >= i.
      double dist_i_k = distance(points[i],points[k]); 
      if(dist_i_k < eps_range[k]) { eps_range[k] = dist_i_k; }
    } );
#else
    for(size_t k=i+1; k<n; ++k) {
    //set eps_range[k] <- d(p_k, P_i) == 
    //                           min{ d(p_k, P_{i-1}), d(p_k, p_i) }  for k >= i.
      double dist_i_k = distance(points[i],points[k]); 
      if(dist_i_k < eps_range[k]) { eps_range[k] = dist_i_k; }
    }
#endif
  //we have now eps_range[k] = d(p_k, P_i) for k > i.
  //to do: implement parallel version by dividing the vector
      //set eps_range[i] <- eps_i = d_H(P_i,P) = max_{k>i} d(p_k, P_i) 
    double eps_i = 0.;
    for(size_t k=i+1; k<n; ++k) {
      if(eps_range[k] > eps_i) { eps_i = eps_range[k]; }
    }
    eps_range[i] = eps_i;
  }
}

/** Re-ordering policies for a set of input points.
  *
  * already_ordered implies no re-ordering, O(1) complexity.
  *
  * farthest_point_ordering re-orders points with a farthest point strategy, 
  *   starting with point [0]. This strategy is the slowest but garantees the best 
  *   sparsification of the point set.
  *
  * random_ordering shuffles randomly the point set.
  */
struct already_ordered { int policy(){return 0;} };
struct farthest_point_ordering {int policy(){return 1;} };
struct random_point_ordering {int policy(){return 2;} };

/** \brief Given a set of points \f$P = \{p_0 ... p_{n-1}\}\f$, 
  * computes the 1-skeleton-filtration corresponding to the oscillating Rips zigzag 
  * filtration.
  *
  * The method first orders the points according to an order policy, with value 
  * already_ordered, farthest_point_ordering, or random_ordering, and computes the 
  * values eps_i = d_H(P_i,P) in 
  * eps_values[i] = eps_i, where P_i consists of the first i points after ordering 
  * P, and d_H denotes the Hausdorff distance. A simplex appearing in the 
  * inclusion  
  * R({p_0, ... , p_{i}}, nu * eps_i) -> R({p_1, ... , p_i, p_{i+1}}, mu * eps_i)
  * is given filtration value eps_i.
  *
  * @param[in] points A set of points in a metric space.
  * @param[in] distance The distance function of the metric space.
  * @param[in] mu, nu Parameters of the zigzag Rips filtration
  * @param[in] edge_filtration Where the 1-skeleton filtration is written.
  */
template< typename PointRange,
          typename Distance, //furnish()
          typename FiltrationValue,
          typename Edge_t,
          typename OrderPolicy >
void zigzag_filtration_one_skeleton(PointRange      const  & points,
                						        Distance                 distance,
                						        FiltrationValue          nu,
                						        FiltrationValue          mu,
                						        std::vector<Edge_t>    & edge_filtration,
                                    OrderPolicy   order_policy = already_ordered() ) 
{
  std::vector<FiltrationValue> eps_values;
  size_t n = points.size();//number of points
  PointRange sorted_points; sorted_points.reserve(n);
  edge_filtration.clear();

  switch(order_policy.policy()) {
  //points already ordered, compute the epsilon values in eps_values  
    case 0:
      sorted_points.assign(points.begin(),points.end());
      compute_epsilon_values(sorted_points, distance, eps_values);
      break;
  //order points according to a farthest point ordering, starting with point 0.
    case 1: 
      eps_values.reserve(n); 
      Gudhi::subsampling::choose_n_farthest_points( 
          distance
        , points.begin()
        , points.end()
        , n //num of input points
        , n //order all points
        , 0//start with point [0]//Gudhi::subsampling::random_starting_point
        , std::back_inserter(sorted_points)
        , std::back_inserter(eps_values) ); 
      //need to shift values output by subsampling:
      for(int i=1; i<n; ++i) { eps_values[i-1] = eps_values[i]; }
      eps_values[n-1] = 0;
      break;  
  //order points randomly.
    case 2: 
      Gudhi::subsampling::pick_n_random_points(points, n, 
                                               std::back_inserter(sorted_points) );
      compute_epsilon_values(sorted_points, distance, eps_values);
      break;
    default:
      std::cerr << "Non valid ordering policy in zigzag_filtration_one_skeleton.\n";
  }

//compute the distance matrix
  std::vector< std::vector<FiltrationValue> > dist_mat(n);
#ifdef GUDHI_USE_TBB
  tbb::parallel_for(size_t(0), n, [&](size_t i) {
    // dist_matrix[i] = std::vector< std::pair<int, FiltrationValue> >();
    dist_mat[i].resize(i);
    for(size_t j=0; j<i; ++j) {
      dist_mat[i][j] = distance(sorted_points[i],sorted_points[j]);
    } 
  } );
#else
  for(size_t i=0; i<n; ++i) {//for all vertices
    // dist_mat[i] = std::vector< std::pair<int, FiltrationValue> >();
    dist_mat[i].resize(i);
    for(size_t j=0; j<i; ++j) {
      dist_mat[i][j] = distance(sorted_points[i],sorted_points[j]);
    } 
  }
#endif
//compute the 1-skeleton filtration on the sorted set of points
  zigzag_filtration_one_skeleton(dist_mat, eps_values, nu, mu, edge_filtration);
}

/** 
  * Iterator over a flag zigzag filtration implicitly 
  * represented by a list of Zigzag_edges. 
  *
  * Given an empty FlagZigzagFiltrationComplex and a range of insertions and 
  * deletion of vertices and edges, the iterator adds/removes on the fly the range 
  * of edges and 
  * expand the complex in consequence to maintain the d-skeleton of the induced 
  * flag complexes filtration. It traverses all the newly added/removed 
  * simplices (induced by the new edge) before doing further modifications.
  */
template< class FlagZigzagFilteredComplex >
class Flagzigzag_simplex_iterator 
: public boost::iterator_facade<
            Flagzigzag_simplex_iterator<FlagZigzagFilteredComplex>
          , typename FlagZigzagFilteredComplex::Simplex_handle
          , boost::forward_traversal_tag >
{
  public:
  typedef typename FlagZigzagFilteredComplex::Simplex_handle   Simplex_handle;
  // typedef typename FlagZigzagFilteredComplex::Edge_type        Edge_type;
  typedef Zigzag_edge<FlagZigzagFilteredComplex>               Edge_type;
  typedef typename FlagZigzagFilteredComplex::Filtration_value Filtration_value;

/** Default constructor also encodes the end() iterator for any 
  * Flagzigzag_simplex_range.
  */
  Flagzigzag_simplex_iterator() 
  : cpx_(nullptr) //checking for end() <=> checking for cpx_ == nullptr
  , zigzag_edge_filtration_()
  , dim_max_(-1) 
  , partial_zzfil_()
  , arrow_direction_(true)
  , counter_insert(0)
  , are_we_done(true)
  , fil_(0)
  {
    sh_it_ = partial_zzfil_.begin(); edge_it_ = zigzag_edge_filtration_.begin();
  }

/** Constructor from a point cloud and a distance function. 
  * Constructs the d-skeleton of the oscillating Rips zigzag filtration on the 
  * ordered set of points with paramters mu and nu.
  */
  template< typename PointRange,
            typename Distance,
            typename OrderPolicy >
  Flagzigzag_simplex_iterator(FlagZigzagFilteredComplex    * cpx,
                              Filtration_value const         nu,
                              Filtration_value const         mu,
                              int                            dim_max,
                              PointRange       const       & points,
                              Distance         const         distance,
                              OrderPolicy order_policy = farthest_point_ordering()) 
  {
    //check that the model of FlagZigzagFilteredComplex is empty
    GUDHI_CHECK(cpx->empty(), "complex must be empty");

    //compute the filtration of the 1-skeleton for the oRzz filtration
    zigzag_edge_filtration_ = std::vector< Edge_type >();

    std::cout << "B\n";

    std::vector<Filtration_value> filtration_values;
    zigzag_filtration_one_skeleton(points, distance, nu, mu,
                                   zigzag_edge_filtration_, 
                                   order_policy );

    std::cout << "C\n";


    dim_max_                = dim_max;
    are_we_done             = false;
    cpx_                    = cpx;
    counter_insert          = 0;
    partial_zzfil_          = std::vector< Simplex_handle >(); //TODO?
    edge_it_                = zigzag_edge_filtration_.begin();

    //if 1-skeleton filtration is empty, set the iterator to end() by cpx_<-nullptr    
    if(edge_it_ == zigzag_edge_filtration_.end()) { cpx_ = nullptr; return; }
    //otherwise, add the first edge and expand the Rips complex
    arrow_direction_ = edge_it_->type(); //must be true, i.e., an insertion
    GUDHI_CHECK(arrow_direction_, "cannot remove a simplex from an empty complex");
    
    cpx_->flag_add_edge( edge_it_->u(), edge_it_->v(), edge_it_->fil()
                       , dim_max_, partial_zzfil_);
    fil_ = edge_it_->fil();
    sh_it_ = partial_zzfil_.begin();
    ++edge_it_;
    for(auto & sh : partial_zzfil_) 
    { cpx_->assign_key(sh,counter_insert); ++counter_insert; } 
  }

/** Constructor from a pre-computed 1-skeleton zigzag filtration.
  */
  Flagzigzag_simplex_iterator( FlagZigzagFilteredComplex * cpx 
                             , std::vector< Edge_type >  & zz_edge_fil
                             , int                         dim_max )
  {
    GUDHI_CHECK(cpx->empty(), "complex must be empty");
    // zigzag_edge_filtration_ = std::vector< Edge_type >();
    zigzag_edge_filtration_ = zz_edge_fil;
    dim_max_                = dim_max;
    are_we_done             = false;
    cpx_                    = cpx;
    counter_insert          = 0;
    partial_zzfil_          = std::vector< Simplex_handle >(); //TODO?
    edge_it_                = zigzag_edge_filtration_.begin();
    if(edge_it_ == zigzag_edge_filtration_.end()) 
    { cpx_ = nullptr; return; } //end() iterator
    
    //add the first edge
    arrow_direction_ = edge_it_->type(); //must be true, i.e., an insertion
    GUDHI_CHECK(arrow_direction_, "cannot remove a simplex from an empty complex");
    cpx_->flag_add_edge( edge_it_->u(), edge_it_->v(), edge_it_->fil()
                       , dim_max_, partial_zzfil_);
    fil_ = edge_it_->fil();
    sh_it_ = partial_zzfil_.begin();
    ++edge_it_;
    for(auto & sh : partial_zzfil_) 
    { cpx_->assign_key(sh,counter_insert); ++counter_insert; } 
  }

  // Flagzigzag_simplex_iterator(const Flagzigzag_simplex_iterator& other) = delete;

  //because the iterator modifies a complex represented by pointer, the iterator 
  //must be non-copiable.
  // Flagzigzag_simplex_iterator(const Flagzigzag_simplex_iterator & ) = delete;
  //move constructor
  Flagzigzag_simplex_iterator(Flagzigzag_simplex_iterator&& other) 
  : cpx_(other.cpx_)
  , dim_max_(other.dim_max_)
  , sh_it_(other.begin())
  , edge_it_(other.edge_it_)
  , arrow_direction_(other.arrow_direction_)
  , counter_insert(other.counter_insert)
  , are_we_done(other.are_we_done)
  , fil_(other.fil_) 
  { 
    zigzag_edge_filtration_.clear();
    zigzag_edge_filtration_ = std::move(other.zigzag_edge_filtration_);
    partial_zzfil_.clear();
    partial_zzfil_ = std::move(other.partial_zzfil_);
  }
/** Copy constructor can only be called if no increment has been called on the 
  * pointer (i.e., right after initialization).
  */
  Flagzigzag_simplex_iterator(const Flagzigzag_simplex_iterator& other) 
  { 
    GUDHI_CHECK( (other.sh_it_ == other.partial_zzfil_.begin()), 
            "Flagzigzag_simplex_iterator copy constructor - Unsafe copy" );

    if(other.zigzag_edge_filtration_.empty()) {//case edge filtration empty
      GUDHI_CHECK(edge_it_ == other.zigzag_edge_filtration_.begin(), 
            "Flagzigzag_simplex_iterator copy constructor - Unsafe copy" );
    }
    else {//general case
      GUDHI_CHECK((other.edge_it_ == ++(other.zigzag_edge_filtration_.begin()) ), 
              "Flagzigzag_simplex_iterator copy constructor - Unsafe copy" );
    }

    cpx_                    = other.cpx_;
    dim_max_                = other.dim_max_;
    arrow_direction_        = other.arrow_direction_;
    counter_insert          = other.counter_insert;
    are_we_done             = other.are_we_done;
    fil_                    = other.fil_; 
    zigzag_edge_filtration_ = other.zigzag_edge_filtration_;
    edge_it_                = zigzag_edge_filtration_.begin();
    if(!zigzag_edge_filtration_.empty()) { ++edge_it_; }//next edge
    partial_zzfil_          = other.partial_zzfil_;
    sh_it_                  = partial_zzfil_.begin();
  }

  // Flagzigzag_simplex_iterator& operator=(const Flagzigzag_simplex_iterator& ) = 
                                                                            // delete;



//move assignement
  Flagzigzag_simplex_iterator& operator=(Flagzigzag_simplex_iterator&& other )
  {
    // std::cout << "call to copy assignement =??\n";
    cpx_                    = other.cpx_;
    zigzag_edge_filtration_.clear();
    zigzag_edge_filtration_ = std::move(other.zigzag_edge_filtration_);
    dim_max_                = other.dim_max_;
    partial_zzfil_.clear();
    partial_zzfil_          = std::move(other.partial_zzfil_);
    sh_it_                  = other.sh_it_;
    edge_it_                = other.edge_it_;
    arrow_direction_        = other.arrow_direction_;
    counter_insert          = other.counter_insert;
    are_we_done             = other.are_we_done;
    fil_ = other.fil_;
    return *this;
  }

// private://prevent copy constructor and assignmenet to be called
  // Flagzigzag_simplex_iterator(const Flagzigzag_simplex_iterator & );
  // Flagzigzag_simplex_iterator& operator=(const Flagzigzag_simplex_iterator& );
  //User-defined copy constructor DANGEROUS
    // Flagzigzag_simplex_iterator(const Flagzigzag_simplex_iterator& other )
    // , cpx_(other.cpx_)
    // , zigzag_edge_filtration_(other.zigzag_edge_filtration_)
    // , dim_max_(other.dim_max_)
    // , partial_zzfil_(other.partial_zzfil_)
    // , sh_it_(partial_zzfil_.begin())
    // , edge_it_(other.edge_it_)
    // , arrow_direction_(other.arrow_direction_)
    // , counter_insert(other.counter_insert)
    // , are_we_done(other.are_we_done)
    // , fil_(other.fil_) {} 

/** Returns true if the Simplex_handle pointed to is an insertion, false if it is a
  * deletion.*/
  bool arrow_direction() { return arrow_direction_; }
/** Returns the filtration value of the simplex pointed to.*/
  Filtration_value filtration() { return fil_; } 
/** Returns an upper bound on the dimension of the simplices of the flag complexes.
  * Flag complexes are expanded up to dimension dim_max().
  */
  int dim_max() { return dim_max_; }

private:
  friend class boost::iterator_core_access;

  bool equal(Flagzigzag_simplex_iterator const& other) const {
    if(cpx_ == nullptr) { return (other.cpx_ == nullptr); }      
    return ( cpx_     == other.cpx_     && 
             edge_it_ == other.edge_it_ &&
             sh_it_   == other.sh_it_ );
  }

  Simplex_handle & dereference() const { return *sh_it_; }

  void increment() 
  {
    ++sh_it_;
    if(sh_it_ == partial_zzfil_.end()) //add or remove the next edge
    { //check if we have reached the end of a sequence of backward arrows, 
      //associated to the removal of an edge. If so, we remove effectively 
      //the simplices from the complex.
      if(!arrow_direction_) //need to effectively remove the simpl. just considered
      { //The simplices in partial_zzfil_ come by decreasing keys, hence they
        //are all maximal when removing from left to right (it's a filtration 
        //read in reverse).  
        //FlagZigzagFilteredComplex::Dictionary must not invalidate iterators
        //when inserting and removing (e.g., std::map<,>).

        //effectively remove all simplices from partial_zzfil_; must be sorted 
        for(auto sh : partial_zzfil_) { cpx_->remove_maximal_simplex(sh); } 
        counter_insert += partial_zzfil_.size();
      }
      partial_zzfil_.clear();//empty the chunk of filtration
      //if all edges have been considered: no edge left to remove, 
      //there may still be simplices in the complex
      if(edge_it_ == zigzag_edge_filtration_.end()) 
      { 
        if(are_we_done) { //no more edges, no more simplices, we are done
          //explicitly empty the complex
          for(auto sh : partial_zzfil_) { cpx_->remove_maximal_simplex(sh); } 
          counter_insert += partial_zzfil_.size();
          partial_zzfil_.clear();
          cpx_ = nullptr; //set iterator to end() 
          return; 
        } 
        else {//no edge left, remove the simplices remaining in the complex 
          fil_ = - std::numeric_limits<Filtration_value>::infinity();
          are_we_done = true;//happens once
          //fills up zz_partial with the remaining simplices in complex, but 
          //does not actually remove them from the complex.
          //comes sorted by decreasing key values
          cpx_->flag_lazy_empty_complex(partial_zzfil_); 
          arrow_direction_ = false; //only backward arrows now, these are removals
          sh_it_ = partial_zzfil_.begin();
          return;
        }
      }
      //partial_zzfil_ is empty and edge_it points to a new edge
      if( edge_it_->type() ) { //forward arrow //modify the complex
        //add the new edge and expand the flag complex. partial_zz_fil_ points to
        //all the newly inserted simplices, in filtration order.
        cpx_->flag_add_edge( edge_it_->u(), edge_it_->v()
                           , edge_it_->fil()
                           , dim_max_, partial_zzfil_ );
        arrow_direction_ = true; //the arrow is forward, these are insertions

        //flag_add_edge output a SORTED sequence of simplices
        for(auto & sh : partial_zzfil_) //set key values
        { cpx_->assign_key(sh,counter_insert); ++counter_insert; }
      }
      else { //backward arrow
        //record all simplices to remove, due to the removal of an edge, 
        //but do not actually remove them from the complex.
        cpx_->flag_lazy_remove_edge( edge_it_->u(), edge_it_->v()
                                   , partial_zzfil_ ); //does not modify cpx

//if partial_zzfil_ is empty after flag_lazy_remove (or flag_lazy_insert) ; in case
//the edge or vertex in not in the complex. then *sh_it_ becomes invalid !!
//throw an exception.
        arrow_direction_ = false; //the arrow is backward, these are removals
//sort the simplices to remove by decreasing key values. This ensures that cofaces 
//come before subfaces AND in case of Morse filtration, paired simplices that are 
//removed together MUST have consecutive keys. Consequently, if a simplex that is 
// not critical is removed, and the next simplex is NOT the simplex it's paired 
//with, then we need to make all simplices in the pair critical.
#ifdef GUDHI_USE_TBB
        tbb::parallel_sort( partial_zzfil_.begin(), partial_zzfil_.end()
            , [](Simplex_handle sh1, Simplex_handle sh2)->bool {
                return sh1->second.key() > sh2->second.key();
            });
#else
        sort( partial_zzfil_.begin(), partial_zzfil_.end()
            , [](Simplex_handle sh1, Simplex_handle sh2)->bool {
                return sh1->second.key() > sh2->second.key();
            });
#endif
        //if a paired Morse simplex s is removed, but not the paired simplex t, we 
        //signal it using break_morse_pair()
      }
      //partial_zzfil_ contains at least the new edge, i.e., is non-empty
      fil_ = edge_it_->fil();
      sh_it_ = partial_zzfil_.begin(); 
      ++edge_it_; 
    }
  }


public:
/* Return true iff sh_it is pointing to a backward arrow, with the removal of a   
 * simplex s, paired to a simplex t that is not removed in the current sequence 
 * encoded in partial_zzfil_. This is used to ``break a Morse pair'' in zigzag 
 * persistent homology.
 * Check simplex after and simplex before.
 */
  bool break_morse_pair() {
    if(arrow_direction_) { return false; }//no problem with insertions
    if(cpx_->critical(*sh_it_)) { return false; }//no problem with critical faces
    //are we paired with the next simplex?
    auto sh_it_next = sh_it_;     ++sh_it_next;
    if(sh_it_next != partial_zzfil_.end() //there's something after
       && cpx_->is_paired_with(*sh_it_,*sh_it_next))//sh_it paired with next, ensured by design if both s and t are removed 
    { return false; }//no problem, we remove both simplices

    //are we paired with the previous simplex?
    if(sh_it_ != partial_zzfil_.begin()) {
      auto sh_it_prev = sh_it_;      --sh_it_prev;
      if(cpx_->is_paired_with(*sh_it_,*sh_it_prev)) 
      { return false; }//no problem, we remove both simplices
    }
    return true;//we are breaking a Morse pair
  }
  
/* If sh_it_ points to a non-critical simplex, make both sh_it_ and the simplex 
 * it's paired with critical. Useful if break_morse_pair == true;
 */
  void make_critical() {
    if(!cpx_->critical(*sh_it_)) {
      auto sh_pair = cpx_->paired_with(*sh_it_);
      cpx_->make_critical(*sh_it_);
      cpx_->make_critical(sh_pair);
    }
  }

public:
/* Complex getting modified by the iterator, must be model of 
 * FlagZigzagFilteredComplex.*/
  FlagZigzagFilteredComplex                      * cpx_; 
/* List of insertion and deletion of vertices and edges representing the 
 * zigzag filtration of the 1-skeleton.*/
  std::vector< Edge_type >                         zigzag_edge_filtration_;
/* Maximal dimension d of the flag complex, i.e., the iterator gives the flag 
 * zigzag filtration induced by the insertion and deletion of vertices and edges 
 * in zigzag_edge_filtration_, restricted to the d-skeleton. */
  int                                              dim_max_;
/* A chunk of of the zz filtration, which we have already computed. 
 * ..constructed by the last edge insertion. 
 * When reaching the end of it, clear it, insert a new edge via the edge_it_++ 
 * and compute a new bunch of simplices of the zz filtration. */
  typename std::vector< Simplex_handle >           partial_zzfil_;
  //current simplex in partial_zzfil_, returned by *()
  typename std::vector< Simplex_handle >::iterator sh_it_;
  //iterator in the range of edges; points to the next edge to insert or remove
  typename std::vector< Edge_type >::iterator      edge_it_;
  //true if the simplices in partial_zzfil_ are insertions, and false if deletions
  bool                                             arrow_direction_;
  //counts the total number of insertions in the zigzag, used to assign keys
  int                                              counter_insert;
  //true iff we are finishing emptying the complex
  bool                                             are_we_done;
  //filtration value attached to the arrow
  Filtration_value                                 fil_;
};

#endif //SIMPLEX_TREE_ZIGZAG_ITERATORS_H_



