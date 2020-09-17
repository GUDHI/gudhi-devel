#include <iostream>
#include <fstream>

#ifdef GUDHI_USE_TBB
#include <tbb/tbb.h>
#endif

/** The oscillating Rips zigzag filtration is a zigzag filtration associated to an 
  * ordered set of points \f$P = \{p_0 ... p_{n-1}\}\f$, 
  * and two parameters \f$\mu\f$ 
  * and \f$\nu\f$, such that \f$0 < \nu \leq \mu\f$. 
  *
  * Specifically, denote by \f$P_i\f$ the subset of points \f$\{p_0 ... p_i\}\f$ 
  * for \f$0\leq i \leq n-1\f$, and \f$\mathcal{R}^{\rho}(P_i)\f$ the Rips complex 
  * of threshold \f$\rho \geq 0\f$ built on points \f$P_i\f$. Associate also 
  * to the ordered set of points the values \f$\varepsilon_i\f$, 
  * \f$0\leq i \leq n-1\f$, where \f$\varepsilon_i\f$ is the Hausdorff distance 
  * between the subset of points \f$P_i\f$ and the entire set of points \f$P\f$. 
  * By definition, the sequence \f$\varepsilon_i\f$ is decreasing.
  *
  * Then the oscillating Rips zigzag filtration associated to \f$P\f$ and 
  * parameters \f$\mu\f$ and \f$\nu\f$ as above is the filtration:
  *
  * \f$
  * \mathcal{R}^{\nu \cdot +\infty}(\emptyset)
  *  \subseteq \mathcal{R}^{\mu \cdot +\infty}(P_0)
  *   \supseteq \mathcal{R}^{\nu \cdot \varepsilon_0}(P_0) 
  *    \subseteq \mathcal{R}^{\mu \cdot \varepsilon_0}(P_1) 
  *     \supseteq \mathcal{R}^{\nu \cdot \varepsilon_1}(P_1) 
  *      \subseteq \mathcal{R}^{\mu \cdot \varepsilon_1}(P_2) 
  *       \supseteq \mathcal{R}^{\nu \cdot \varepsilon_2}(P_2) 
  *        \subseteq \cdots 
  *         \supseteq \mathcal{R}^{\nu \cdot \varepsilon_n}(P_n).\f$
  *
  * By convention, we start with the empty complex, as a Rips complex of infinite 
  * threshold on an empty set of vertices.
  *
  * A simplex inserted in the inclusion 
  * \f$\mathcal{R}^{\nu \cdot \varepsilon_i}(P_i) 
  *     \subseteq \mathcal{R}^{\mu \cdot \varepsilon_i}(P_{i+1})\f$ 
  * is given filtration value \f$\varepsilon_i\f$. The vertex $p_0$ inserted as 
  * first simplex is given filtration value \f$+\infty\f$. 
  *
  * A simplex removed in the inclusion 
  * \f$\mathcal{R}^{\mu \cdot \varepsilon_i}(P_{i+1}) 
  *     \supseteq \mathcal{R}^{\nu \cdot \varepsilon_{i+1}}(P_{i+1})\f$ 
  * is given filtration value \f$\varepsilon_{i+1}\f$.  
  */


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

/* Returns vertex with smaller label. */
  typename FilteredComplex::Vertex_handle    u()    { return u_; }
/* Returns vertex with bigger label. */
  typename FilteredComplex::Vertex_handle    v()    { return v_; }
/* Returns the filtration value of the edge. */
  typename FilteredComplex::Filtration_value fil()  { return fil_; }
/* Returns true if insertion of the edge, false if removal. */
  bool                                       type() { return type_; }

  bool operator==(const Zigzag_edge &e) const {
    return ( (e.u_ == u_) && (e.v_ == v_) && 
             (e.fil_ == fil_) && (e.type_ == type_) );
  }

  void assign_fil(typename FilteredComplex::Filtration_value fil) { fil_ = fil; }

private:
  typename FilteredComplex::Vertex_handle    u_;
  typename FilteredComplex::Vertex_handle    v_;
  typename FilteredComplex::Filtration_value fil_;
  bool                                       type_;
};


/** Given an ordered set of points \f$P = \{p_0 ... p_{n-1}\}\f$, 
  * computes the 1-skeleton-filtration corresponding to the oscillating Rips zigzag 
  * filtration of the set of points, i.e., the sequence of insertions and removals 
  * of the oscillating Rips zigzag filtration restricted to vertices and edges. 
  * Because all complexes of the filtration are flag complexes, the 
  * 1-skeleton-filtration encodes fully the entire filtration, in a much more 
  * compact way.
  *
  * The method ordered_points_to_one_skeleton_zigzag_filtration computes the 
  * eps_i = d_H(P_i,P) in 
  * filtration_value[i] = eps_i. A simplex appearing in the inclusion  
  * R({p_0, ... , p_{i}}, nu * eps_i) -> R({p_1, ... , p_i, p_{i+1}}, mu * eps_i)
  * is given filtration value eps_i.
  *
  * filtration_values must be empty, receives the filtration values.
  * edge_filtration must be empty, receives the edge filtration.

  * Edge_t must be of type Zigzag_edge
  */
template< typename Kernel,
		      typename Point_container,
          typename Distance, //furnish()
          typename FiltrationValue,
          typename Edge_t >
void ordered_points_to_one_skeleton_zigzag_filtration(
              							    	Point_container const        & points,
              						        Distance        const          distance,
              						        FiltrationValue const          nu,
              						        FiltrationValue const          mu,
              						        std::vector<FiltrationValue> & filtration_values,
              						        std::vector<Edge_t>          & edge_filtration   )
{
  size_t n = points.size();//number of points

  filtration_values.clear(); filtration_values.resize(n);
  for(size_t i=0; i<n; ++i) 
  { filtration_values[i] = std::numeric_limits<double>::infinity(); }

  //total number of insertion and removal of vertices and edges in the zigzag 
  //filtration
  size_t number_of_arrows = 0;

  //compute all \f$\varepsilon_i\f$ values, such that filtration_values[i] == 
  //eps_i, for i=0 ... n-1:
  for(size_t i=0; i<n; ++i) {
    //entering step i, maintain filtration_value[j] = eps_j for j<i, and
    //filtration_value[k] = d(p_k, P_{i-1}) for k >= i.
#ifdef GUDHI_USE_TBB
    tbb::parallel_for(size_t(i+1), n, [&](size_t k) {
      //set filtration_value[k] <- d(p_k, P_i) == 
      //                           min{ d(p_k, P_{i-1}), d(p_k, p_i) }  for k >= i.
      double dist_i_k = distance(points[i],points[k]); 
      if(dist_i_k < filtration_value[k]) { filtration_value[k] = dist_i_k; }
    } );
#else
    for(size_t k=i+1; k<n; ++k) {
      //set filtration_value[k] <- d(p_k, P_i) == 
      //                           min{ d(p_k, P_{i-1}), d(p_k, p_i) }  for k >= i.
      double dist_i_k = distance(points[i],points[k]); 
      if(dist_i_k < filtration_value[k]) { filtration_value[k] = dist_i_k; }
    }
#endif
//to do: implement parallel version by dividing the vector
    //set filtration_value[i] <- eps_i = d_h(P_i,P) = max_{k>i} d(p_k, P_i) 
    double eps_i = 0.;
    for(size_t k=i+1; k<n; ++k) {
      if(filtration_value[k] > eps_i) { eps_i = filtration_value[k]; }
    }
    filtration_value[i] = eps_i;
  }

  //pre-compute distances (in parallel if possible). Used to accelerate computation
  //of edges according to their lengths.
  //initialize a distance matrix where dist_matrix[i][j] contains the pair 
  //(j, d(p_i,p_j)) for j < i. We sort edges according to length later.
  std::vector< std::vector< std::pair<int, FiltrationValue> > > dist_matrix;
  dist_matrix.resize(n);

#ifdef GUDHI_USE_TBB
  tbb::parallel_for(size_t(0), n, [&](size_t i) {
    // dist_matrix[i] = std::vector< std::pair<int, FiltrationValue> >();
    dist_matrix[i].resize(i);
    for(size_t j=0; j<i; ++j) {
      dist_matrix[i][j] 
                = std::pair<int, FiltrationValue>(j, distance(points[j],points[i]));
    } 
  } );
#else
  for(size_t i=0; i<n; ++i) {//for all vertices
    dist_matrix[i] = std::vector< std::pair<int, FiltrationValue> >();
    dist_matrix[i].resize(i);
    for(size_t j=0; j<i; ++j) {
      dist_matrix[i][j] 
                = std::pair<int, FiltrationValue>(j, distance(points[j],points[i]));
    } 
  }
#endif

/** The two input types std::pair<int, FiltrationValue> encode pairs 
  * (j, d(p_i,p_j)) and (k, d(p_i,p_k)) for some fixed point p_i. 
  * The operator() orders edges by length. By convention, if lengths are equal, 
  * it orders pairs by taking the 
  * smaller vertex label between j and k.
  */
// template<typename FiltrationValue>
  struct point_distance_cmp {
    bool operator()( std::pair<int, FiltrationValue> p1
                   , std::pair<int, FiltrationValue> p2 ) {
      { 
        if(p1.second != p2.second) {return p1.second < p2.second;} //shorter first
        return p1.first < p2.first; 
      }
    }
  };  point_distance_cmp<FiltrationValue> cmp;
  //dist_matrix[i] is now sorted by (j, d(p_i,p_j)) < (k, d(p_i,p_k)) iff 
  //d(p_i,p_j) < d(p_i,p_k) or (j<k in case d(p_i,p_j) == d(p_i,p_k)).
#ifdef GUDHI_USE_TBB
  tbb::parallel_for(size_t(0), n, [&](size_t i) {
    std::stable_sort(dist_matrix[i].begin(), dist_matrix[i].end(), cmp);
  } );
#else
  for(size_t i=0; i<n; ++i) { //all vertices
    std::stable_sort(dist_matrix[i].begin(), dist_matrix[i].end(), cmp);
  }//shortest distance order
#endif

  typename std::vector< std::pair<int, FiltrationValue> >::iterator it;
//edges_added[i] (resp. adges_removed[i]) == list of edges (i,j), with j<i, added (resp. removed) at eps_i
//we also put there (later) vertices that are added
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
  //case i=n-1 has eps_{n-1} = 0 -> no edges
  tbb::parallel_for(size_t(0), n-1, [&](size_t i) {
  //edges_added[i]:
    //consider first all edges added in inclusion:
    //R({p_0, ... , p_i}, nu * eps_i) -> R({p_0, ... , p_i}, mu * eps_i),
    //i.e., all (p_j,p_k) with 0 <= k < j <= i with
    //                                           nu eps_i < d(p_j,p_k) <= mu eps_i  
    for(size_t j = 1; j <= i ; ++j) { 
      //get very first edge (k,j), over all k<j, strictly longer than  mu * eps_i
      //dist_matrix[j][k] = d(p_j,p_k) with k<j
      it = std::upper_bound(dist_matrix[j].begin(), dist_matrix[j].end()
                   , std::pair<int, FiltrationValue>(n, mu * filtration_values[i])
                   , cmp);
      //check
      if(it != dist_matrix[j].end() && it->second <= mu * filtration_values[i]) {
        std::cout << "Edge is not STRICTLY longer.\n";
      }

      while(it != dist_matrix[j].begin()) {
        --it;
        //if edge already in R({p_0, ... , p_i}, nu * eps_i), stop
        if(it->second <= nu * filtration_values[i]) { break; }
        edges_added[i].emplace_back(it->first, j, filtration_values[i], true);
        ++number_of_arrows;
      }
    }
    //now consider all edges added in inclusion:
    //R({p_0, ... , p_i}, mu * eps_i) -> R({p_0, ... , p_i, p_i+1}, mu * eps_i)
    //i.e., all (p_j,p_i+1) with 0 <= j <= i with       d(p_j,p_i+1) <= mu eps_i  
    //first striclty longer edge
    it = std::upper_bound(dist_matrix[i+1].begin(), dist_matrix[i+1].end(), 
            std::pair<int, FiltrationValue>(n, mu * filtration_values[i]), cmp); 
    //check
    if(it != dist_matrix[i+1].end() && it->second <= mu * filtration_values[i]) {
      std::cout << "Edge is not STRICTLY longer.\n";
    }
    while(it != dist_matrix[i+1].begin()) {
      --it;
      edges_added[i].emplace_back(it->first, i+1, filtration_values[i], true);
      ++number_of_arrows;
    }
  //edges_removed[i]:
    //consider all edges removed in
    //R({p_0, ... , p_{i+1}}, mu * eps_i) <- R({p_0, ... , p_{i+1}}, nu * eps_{i+1})
    //i.e., all edges (p_k,p_j), 0<=k<j<=i+1, such that 
    // nu eps_{i+1} < d(p_k,p_j) <= mu eps_i
    for(size_t j = 1; j <= i+1; ++j) {
      //get very first edge (k,j), over all k<j, strictly longer than  mu * eps_i
      //dist_matrix[j][k] = d(p_j,p_k) with k<j
      it = std::upper_bound(dist_matrix[j].begin(), dist_matrix[j].end(), 
             std::pair<int, FiltrationValue>(n, mu * filtration_values[i]), cmp ); 
      //check
      if(it != dist_matrix[j].end() && it->second <= mu * filtration_values[i]) {
        std::cout << "Edge is not STRICTLY longer.\n";
      }      

      while(it != dist_matrix[j].begin()) {
        --it;
        //when reading an edge in R({p_0, ... , p_{i+1}}, nu * eps_{i+1}), stop
        if(it->second <= nu * filtration_values[i+1]) { break; }
        edges_removed[i].emplace_back(it->first, j, filtration_values[i+1], false);
        ++number_of_arrows;
      }
    }
  } );
#else //GUDHI_USE_TBB not defined
  for(size_t i=0; i<n-1; ++i) {
  //edges_added[i]:
    //consider first all edges added in inclusion:
    //R({p_0, ... , p_i}, nu * eps_i) -> R({p_0, ... , p_i}, mu * eps_i),
    //i.e., all (p_j,p_k) with 0 <= k < j <= i with
    //                                           nu eps_i < d(p_j,p_k) <= mu eps_i  
    for(size_t j = 1; j <= i ; ++j) { 
      //get very first edge (k,j), over all k<j, strictly longer than  mu * eps_i
      //dist_matrix[j][k] = d(p_j,p_k) with k<j
      it = std::upper_bound(dist_matrix[j].begin(), dist_matrix[j].end()
                   , std::pair<int, FiltrationValue>(n, mu * filtration_values[i])
                   , cmp);
      //check
      if(it != dist_matrix[j].end() && it->second <= mu * filtration_values[i]) {
        std::cout << "Edge is not STRICTLY longer.\n";
      }

      while(it != dist_matrix[j].begin()) {
        --it;
        //if edge already in R({p_0, ... , p_i}, nu * eps_i), stop
        if(it->second <= nu * filtration_values[i]) { break; }
        edges_added[i].emplace_back(it->first, j, filtration_values[i], true);
        ++number_of_arrows;
      }
    }
    //now consider all edges added in inclusion:
    //R({p_0, ... , p_i}, mu * eps_i) -> R({p_0, ... , p_i, p_i+1}, mu * eps_i)
    //i.e., all (p_j,p_i+1) with 0 <= j <= i with       d(p_j,p_i+1) <= mu eps_i  
    //first striclty longer edge
    it = std::upper_bound(dist_matrix[i+1].begin(), dist_matrix[i+1].end(), 
            std::pair<int, FiltrationValue>(n, mu * filtration_values[i]), cmp); 
    //check
    if(it != dist_matrix[i+1].end() && it->second <= mu * filtration_values[i]) {
      std::cout << "Edge is not STRICTLY longer.\n";
    }
    while(it != dist_matrix[i+1].begin()) {
      --it;
      edges_added[i].emplace_back(it->first, i+1, filtration_values[i], true);
      ++number_of_arrows;
    }
  //edges_removed[i]:
    //consider all edges removed in
    //R({p_0, ... , p_{i+1}}, mu * eps_i) <- R({p_0, ... , p_{i+1}}, nu * eps_{i+1})
    //i.e., all edges (p_k,p_j), 0<=k<j<=i+1, such that 
    // nu eps_{i+1} < d(p_k,p_j) <= mu eps_i
    for(size_t j = 1; j <= i+1; ++j) {
      //get very first edge (k,j), over all k<j, strictly longer than  mu * eps_i
      //dist_matrix[j][k] = d(p_j,p_k) with k<j
      it = std::upper_bound(dist_matrix[j].begin(), dist_matrix[j].end(), 
             std::pair<int, FiltrationValue>(n, mu * filtration_values[i]), cmp ); 
      //check
      if(it != dist_matrix[j].end() && it->second <= mu * filtration_values[i]) {
        std::cout << "Edge is not STRICTLY longer.\n";
      }      

      while(it != dist_matrix[j].begin()) {
        --it;
        //when reading an edge in R({p_0, ... , p_{i+1}}, nu * eps_{i+1}), stop
        if(it->second <= nu * filtration_values[i+1]) { break; }
        edges_removed[i].emplace_back(it->first, j, filtration_values[i+1], false);
        ++number_of_arrows;
      }
    }
  }  
#endif

//Now, sort edges according to their lengths, and put everything in edge_filtration
  edge_filtration.clear();
  edge_filtration.reserve(number_of_arrows + n); //count edges + vertices additions

// Compare edges by distance first (shorter is smaller), and lexicographic ordering
// otherwise. This agrees with the useful Rips filtration in standard persistence.
struct edge_cmp {
  edge_cmp(Point_container &points, Distance distance) 
  : points_(&points), distance_(distance) {}

  bool operator()(Edge_t e1, Edge_t e2) 
  { //lengths of edges e1 and e2
    FiltrationValue dist1 = distance_((*points_)[e1.u()], (*points_)[e1.v()]);
    FiltrationValue dist2 = distance_((*points_)[e2.u()], (*points_)[e2.v()]);
    if(dist1  != dist2)  {return dist1 < dist2;}
    if(e1.u() != e2.u()) {return e1.u() < e2.u();}
    return e1.v() < e2.v(); 
  }
private:
  Point_container  *points_;
  Distance          distance_;
};
//sort insertions and deletion by edge length, then lexicographic order
  edge_cmp<Point_container, Distance, FiltrationValue, Edge_t > e_cmp(points, distance);

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

  //initialise R({p_0}, \nu * eps_0)
  edge_filtration.emplace_back(0, 0, filtration_values[0], true);//add vertex p_0
  for(size_t i = 0; i < n-1; ++i) {//all ascending arrows eps_i
    edge_filtration.emplace_back(i+1, i+1, filtration_values[i], true);//add p_{i+1}
    for(auto edg_it = edges_added[i].begin(); 
            edg_it != edges_added[i].end(); ++edg_it) {
      edge_filtration.push_back(*edg_it);
    }
    for(auto edg_it = edges_removed[i].rbegin(); 
             edg_it != edges_removed[i].rend(); ++edg_it) {
      edge_filtration.push_back(*edg_it);
    }
  }
}

