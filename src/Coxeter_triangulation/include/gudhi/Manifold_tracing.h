#ifndef MANIFOLD_TRACING_H_
#define MANIFOLD_TRACING_H_

#include <queue>
#include <unordered_set>
#include <unordered_map>

namespace Gudhi {

class Manifold_tracing {

public:  
  Manifold_tracing() {}


  /* \internal \brief Computes the set of (d-m)-simplices that intersect
   * the m-dimensional manifold given by an intersection oracle.
   * The computation is based on the seed propagation - it starts at the 
   * given seed points and then propagates along the manifold.
   * \param[in] seed_points The range of points on the manifold from which 
   * the computation begins.
   * \param[in] triangulation The ambient triangulation.
   * \param[in] level The scale parameter for the ambient triangulation.
   * \param[in] oracle The intersection oracle for the manifold.
   * The ambient dimension needs to match the dimension of the
   * triangulation trian.
   * \param[out] out_simplices The output set of the (d-m)-simplices in
   * the input triangulation that intersect the input manifold.
   */
  template <class Point_range,
            class Triangulation,
	    class Intersection_oracle >
  void compute_complex(const Point_range& seed_points,
		       const Triangulation& triangulation,
		       double level,
		       const Intersection_oracle& oracle,
		       std::unordered_set<typename Triangulation::Simplex_handle>& out_simplices)
  {
    /* initialization */
    using Simplex_handle = typename Triangulation::Simplex_handle;
    std::size_t amb_d = oracle.amb_d();
    std::size_t cod_d = oracle.cod_d();
    std::unordered_set<Simplex_handle> facet_cells;
    std::unordered_map<Simplex_handle, std::size_t> order_map;
    std::queue<Simplex_handle> queue;

    for (const auto& p: seed_points) {
      Simplex_handle c = triangulation.locate_point(p, level);
      for (auto f: triangulation.face_range(c, cod_d))
	if (oracle.intersects(f, triangulation) &&
	    out_simplices.emplace(f).second &&
	    order_map.emplace(std::make_pair(f, order_map.size())).second)
	  queue.emplace(f);
    }
    
    /* propagation */
    while (!queue.empty()) {
      Simplex_handle s = queue.front();
      queue.pop();
      for (auto cof: triangulation.coface_range(s, cod_d+1))
	if (facet_cells.emplace(cof).second)
	  for (auto f: triangulation.face_range(cof, cod_d))
	    if (oracle.intersects(f, triangulation) &&
		out_simplices.emplace(f).second &&
		order_map.emplace(std::make_pair(f, order_map.size())).second)
	      queue.emplace(f);
    }

  }
};

}

#endif
