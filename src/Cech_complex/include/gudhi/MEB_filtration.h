/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Marc Glisse
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef MEB_FILTRATION_H_
#define MEB_FILTRATION_H_

namespace Gudhi::cech_complex {

/**
 * \ingroup cech_complex
 *
 * Given a simplicial complex and an embedding of its vertices, this assigns to
 * each simplex a filtration value equal to the squared radius of its minimal
 * enclosing ball.
 *
 * Applied on a Cech complex, it recomputes the same values (squared). Applied on a Delaunay triangulation, it computes the Delaunay-Cech filtration.
 *
 * \tparam Kernel CGAL kernel: either Epick_d or Epeck_d.
 * \tparam SimplicialComplexForMEB Simplex_tree with Simplex_key an integer type large enough to store the number of simplices.
 * \tparam PointRange Random access range of `Kernel::Point_d`.
 *
 * @param[in] points Embedding of the vertices of the complex.
 * @param[in] complex The simplicial complex.
 * @param[in] exact If true and `Kernel` is <a href="https://doc.cgal.org/latest/Kernel_d/structCGAL_1_1Epeck__d.html">CGAL::Epeck_d</a>, the filtration values are computed exactly. Default is false.
 */

template<typename Kernel, typename SimplicialComplexForMEB, typename PointRange>
void assign_MEB_filtration(Kernel&&k, SimplicialComplexForMEB& complex, PointRange points, bool exact = false) {
  using Point_d = typename Kernel::Point_d;
  using FT = typename Kernel::FT;
  using Sphere = std::pair<Point_d, FT>;

  using Vertex_handle = typename SimplicialComplexForMEB::Vertex_handle;
  using Simplex_handle = typename SimplicialComplexForMEB::Simplex_handle;
  using Filtration_value = typename SimplicialComplexForMEB::Filtration_value;

  std::vector<Sphere> cache_;
  std::vector<Point_d> pts;
  CGAL::NT_converter<FT, Filtration_value> cvt;

  auto fun = [&](Simplex_handle sh, int dim){
    if (dim == 0) complex.assign_filtration(sh, 0);
    else if (dim == 1) {
      // Vertex_handle u = sh->first;
      // Vertex_handle v = self_siblings(sh)->parent();
      auto verts = complex.simplex_vertex_range(sh);
      auto vert_it = verts.begin();
      Vertex_handle u = *vert_it;
      Vertex_handle v = *++vert_it;
      auto&& pu = points[u];
      Point_d m = k.midpoint_d_object()(pu, points[v]);
      FT r = k.squared_distance_d_object()(m, pu);
      if (exact) CGAL::exact(r);
      complex.assign_key(sh, cache_.size()); // TODO: use a map if complex does not provide key?
      complex.assign_filtration(sh, std::max(cvt(r), Filtration_value(0)));
      cache_.emplace_back(std::move(m), std::move(r));
    } else {
      for (auto face_opposite_vertex : complex.boundary_opposite_vertex_simplex_range(sh)) {
        auto key = complex.key(face_opposite_vertex.first);
        Sphere const& sph = cache_[key];
        if (k.squared_distance_d_object()(sph.first, points[face_opposite_vertex.second]) > sph.second) continue;
        complex.assign_key(sh, key);
        complex.assign_filtration(sh, complex.filtration(face_opposite_vertex.first));
        return;
      }
      // None of the faces are good enough, MEB must be the circumsphere.
      pts.clear();
      for (auto vertex : complex.simplex_vertex_range(sh))
        pts.push_back(points[vertex]);
      Point_d c = k.construct_circumcenter_d_object()(pts.begin(), pts.end());
      FT r = k.squared_distance_d_object()(c, pts.front());
      if (exact) CGAL::exact(r);
      complex.assign_key(sh, cache_.size()); // TODO: use a map if complex does not provide key?
      complex.assign_filtration(sh, cvt(r));
      cache_.emplace_back(std::move(c), std::move(r));
    }
  };
  complex.for_each_simplex_with_dim(fun);
  // TODO: when !exact, how do we ensure that the value is indeed larger
  // than for the faces? Cech_complex also looks at the max of the faces,
  // given by expansion. We would need to do compute this max as well, and
  // in particular not short-circuit when we find a good face. For now, the
  // simplest:
  if (!exact) complex.make_filtration_non_decreasing();
}
}  // namespace Gudhi::cech_complex

#endif  // MEB_FILTRATION_H_
