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
 * \brief
 * Given a simplicial complex and an embedding of its vertices, this assigns to
 * each simplex a filtration value equal to the squared radius of its minimal
 * enclosing ball (MEB).
 *
 * Applied on a Čech complex, it recomputes the same values (squared). Applied on a Delaunay triangulation, it computes the Delaunay-Čech filtration.
 *
 * \tparam Kernel CGAL kernel: either Epick_d or Epeck_d.
 * \tparam PointRange Random access range of `Kernel::Point_d`.
 *
 * @param[in] k The geometric kernel.
 * @param[in] complex The simplicial complex.
 * @param[in] points Embedding of the vertices of the complex.
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
      // For a Simplex_tree, this would be a bit faster, but that's probably negligible
      // Vertex_handle u = sh->first; Vertex_handle v = self_siblings(sh)->parent();
      auto verts = complex.simplex_vertex_range(sh);
      auto vert_it = verts.begin();
      Vertex_handle u = *vert_it;
      Vertex_handle v = *++vert_it;
      auto&& pu = points[u];
      Point_d m = k.midpoint_d_object()(pu, points[v]);
      FT r = k.squared_distance_d_object()(m, pu);
      if (exact) CGAL::exact(r);
      complex.assign_key(sh, cache_.size());
      complex.assign_filtration(sh, std::max(cvt(r), Filtration_value(0)));
      cache_.emplace_back(std::move(m), std::move(r));
    } else {
      Filtration_value maxf = 0; // max filtration of the faces
      bool found = false;
      using std::max;
      for (auto face_opposite_vertex : complex.boundary_opposite_vertex_simplex_range(sh)) {
        maxf = max(maxf, complex.filtration(face_opposite_vertex.first));
        if (!found) {
          auto key = complex.key(face_opposite_vertex.first);
          Sphere const& sph = cache_[key];
          if (k.squared_distance_d_object()(sph.first, points[face_opposite_vertex.second]) > sph.second) continue;
          found = true;
          complex.assign_key(sh, key);
          // With exact computations, we could stop here
          // complex.assign_filtration(sh, complex.filtration(face_opposite_vertex.first)); return;
          // but because of possible rounding errors, we continue with the equivalent of make_filtration_non_decreasing
        }
      }
      if (!found) {
        // None of the faces are good enough, MEB must be the circumsphere.
        pts.clear();
        for (auto vertex : complex.simplex_vertex_range(sh))
          pts.push_back(points[vertex]);
        Point_d c = k.construct_circumcenter_d_object()(pts.begin(), pts.end());
        FT r = k.squared_distance_d_object()(c, pts.front());
        if (exact) CGAL::exact(r);
        maxf = max(maxf, cvt(r)); // maxf = cvt(r) except for rounding errors
        complex.assign_key(sh, cache_.size());
        // We could check if the simplex is maximal and avoiding adding it to the cache in that case.
        cache_.emplace_back(std::move(c), std::move(r));
      }
      complex.assign_filtration(sh, maxf);
    }
  };
  complex.for_each_simplex(fun);

  // We could avoid computing maxf, but when !exact rounding errors may cause
  // the filtration values to be non-monotonous, so we would need to call
  // if (!exact) complex.make_filtration_non_decreasing();
  // which is way more costly than computing maxf. The exact case is already so
  // costly that it isn't worth maintaining code without maxf just for it.
  // Cech_complex has "free" access to the max of the faces, because
  // expansion_with_blockers computes it before the callback.

  // TODO: use a map if complex does not provide key?
}
}  // namespace Gudhi::cech_complex

#endif  // MEB_FILTRATION_H_
