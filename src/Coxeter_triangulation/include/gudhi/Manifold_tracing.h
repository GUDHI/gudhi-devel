/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2019 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef MANIFOLD_TRACING_H_
#define MANIFOLD_TRACING_H_

#include <gudhi/IO/output_debug_traces_to_html.h>  // for DEBUG_TRACES
#include <gudhi/Query_result.h>

#include <boost/functional/hash.hpp>

#include <Eigen/Dense>

#include <queue>
#include <unordered_map>

namespace Gudhi {

namespace coxeter_triangulation {

/**
 *  \ingroup coxeter_triangulation
 */

/** \class Manifold_tracing
 *  \brief A class that assembles methods for manifold tracing algorithm.
 *
 *  \tparam Triangulation_ The type of the ambient triangulation.
 *   Needs to be a model of the concept TriangulationForManifoldTracing.
 */
template <class Triangulation_>
class Manifold_tracing {
  typedef typename Triangulation_::Simplex_handle Simplex_handle;

  struct Simplex_hash {
    typedef Simplex_handle argument_type;
    typedef std::size_t result_type;
    result_type operator()(const argument_type& s) const noexcept {
      return boost::hash<typename Simplex_handle::Vertex>()(s.vertex());
    }
  };

 public:
  /** \brief Type of the output simplex map with keys of type Triangulation_::Simplex_handle
   *   and values of type Eigen::VectorXd.
   *   This type should be used for the output in the method manifold_tracing_algorithm.
   */
  typedef std::unordered_map<Simplex_handle, Eigen::VectorXd, Simplex_hash> Out_simplex_map;

  /**
   * \brief Computes the set of k-simplices that intersect
   * a boundaryless implicit manifold given by an intersection oracle, where k
   * is the codimension of the manifold.
   * The computation is based on the seed propagation --- it starts at the
   * given seed points and then propagates along the manifold.
   *
   * \tparam Point_range Range of points of type Eigen::VectorXd.
   * \tparam Intersection_oracle Intersection oracle that represents the manifold.
   *  Needs to be a model of the concept IntersectionOracle.
   *
   * \param[in] seed_points The range of points on the manifold from which
   * the computation begins.
   * \param[in] triangulation The ambient triangulation.
   * \param[in] oracle The intersection oracle for the manifold.
   * The ambient dimension needs to match the dimension of the
   * triangulation.
   * \param[out] out_simplex_map The output map, where the keys are k-simplices in
   * the input triangulation that intersect the input manifold and the mapped values
   * are the intersection points.
   */
  template <class Point_range, class Intersection_oracle>
  void manifold_tracing_algorithm(const Point_range& seed_points, const Triangulation_& triangulation,
                                  const Intersection_oracle& oracle, Out_simplex_map& out_simplex_map) {
    std::size_t cod_d = oracle.cod_d();
    std::queue<Simplex_handle> queue;

    for (const auto& p : seed_points) {
      Simplex_handle full_simplex = triangulation.locate_point(p);
      for (Simplex_handle face : full_simplex.face_range(cod_d)) {
        Query_result<Simplex_handle> qr = oracle.intersects(face, triangulation);
        if (qr.success && out_simplex_map.emplace(std::make_pair(face, qr.intersection)).second) {
#ifdef DEBUG_TRACES
          mt_seed_inserted_list.push_back(MT_inserted_info(qr, face, false));
#endif
          queue.emplace(face);
          break;
        }
      }
    }

    while (!queue.empty()) {
      Simplex_handle s = queue.front();
      queue.pop();
      for (auto cof : s.coface_range(cod_d + 1)) {
        for (auto face : cof.face_range(cod_d)) {
          Query_result<Simplex_handle> qr = oracle.intersects(face, triangulation);
          if (qr.success && out_simplex_map.emplace(std::make_pair(face, qr.intersection)).second) queue.emplace(face);
        }
      }
    }
  }

  /**
   * \brief Computes the set of k-simplices that intersect
   * the dimensional manifold given by an intersection oracle, where k
   * is the codimension of the manifold.
   * The computation is based on the seed propagation --- it starts at the
   * given seed points and then propagates along the manifold.
   *
   * \tparam Point_range Range of points of type Eigen::VectorXd.
   * \tparam Intersection_oracle Intersection oracle that represents the manifold.
   *  Needs to be a model of the concept IntersectionOracle.
   *
   * \param[in] seed_points The range of points on the manifold from which
   * the computation begins.
   * \param[in] triangulation The ambient triangulation.
   * \param[in] oracle The intersection oracle for the manifold.
   * The ambient dimension needs to match the dimension of the
   * triangulation.
   * \param[out] interior_simplex_map The output map, where the keys are k-simplices in
   * the input triangulation that intersect the relative interior of the input manifold
   * and the mapped values are the intersection points.
   * \param[out] boundary_simplex_map The output map, where the keys are k-simplices in
   * the input triangulation that intersect the boundary of the input manifold
   * and the mapped values are the intersection points.
   */
  template <class Point_range, class Intersection_oracle>
  void manifold_tracing_algorithm(const Point_range& seed_points, const Triangulation_& triangulation,
                                  const Intersection_oracle& oracle, Out_simplex_map& interior_simplex_map,
                                  Out_simplex_map& boundary_simplex_map) {
    std::size_t cod_d = oracle.cod_d();
    std::queue<Simplex_handle> queue;

    for (const auto& p : seed_points) {
      Simplex_handle full_simplex = triangulation.locate_point(p);
      for (Simplex_handle face : full_simplex.face_range(cod_d)) {
        auto qr = oracle.intersects(face, triangulation);
#ifdef DEBUG_TRACES
        mt_seed_inserted_list.push_back(MT_inserted_info(qr, face, false));
#endif
        if (qr.success) {
          if (oracle.lies_in_domain(qr.intersection, triangulation)) {
            if (interior_simplex_map.emplace(std::make_pair(face, qr.intersection)).second) queue.emplace(face);
          } else {
            for (Simplex_handle cof : face.coface_range(cod_d + 1)) {
              auto qrb = oracle.intersects_boundary(cof, triangulation);
#ifdef DEBUG_TRACES
              mt_seed_inserted_list.push_back(MT_inserted_info(qrb, cof, true));
#endif
              if (qrb.success) boundary_simplex_map.emplace(cof, qrb.intersection);
            }
          }
          // break;
        }
      }
    }

    while (!queue.empty()) {
      Simplex_handle s = queue.front();
      queue.pop();
      for (auto cof : s.coface_range(cod_d + 1)) {
        for (auto face : cof.face_range(cod_d)) {
          auto qr = oracle.intersects(face, triangulation);
#ifdef DEBUG_TRACES
          mt_inserted_list.push_back(MT_inserted_info(qr, face, false));
#endif
          if (qr.success) {
            if (oracle.lies_in_domain(qr.intersection, triangulation)) {
              if (interior_simplex_map.emplace(face, qr.intersection).second) queue.emplace(face);
            } else {
              auto qrb = oracle.intersects_boundary(cof, triangulation);
#ifdef DEBUG_TRACES
              mt_inserted_list.push_back(MT_inserted_info(qrb, cof, true));
#endif
              // assert (qrb.success); // always a success
              if (qrb.success) boundary_simplex_map.emplace(cof, qrb.intersection);
            }
          }
        }
      }
    }
  }

  /** \brief Empty constructor */
  Manifold_tracing() {}
};

/**
 * \brief Static method for Manifold_tracing<Triangulation_>::manifold_tracing_algorithm
 * that computes the set of k-simplices that intersect
 * a boundaryless implicit manifold given by an intersection oracle, where k
 * is the codimension of the manifold.
 * The computation is based on the seed propagation --- it starts at the
 * given seed points and then propagates along the manifold.
 *
 * \tparam Point_range Range of points of type Eigen::VectorXd.
 *  \tparam Triangulation_ The type of the ambient triangulation.
 *   Needs to be a model of the concept TriangulationForManifoldTracing.
 * \tparam Intersection_oracle Intersection oracle that represents the manifold.
 *  Needs to be a model of the concept IntersectionOracle.
 * \tparam Out_simplex_map Needs to be Manifold_tracing<Triangulation_>::Out_simplex_map.
 *
 * \param[in] seed_points The range of points on the manifold from which
 * the computation begins.
 * \param[in] triangulation The ambient triangulation.
 * \param[in] oracle The intersection oracle for the manifold.
 * The ambient dimension needs to match the dimension of the
 * triangulation.
 * \param[out] out_simplex_map The output map, where the keys are k-simplices in
 * the input triangulation that intersect the input manifold and the mapped values
 * are the intersection points.
 *
 * \ingroup coxeter_triangulation
 */
template <class Point_range, class Triangulation, class Intersection_oracle, class Out_simplex_map>
void manifold_tracing_algorithm(const Point_range& seed_points, const Triangulation& triangulation,
                                const Intersection_oracle& oracle, Out_simplex_map& out_simplex_map) {
  Manifold_tracing<Triangulation> mt;
  mt.manifold_tracing_algorithm(seed_points, triangulation, oracle, out_simplex_map);
}

/**
 * \brief Static method for Manifold_tracing<Triangulation_>::manifold_tracing_algorithm
 * the dimensional manifold given by an intersection oracle, where k
 * is the codimension of the manifold.
 * The computation is based on the seed propagation --- it starts at the
 * given seed points and then propagates along the manifold.
 *
 * \tparam Point_range Range of points of type Eigen::VectorXd.
 *  \tparam Triangulation_ The type of the ambient triangulation.
 *   Needs to be a model of the concept TriangulationForManifoldTracing.
 * \tparam Intersection_oracle Intersection oracle that represents the manifold.
 *  Needs to be a model of the concept IntersectionOracle.
 * \tparam Out_simplex_map Needs to be Manifold_tracing<Triangulation_>::Out_simplex_map.
 *
 * \param[in] seed_points The range of points on the manifold from which
 * the computation begins.
 * \param[in] triangulation The ambient triangulation.
 * \param[in] oracle The intersection oracle for the manifold.
 * The ambient dimension needs to match the dimension of the
 * triangulation.
 * \param[out] interior_simplex_map The output map, where the keys are k-simplices in
 * the input triangulation that intersect the relative interior of the input manifold
 * and the mapped values are the intersection points.
 * \param[out] boundary_simplex_map The output map, where the keys are k-simplices in
 * the input triangulation that intersect the boundary of the input manifold
 * and the mapped values are the intersection points.
 *
 * \ingroup coxeter_triangulation
 */
template <class Point_range, class Triangulation, class Intersection_oracle, class Out_simplex_map>
void manifold_tracing_algorithm(const Point_range& seed_points, const Triangulation& triangulation,
                                const Intersection_oracle& oracle, Out_simplex_map& interior_simplex_map,
                                Out_simplex_map& boundary_simplex_map) {
  Manifold_tracing<Triangulation> mt;
  mt.manifold_tracing_algorithm(seed_points, triangulation, oracle, interior_simplex_map, boundary_simplex_map);
}

}  // namespace coxeter_triangulation

}  // namespace Gudhi

#endif
