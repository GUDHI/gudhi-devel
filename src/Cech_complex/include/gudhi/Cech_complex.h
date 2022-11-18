/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau, Hind Montassif
 *
 *    Copyright (C) 2018 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 *      - 2022/02 Hind Montassif : Replace MiniBall with Sphere_circumradius
 */

#ifndef CECH_COMPLEX_H_
#define CECH_COMPLEX_H_

#include <gudhi/Sphere_circumradius.h>       // for Gudhi::cech_complex::Sphere_circumradius
#include <gudhi/graph_simplicial_complex.h>  // for Gudhi::Proximity_graph
#include <gudhi/Debug_utils.h>               // for GUDHI_CHECK
#include <gudhi/Cech_complex_blocker.h>      // for Gudhi::cech_complex::Cech_blocker

#include <iostream>
#include <stdexcept>  // for exception management

namespace Gudhi {

namespace cech_complex {

/**
 * \class Cech_complex
 * \brief Cech complex class.
 *
 * \ingroup cech_complex
 *
 * \details
 * Cech complex is a simplicial complex where the set of all simplices is filtered
 * by the radius of their minimal enclosing ball and bounded by the given max_radius.
 *
 * \tparam Kernel CGAL kernel: either Epick_d or Epeck_d.
 *
 * \tparam SimplicialComplexForCechComplex furnishes `Vertex_handle` and `Filtration_value` type definition required
 * by `Gudhi::Proximity_graph` and Cech blocker.
 *
 */
template <typename Kernel, typename SimplicialComplexForCechComplex>
class Cech_complex {
 private:
  // Required by compute_proximity_graph
  using Vertex_handle = typename SimplicialComplexForCechComplex::Vertex_handle;
  using Filtration_value = typename SimplicialComplexForCechComplex::Filtration_value;
  using Proximity_graph = Gudhi::Proximity_graph<SimplicialComplexForCechComplex>;

  using cech_blocker = Cech_blocker<SimplicialComplexForCechComplex, Cech_complex, Kernel>;

  using Point_d = typename cech_blocker::Point_d;
  using Point_cloud = std::vector<Point_d>;

  // Numeric type of coordinates in the kernel
  using FT = typename cech_blocker::FT;
  // Sphere is a pair of point and squared radius.
  using Sphere = typename cech_blocker::Sphere;

  public:
  /** \brief Cech_complex constructor from a range of points.
   *
   * @param[in] points Range of points where each point is defined as `kernel::Point_d`.
   * @param[in] max_radius Maximal radius value.
   * @param[in] exact Exact filtration values computation. Not exact if `Kernel` is not <a target="_blank"
   * href="https://doc.cgal.org/latest/Kernel_d/structCGAL_1_1Epeck__d.html">CGAL::Epeck_d</a>.
   * Default is false.
   *
   */
  template<typename InputPointRange >
  Cech_complex(const InputPointRange & points, Filtration_value max_radius, const bool exact = false) : max_radius_(max_radius), exact_(exact) {

    point_cloud_.assign(std::begin(points), std::end(points));

    cech_skeleton_graph_ = Gudhi::compute_proximity_graph<SimplicialComplexForCechComplex>(
        point_cloud_, max_radius_, Sphere_circumradius<Kernel, Filtration_value>(exact));
  }

  /** \brief Initializes the simplicial complex from the proximity graph and expands it until a given maximal
   * dimension, using the Cech blocker oracle.
   *
   * @param[in] complex SimplicialComplexForCech to be created.
   * @param[in] dim_max graph expansion until this given maximal dimension.
   * @exception std::invalid_argument In debug mode, if `complex.num_vertices()` does not return 0.
   *
   */
  void create_complex(SimplicialComplexForCechComplex& complex, int dim_max) {
    GUDHI_CHECK(complex.num_vertices() == 0,
                std::invalid_argument("Cech_complex::create_complex - simplicial complex is not empty"));

    // insert the proximity graph in the simplicial complex
    complex.insert_graph(cech_skeleton_graph_);
    // expand the graph until dimension dim_max
    complex.expansion_with_blockers(dim_max, cech_blocker(&complex, this));
  }

  /** @return max_radius value given at construction. */
  Filtration_value max_radius() const { return max_radius_; }

  /** @param[in] vertex Point position in the range.
   * @return The point.
   */
  const Point_d& get_point(Vertex_handle vertex) const { return point_cloud_[vertex]; }

  /**
   * @return Vector of cached spheres.
   */
  std::vector<Sphere> & get_cache() { return cache_; }

  /** \brief Check exact option
   * @return Exact option.
   */
  const bool is_exact() { return exact_; }

 private:
  Proximity_graph cech_skeleton_graph_;
  Filtration_value max_radius_;
  Point_cloud point_cloud_;
  std::vector<Sphere> cache_;
  const bool exact_;
};

}  // namespace cech_complex

}  // namespace Gudhi

#endif  // CECH_COMPLEX_H_
