/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2018 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CECH_COMPLEX_BLOCKER_H_
#define CECH_COMPLEX_BLOCKER_H_

#include <CGAL/NT_converter.h> // for casting from FT to double

#include <iostream>
#include <vector>
#include <set>
#include <cmath>  // for std::sqrt

namespace Gudhi {

namespace cech_complex {

/** \internal
 * \class Cech_blocker
 * \brief Čech complex blocker.
 *
 * \ingroup cech_complex
 *
 * \details
 * Čech blocker is an oracle constructed from a Cech_complex and a simplicial complex.
 *
 * \tparam SimplicialComplexForProximityGraph furnishes `Simplex_handle` and `Filtration_value` type definition,
 * `simplex_vertex_range(Simplex_handle sh)`and `assign_filtration(Simplex_handle sh, Filtration_value filt)` methods.
 *
 * \tparam Chech_complex is required by the blocker.
 */
template <typename SimplicialComplexForCech, typename Cech_complex, typename Kernel>
class Cech_blocker {
 private:

  using Simplex_handle = typename SimplicialComplexForCech::Simplex_handle;
  using Filtration_value = typename SimplicialComplexForCech::Filtration_value;

 public:

  using Point_d = typename Kernel::Point_d;
  // Numeric type of coordinates in the kernel
  using FT = typename Kernel::FT;
  // Sphere is a pair of point and squared radius.
  using Sphere = typename std::pair<Point_d, FT>;

  template<class PointIterator>
  FT get_squared_radius(PointIterator begin, PointIterator end) const {
    return kernel_.compute_squared_radius_d_object()(begin, end);
  }

  template<class PointIterator>
  Sphere get_sphere(PointIterator begin, PointIterator end) const {
    Point_d c = kernel_.construct_circumcenter_d_object()(begin, end);
    FT r = kernel_.squared_distance_d_object()(c, *begin);
    return std::make_pair(std::move(c), std::move(r));
  }


  /** \internal \brief Čech complex blocker operator() - the oracle - assigns the filtration value from the simplex
   * radius and returns if the simplex expansion must be blocked.
   *  \param[in] sh The Simplex_handle.
   *  \return true if the simplex radius is greater than the Cech_complex max_radius*/
  bool operator()(Simplex_handle sh) {
    using Point_cloud =  std::vector<Point_d>;
    CGAL::NT_converter<FT, double> cast_to_double;
    Filtration_value radius = 0.;

    // for each face of simplex sh, test outsider point is indeed inside enclosing ball, if yes, take it and exit loop, otherwise, new sphere is circumsphere of all vertices
    Sphere min_enclos_ball;
    CGAL::NT_converter<double, FT> cast_to_FT;
    min_enclos_ball.second = cast_to_FT(std::numeric_limits<double>::max());
    Point_cloud face_points;
    for (auto face : sc_ptr_->boundary_simplex_range(sh)) {
        // Find which vertex of sh is missing in face. We rely on the fact that simplex_vertex_range is sorted.
        auto longlist = sc_ptr_->simplex_vertex_range(sh);
        auto shortlist = sc_ptr_->simplex_vertex_range(face);

        auto longiter = std::begin(longlist);
        auto shortiter = std::begin(shortlist);
        auto enditer = std::end(shortlist);
        while(shortiter != enditer && *longiter == *shortiter) { ++longiter; ++shortiter; }
        auto extra = *longiter; // Vertex_handle

        for (auto vertex : sc_ptr_->simplex_vertex_range(face)) {
            face_points.push_back(cc_ptr_->get_point(vertex));
    #ifdef DEBUG_TRACES
        std::clog << "#(" << vertex << ")#";
    #endif  // DEBUG_TRACES
        }
        Sphere sph;
        auto k = sc_ptr_->key(face);
        if(k != sc_ptr_->null_key()) {
            sph = cc_ptr_->get_cache().at(k);
        }
        else {
            sph = get_sphere(face_points.cbegin(), face_points.cend());
        }
        face_points.clear();

        if (kernel_.squared_distance_d_object()(sph.first, cc_ptr_->get_point(extra)) <= sph.second) {
            radius = std::sqrt(cast_to_double(sph.second));
            #ifdef DEBUG_TRACES
                std::clog << "circumcenter: " << sph.first << ", radius: " <<  radius << std::endl;
            #endif  // DEBUG_TRACES
            if (cast_to_double(sph.second) < cast_to_double(min_enclos_ball.second))
                min_enclos_ball = sph;
        }
    }
    // Get the minimal radius of all faces enclosing balls if exists
    if(cast_to_double(min_enclos_ball.second) != std::numeric_limits<double>::max()) {
        radius = std::sqrt(cast_to_double(min_enclos_ball.second));

        sc_ptr_->assign_key(sh, cc_ptr_->get_cache().size());
        cc_ptr_->get_cache().push_back(min_enclos_ball);
    }

    if (radius == 0.) { // Spheres of each face don't contain the whole simplex
        Point_cloud points;
        for (auto vertex : sc_ptr_->simplex_vertex_range(sh)) {
            points.push_back(cc_ptr_->get_point(vertex));
        }
        Sphere sph = get_sphere(points.cbegin(), points.cend());
        radius = std::sqrt(cast_to_double(sph.second));

        sc_ptr_->assign_key(sh, cc_ptr_->get_cache().size());
        cc_ptr_->get_cache().push_back(sph);
    }

#ifdef DEBUG_TRACES
    if (radius > cc_ptr_->max_radius()) std::clog << "radius > max_radius => expansion is blocked\n";
#endif  // DEBUG_TRACES
    sc_ptr_->assign_filtration(sh, radius);
    return (radius > cc_ptr_->max_radius());
  }

  /** \internal \brief Čech complex blocker constructor. */
  Cech_blocker(SimplicialComplexForCech* sc_ptr, Cech_complex* cc_ptr) : sc_ptr_(sc_ptr), cc_ptr_(cc_ptr) {}

 private:
  SimplicialComplexForCech* sc_ptr_;
  Cech_complex* cc_ptr_;
  Kernel kernel_;
};

}  // namespace cech_complex

}  // namespace Gudhi

#endif  // CECH_COMPLEX_BLOCKER_H_
