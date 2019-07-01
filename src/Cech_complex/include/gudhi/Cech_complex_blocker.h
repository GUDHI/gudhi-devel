/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2018 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CECH_COMPLEX_BLOCKER_H_
#define CECH_COMPLEX_BLOCKER_H_

#include <gudhi/distance_functions.h>  // for Gudhi::Minimal_enclosing_ball_radius

#include <iostream>
#include <vector>
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
template <typename SimplicialComplexForCech, typename Cech_complex>
class Cech_blocker {
 private:
  using Point_cloud = typename Cech_complex::Point_cloud;

  using Simplex_handle = typename SimplicialComplexForCech::Simplex_handle;
  using Filtration_value = typename SimplicialComplexForCech::Filtration_value;

 public:
  /** \internal \brief Čech complex blocker operator() - the oracle - assigns the filtration value from the simplex
   * radius and returns if the simplex expansion must be blocked.
   *  \param[in] sh The Simplex_handle.
   *  \return true if the simplex radius is greater than the Cech_complex max_radius*/
  bool operator()(Simplex_handle sh) {
    Point_cloud points;
    for (auto vertex : sc_ptr_->simplex_vertex_range(sh)) {
      points.push_back(cc_ptr_->get_point(vertex));
#ifdef DEBUG_TRACES
      std::cout << "#(" << vertex << ")#";
#endif  // DEBUG_TRACES
    }
    Filtration_value radius = Gudhi::Minimal_enclosing_ball_radius()(points);
#ifdef DEBUG_TRACES
    if (radius > cc_ptr_->max_radius()) std::cout << "radius > max_radius => expansion is blocked\n";
#endif  // DEBUG_TRACES
    sc_ptr_->assign_filtration(sh, radius);
    return (radius > cc_ptr_->max_radius());
  }

  /** \internal \brief Čech complex blocker constructor. */
  Cech_blocker(SimplicialComplexForCech* sc_ptr, Cech_complex* cc_ptr) : sc_ptr_(sc_ptr), cc_ptr_(cc_ptr) {}

 private:
  SimplicialComplexForCech* sc_ptr_;
  Cech_complex* cc_ptr_;
};

}  // namespace cech_complex

}  // namespace Gudhi

#endif  // CECH_COMPLEX_BLOCKER_H_
