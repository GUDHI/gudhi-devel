/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2018 Inria
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CECH_COMPLEX_BLOCKER_H_
#define CECH_COMPLEX_BLOCKER_H_

#include <gudhi/Cech_complex.h>        // Cech_blocker is using a pointer on Gudhi::cech_complex::Cech_complex
#include <gudhi/distance_functions.h>  // for Gudhi::Minimal_enclosing_ball_radius

#include <iostream>
#include <vector>
#include <cmath>  // for std::sqrt

namespace Gudhi {

namespace cech_complex {

// Just declaring Cech_complex class because used and not yet defined.
template<typename SimplicialComplexForCechComplex, typename InputPointRange>
class Cech_complex;

/** \internal
 * \class Cech_blocker
 * \brief Cech complex blocker.
 *
 * \ingroup cech_complex
 *
 * \details
 * Cech blocker is an oracle constructed from a Cech_complex and a simplicial complex.
 *
 * \tparam SimplicialComplexForProximityGraph furnishes `Simplex_handle` and `Filtration_value` type definition,
 * `simplex_vertex_range(Simplex_handle sh)`and `assign_filtration(Simplex_handle sh, Filtration_value filt)` methods.
 *
 * \tparam InputPointRange is required by the pointer on Chech_complex for type definition.
 */
template <typename SimplicialComplexForCech, typename InputPointRange>
class Cech_blocker {
 private:
  using Point = std::vector<double>;
  using Point_cloud = std::vector<Point>;
  using Simplex_handle = typename SimplicialComplexForCech::Simplex_handle;
  using Filtration_value = typename SimplicialComplexForCech::Filtration_value;
  using Cech_complex = Gudhi::cech_complex::Cech_complex<SimplicialComplexForCech, InputPointRange>;

 public:
  /** \internal \brief Cech complex blocker operator() - the oracle - assigns the filtration value from the simplex
   * radius and returns if the simplex expansion must be blocked.
   *  \param[in] sh The Simplex_handle.
   *  \return true if the simplex radius is greater than the Cech_complex max_radius*/
  bool operator()(Simplex_handle sh) {
    Point_cloud points;
    for (auto vertex : sc_ptr_->simplex_vertex_range(sh)) {
      points.push_back(Point(cc_ptr_->point_iterator(vertex)->begin(),
                             cc_ptr_->point_iterator(vertex)->end()));
#ifdef DEBUG_TRACES
      std::cout << "#(" << vertex << ")#";
#endif  // DEBUG_TRACES
    }
    Filtration_value radius = Gudhi::Minimal_enclosing_ball_radius()(points);
#ifdef DEBUG_TRACES
    if (radius > cc_ptr_->max_radius())
      std::cout << "radius > max_radius => expansion is blocked\n";
#endif  // DEBUG_TRACES
    sc_ptr_->assign_filtration(sh, radius);
    return (radius > cc_ptr_->max_radius());
  }

  /** \internal \brief Cech complex blocker constructor. */
  Cech_blocker(SimplicialComplexForCech* sc_ptr, Cech_complex* cc_ptr)
    : sc_ptr_(sc_ptr),
      cc_ptr_(cc_ptr) {
  }
 private:
  SimplicialComplexForCech* sc_ptr_;
  Cech_complex* cc_ptr_;
};

}  // namespace cech_complex

}  // namespace Gudhi

#endif  // CECH_COMPLEX_BLOCKER_H_
