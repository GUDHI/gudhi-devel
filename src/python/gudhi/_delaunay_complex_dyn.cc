/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - 2024/03 Vincent Rouvreau: Renamed Alpha_complex_factory as Delaunay_complex_factory for DelaunayCechComplex.
 *                                  Factorize create_complex
 *      - 2024/10 Vincent Rouvreau: Add square root filtration values interface
 *      - 2025/03 Thibaud Kloczko: Use nanobind instead of Cython for python bindings
 *      - 2025/04 Hannah Schreiber: Re-add possibility of tensors (numpy, torch etc.) as input
 *      - YYYY/MM Author: Description of the modification
 */

#include <memory>  // for std::unique_ptr, std::make_unique

#include <CGAL/Epeck_d.h>
#include <CGAL/Epick_d.h>

#include <python_interfaces/delaunay_complex_interface.h>
#include <python_interfaces/points_utils.h>

namespace Gudhi {
namespace delaunay_complex {

std::unique_ptr<Abstract_delaunay_complex> make_dynamic_d_epick_complex_ptr(const Sequence2D& points,
                                                                            const Sequence1D& weights,
                                                                            bool exact_version)
{
  return std::make_unique<Delaunay_complex_t<CGAL::Epick_d<CGAL::Dynamic_dimension_tag>, true>>(
      points, weights, exact_version);
}

std::unique_ptr<Abstract_delaunay_complex> make_dynamic_d_epeck_complex_ptr(const Sequence2D& points,
                                                                            const Sequence1D& weights,
                                                                            bool exact_version)
{
  return std::make_unique<Delaunay_complex_t<CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>, true>>(
      points, weights, exact_version);
}

std::unique_ptr<Abstract_delaunay_complex> make_dynamic_d_epick_complex_ptr(const Sequence2D& points,
                                                                            bool exact_version)
{
  return std::make_unique<Delaunay_complex_t<CGAL::Epick_d<CGAL::Dynamic_dimension_tag>, false>>(points, exact_version);
}

std::unique_ptr<Abstract_delaunay_complex> make_dynamic_d_epeck_complex_ptr(const Sequence2D& points,
                                                                            bool exact_version)
{
  return std::make_unique<Delaunay_complex_t<CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>, false>>(points, exact_version);
}

}  // namespace delaunay_complex
}  // namespace Gudhi
