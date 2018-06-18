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

#ifndef ALPHA_COMPLEX_3D_H_
#define ALPHA_COMPLEX_3D_H_


#include <gudhi/Debug_utils.h>
#include <gudhi/Alpha_complex_3d_options.h>

#include <CGAL/Object.h>
#include <CGAL/tuple.h>
#include <CGAL/iterator.h>


namespace Gudhi {

namespace alpha_complex {

template<typename AlphaComplex3dOptions>
class Alpha_complex_3d {
private:
  using Alpha_shape_3 = typename AlphaComplex3dOptions::Alpha_shape_3;
  // filtration with alpha values needed type definition
  using Alpha_value_type = typename Alpha_shape_3::FT;
  using Object = CGAL::Object;
  using Dispatch =
  CGAL::Dispatch_output_iterator<CGAL::cpp11::tuple<Object, Alpha_value_type>,
  CGAL::cpp11::tuple<std::back_insert_iterator<std::vector<Object> >,
      std::back_insert_iterator<std::vector<Alpha_value_type> > > >;
  using Cell_handle = typename Alpha_shape_3::Cell_handle;
  using Facet = typename Alpha_shape_3::Facet;
  using Edge_3 = typename Alpha_shape_3::Edge;
  using Vertex_handle = typename Alpha_shape_3::Vertex_handle;

public:
  /** \brief Alpha_complex constructor from a list of points.
  *
  * Duplicate points are inserted once in the Alpha_complex. This is the reason why the vertices may be not contiguous.
  *
  * @param[in] points Range of points to triangulate. Points must be in AlphaComplex3dOptions::Point_3
  *
  * The type InputPointRange must be a range for which std::begin and
  * std::end return input iterators on a AlphaComplex3dOptions::Point_3.
  */
  template<typename InputPointRange >
  Alpha_complex_3d(const InputPointRange& points) {
    static_assert(!AlphaComplex3dOptions::weighted, "This constructor is not available for weighted versions of Alpha_complex_3d");
    static_assert(!AlphaComplex3dOptions::periodic, "This constructor is not available for periodic versions of Alpha_complex_3d");
    std::cout << points[0] << std::endl;
    Alpha_shape_3 alpha_shape_3(std::begin(points), std::end(points), 0, Alpha_shape_3::GENERAL);

    Dispatch disp = CGAL::dispatch_output<Object, Alpha_value_type>(std::back_inserter(the_objects),
                                                                    std::back_inserter(the_alpha_values));

    alpha_shape_3.filtration_with_alpha_values(disp);
#ifdef DEBUG_TRACES
    std::cout << "filtration_with_alpha_values returns : " << the_objects.size() << " objects" << std::endl;
#endif  // DEBUG_TRACES

  }

  /** \brief Alpha_complex constructor from a list of points.
*
* Duplicate points are inserted once in the Alpha_complex. This is the reason why the vertices may be not contiguous.
*
* @param[in] points Range of points to triangulate. Points must be in Kernel::Point_d
*
* The type InputPointRange must be a range for which std::begin and
* std::end return input iterators on a Kernel::Point_d.
*/
  template<typename InputPointRange , typename WeightRange>
  Alpha_complex_3d(const InputPointRange& points, WeightRange weights) {
    static_assert(AlphaComplex3dOptions::weighted,
                  "This constructor is not available for non-weighted versions of Alpha_complex_3d");
    static_assert(!AlphaComplex3dOptions::periodic,
                  "This constructor is not available for periodic versions of Alpha_complex_3d");
    GUDHI_CHECK((weights.size() == points.size()),
                std::invalid_argument("Alpha_complex_3d constructor with weights requires points number to be equal with points number"));

    using Weighted_point_3 = typename AlphaComplex3dOptions::Weighted_point_3;
    std::vector<Weighted_point_3> weighted_points_3;

    std::size_t index = 0;
    weighted_points_3.reserve(points.size());
    while ((index < weights.size()) && (index < points.size())) {
      weighted_points_3.push_back(Weighted_point_3(points[index], weights[index]));
      index++;
    }

    Alpha_shape_3 alpha_shape_3(std::begin(weighted_points_3), std::end(weighted_points_3), 0, Alpha_shape_3::GENERAL);

    Dispatch disp = CGAL::dispatch_output<Object, Alpha_value_type>(std::back_inserter(the_objects),
                                                                    std::back_inserter(the_alpha_values));

    alpha_shape_3.filtration_with_alpha_values(disp);
#ifdef DEBUG_TRACES
    std::cout << "filtration_with_alpha_values returns : " << the_objects.size() << " objects" << std::endl;
#endif  // DEBUG_TRACES
  }

private:
  std::vector<Object> the_objects;
  std::vector<Alpha_value_type> the_alpha_values;

};

}  // namespace alpha_complex

}  // namespace Gudhi

#endif  // ALPHA_COMPLEX_3D_H_
