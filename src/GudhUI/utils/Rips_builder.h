/* This file is part of the Gudhi Library. The Gudhi library 
 *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
 *    library for computational topology.
 *
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
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
 * 
 */

#ifndef UTILS_RIPS_BUILDER_H_
#define UTILS_RIPS_BUILDER_H_

#include <boost/iterator/iterator_facade.hpp>

#include <CGAL/Euclidean_distance.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_d.h>

#include "utils/UI_utils.h"
#include "model/Complex_typedefs.h"

template<typename SkBlComplex> class Rips_builder {
 private:
  SkBlComplex& complex_;

 public:
  /**
   * @brief Modify complex to be the Rips complex
   * of its points with offset alpha.
   */
  Rips_builder(SkBlComplex& complex, double alpha) : complex_(complex) {
    complex.keep_only_vertices();
    if (alpha <= 0) return;
    compute_edges(alpha);
  }

 private:
  double squared_eucl_distance(const Point& p1, const Point& p2) const {
    return Geometry_trait::Squared_distance_d()(p1, p2);
  }

  void compute_edges(double alpha) {
    auto vertices = complex_.vertex_range();
    for (auto p = vertices.begin(); p != vertices.end(); ++p) {
      std::cout << *p << " ";
      std::cout.flush();
      for (auto q = p; ++q != vertices.end(); /**/)
        if (squared_eucl_distance(complex_.point(*p), complex_.point(*q)) < 4 * alpha * alpha)
          complex_.add_edge_without_blockers(*p, *q);
    }
    std::cout << std::endl;
  }
};

#endif  // UTILS_RIPS_BUILDER_H_
