/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
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
      std::clog << *p << " ";
      std::clog.flush();
      for (auto q = p; ++q != vertices.end(); /**/)
        if (squared_eucl_distance(complex_.point(*p), complex_.point(*q)) < 4 * alpha * alpha)
          complex_.add_edge_without_blockers(*p, *q);
    }
    std::clog << std::endl;
  }
};

#endif  // UTILS_RIPS_BUILDER_H_
