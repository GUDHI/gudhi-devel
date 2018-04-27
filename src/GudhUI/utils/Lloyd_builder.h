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

#ifndef UTILS_LLOYD_BUILDER_H_
#define UTILS_LLOYD_BUILDER_H_

#include <vector>
#include <list>

/**
 * Iteratively puts every vertex at the center of its neighbors
 */
template<typename SkBlComplex> class Lloyd_builder {
 private:
  SkBlComplex& complex_;
  int dim;

 public:
  typedef typename SkBlComplex::Vertex_handle Vertex_handle;

  /**
   * @brief Modify complex to be the expansion of the k-nearest neighbor
   * symetric graph.
   */
  Lloyd_builder(SkBlComplex& complex, unsigned num_iterations) : complex_(complex), dim(-1) {
    if (!complex_.empty()) {
      dim = get_dimension();
      while (num_iterations--) {
        std::list<Point> new_points;
        for (auto v : complex.vertex_range())
          new_points.push_back(barycenter_neighbors(v));

        auto new_points_it = new_points.begin();
        for (auto v : complex.vertex_range())
          complex_.point(v) = *(new_points_it++);
      }
    }
  }

 private:
  int get_dimension() {
    assert(!complex_.empty());
    for (auto v : complex_.vertex_range())
      return complex_.point(v).dimension();
    return -1;
  }

  Point barycenter_neighbors(Vertex_handle v) const {
    if (complex_.degree(v) == 0)
      return complex_.point(v);

    std::vector<double> res(dim, 0);
    unsigned num_points = 0;
    for (auto nv : complex_.vertex_range(v)) {
      ++num_points;
      const Point& point = complex_.point(nv);
      assert(point.dimension() == dim);
      for (int i = 0; i < point.dimension(); ++i)
        res[i] += point[i];
    }
    for (auto& x : res)
      x /= num_points;
    return Point(dim, res.begin(), res.end());
  }

  double squared_eucl_distance(const Point& p1, const Point& p2) const {
    return Geometry_trait::Squared_distance_d()(p1, p2);
  }
};

#endif  // UTILS_LLOYD_BUILDER_H_
