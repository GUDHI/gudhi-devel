/*
 * Is_manifold.h
 *  Created on: Jan 28, 2015
 * This file is part of the Gudhi Library. The Gudhi library 
 *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
 *    library for computational topology.
 *
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014  INRIA Sophia Antipolis-Mediterranee (France)
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


#ifndef UTILS_IS_MANIFOLD_H_
#define UTILS_IS_MANIFOLD_H_

#include "utils/UI_utils.h"
#include "utils/Edge_contractor.h"

/**
 * Iteratively tries to anticollapse smallest edge non added so far.
 * If its link is contractible then no topological change and else possible topological change.
 *
 * todo do a sparsification with some parameter eps while growing
 */
template<typename SkBlComplex> class Is_manifold {
 private:
  const SkBlComplex& input_complex_;
  typedef typename SkBlComplex::Vertex_handle Vertex_handle;

 public:
  /*
   * return dim the maximum dimension around one simplex and res which is true if the complex is a manifold.
   * If the complex has dimension <= 3 then if res is false, the complex is not a manifold.
   * For d-manifold with d>=4, res may be false while the complex is a manifold.
   */
  Is_manifold(const SkBlComplex& input_complex, unsigned& dim, bool & res) : input_complex_(input_complex) {
    res = true;
    dim = -1;
    if (!input_complex_.empty()) {
      for (auto v : input_complex_.vertex_range()) {
        dim = local_dimension(v);
        break;
      }
      // check that the link of every vertex is a dim-1 sphere
      for (auto v : input_complex_.vertex_range()) {
        if (!is_k_sphere(v, dim - 1)) {
          res = false;
          break;
        }
      }
    }
  }

 private:
  unsigned local_dimension(Vertex_handle v) {
    unsigned dim = 0;
    for (const auto& s : input_complex_.star_simplex_range(v))
      dim = (std::max)(dim, (unsigned) s.dimension());
    return dim;
  }

  bool is_k_sphere(Vertex_handle v, int k) {
    auto link = input_complex_.link(v);
    Edge_contractor<Complex> contractor(link, link.num_vertices() - 1);
    (void)contractor;
    return (is_sphere_simplex(link) == k);
  }

  // A minimal sphere is a complex that contains vertices v1...vn and all faces
  // made upon this set except the face {v1,...,vn}
  // return -2 if not a minimal sphere
  // and d otherwise if complex is a d minimal sphere

  template<typename SubComplex>
  int is_sphere_simplex(const SubComplex& complex) {
    if (complex.empty()) return -1;
    if (complex.num_blockers() != 1) return -2;

    // necessary and sufficient condition : there exists a unique blocker that passes through all vertices
    auto first_blocker = *(complex.const_blocker_range().begin());

    if (first_blocker->dimension() + 1 != complex.num_vertices())
      return -2;
    else
      return (first_blocker->dimension() - 1);
  }
};

#endif  // UTILS_IS_MANIFOLD_H_
