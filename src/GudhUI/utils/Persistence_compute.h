/* This file is part of the Gudhi Library. The Gudhi library 
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


#ifndef UTILS_PERSISTENCE_COMPUTE_H_
#define UTILS_PERSISTENCE_COMPUTE_H_

#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Rips_complex.h>

#include <vector>

struct Persistence_params {
  int p;
  double threshold;
  int max_dim;
  double min_pers;

  Persistence_params(int p_, double th_, int max_dim_ = 10, double min_pers_ = 0)
      : p(p_), threshold(th_), max_dim(max_dim_), min_pers(min_pers_) { }
};

/**
 * Show persistence into output stream
 */
template<typename SkBlComplex> class Persistence_compute {
 public:
  typedef typename SkBlComplex::Vertex_handle Vertex_handle;
  typedef typename SkBlComplex::Edge_handle Edge_handle;

  /**
   * @brief Compute persistence
   * parameters :
   * unsigned dim_max
   * double threshold
   * int p for coefficient Z_p
   */
  Persistence_compute(SkBlComplex& complex, std::ostream& stream, const Persistence_params& params) {
    // for now everything is copied, todo boost adapt iterators to points of SkBlComplex instead of copying to an
    // initial vector
    typedef std::vector<double> Point_t;
    std::vector< Point_t > points;
    points.reserve(complex.num_vertices());
    for (auto v : complex.vertex_range()) {
      const auto & pt = complex.point(v);
      Point_t pt_to_add(pt.cartesian_begin(), pt.cartesian_end());
      points.emplace_back(std::move(pt_to_add));
    }

    using Simplex_tree = Gudhi::Simplex_tree<>;
    using Filtration_value = Simplex_tree::Filtration_value;
    using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;
    using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
    using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp>;

    Rips_complex rips_complex(points, params.threshold, Euclidean_distance());

    Simplex_tree st;
    if (rips_complex.create_complex(st, params.max_dim)) {
      Persistent_cohomology pcoh(st);
      // initializes the coefficient field for homology
      pcoh.init_coefficients(params.p);
      // put params.min_pers
      pcoh.compute_persistent_cohomology(params.min_pers);
      stream << "persistence: \n";
      stream << "p dimension birth death: \n";
      pcoh.output_diagram(stream);
    }
  }
};

#endif  // UTILS_PERSISTENCE_COMPUTE_H_
