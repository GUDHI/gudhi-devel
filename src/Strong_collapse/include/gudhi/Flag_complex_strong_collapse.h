/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siddharth Pritam
 *
 *    Copyright (C) 2019 INRIA Sophia Antipolis (France)
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

#ifndef FLAG_COMPLEX_STRONG_COLLAPSE_H_
#define FLAG_COMPLEX_STRONG_COLLAPSE_H_

#include <gudhi/Flag_complex_sparse_matrix.h>
#include <gudhi/Flag_complex_tower_assembler.h>

#include <set>
#include <fstream>
#include <string>
#include <limits>
#include <algorithm>

namespace Gudhi {

namespace strong_collapse {

class Flag_complex_strong_collapse {
 private:
  Gudhi::strong_collapse::Flag_complex_tower_assembler tower_assembler_;
 public:
  template<typename FilteredEdgesVector, class InputStepRange = std::initializer_list<double>>
  Flag_complex_strong_collapse(std::size_t number_of_points,
                               const FilteredEdgesVector& edge_graph,
                               const InputStepRange& step_range)
  : tower_assembler_(number_of_points) {
    GUDHI_CHECK(std::begin(step_range) != std::end(step_range),
                std::invalid_argument("Flag_complex_strong_collapse::ctor - At least one step_range is mandatory"));

    Flag_complex_sparse_matrix matrix_before_collapse(number_of_points);

    // Empty step_range means exact version of the Strong Collapse
    /*if (std::begin(step_range + 1) == std::end(step_range)) {
      Flag_complex_sparse_matrix collapsed_matrix(number_of_points,
                                                  edge_graph.sub_filter_edges_by_filtration(*step_range));

      collapsed_matrix.strong_collapse();
      Reduction_map reduction_matrix = collapsed_matrix.reduction_map();

      tower_assembler_.build_tower_for_two_complexes(matrix_before_collapse, collapsed_matrix, reduction_matrix,
                                                     edge_graph.get_filtration_min());
    } else {*/
#ifdef GUDHI_USE_TBB

#else
      for (auto threshold : step_range) {
        std::cout << "threshold=" << threshold << std::endl;
        Flag_complex_sparse_matrix collapsed_matrix(number_of_points,
                                                    edge_graph.sub_filter_edges_by_filtration(threshold));

        collapsed_matrix.strong_collapse();
        Reduction_map reduction_matrix = collapsed_matrix.reduction_map();

        tower_assembler_.build_tower_for_two_complexes(matrix_before_collapse, collapsed_matrix,
                                                       reduction_matrix, threshold);
        matrix_before_collapse = collapsed_matrix;
      }
#endif
    //}

  }

  Distance_matrix get_distance_matrix() {
    return tower_assembler_.distance_matrix();
  }

};

}  // namespace strong_collapse

}  // namespace Gudhi

#endif  // FLAG_COMPLEX_STRONG_COLLAPSE_H_
