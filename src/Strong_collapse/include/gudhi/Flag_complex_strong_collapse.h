/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siddharth Pritam
 *
 *    Copyright (C) 2019 Inria
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

#include <gudhi/Strong_collapse/Flag_complex_sparse_matrix.h>
#include <gudhi/Strong_collapse/Flag_complex_tower_assembler.h>
#include <gudhi/graph_simplicial_complex.h>

#ifdef GUDHI_USE_TBB
#include <gudhi/Strong_collapse/flag_complex_strong_collapse_tbb_filters.h>

#include <tbb/task_scheduler_init.h>  // for tbb::task_scheduler_init
#include <tbb/pipeline.h>  // for tbb::pipeline
#endif  // GUDHI_USE_TBB

#include <set>
#include <fstream>
#include <string>
#include <limits>
#include <algorithm>

namespace Gudhi {

namespace strong_collapse {

/**
 * \class Flag_complex_strong_collapse
 * \brief Flag complex for strong collapse data structure.
 *
 * \ingroup strong_collapse
 *
 * \details
 * Strong collapse is constructing in loop Flag complexes, stored in a `Flag_complex_sparse_matrix`, collapse them,
 * and finally fill a Tower assembler, aka. `Flag_complex_tower_assembler`.
 *
 * \remark `Flag_complex_strong_collapse` is parallelized using \ref tbb.
 * Performances are quite improved when it is installed.
 */
class Flag_complex_strong_collapse {
 private:
  std::size_t number_of_points_;
  Gudhi::strong_collapse::Flag_complex_tower_assembler tower_assembler_;
 public:
  /** \brief Flag_complex_strong_collapse constructor.
   *
   * Initializes the Tower assembler, aka. `Flag_complex_tower_assembler`.
   *
   * @param[in] number_of_points Number of points.
   */
  Flag_complex_strong_collapse(std::size_t number_of_points)
      : number_of_points_(number_of_points)
      , tower_assembler_(number_of_points) {}

  /** \brief Strong collapse approximate version initialization.
   *
   * It constructs in loop Flag complexes but only for the given filtration values given, collapse them,
   * and finally fill a Tower assembler, aka. `Flag_complex_tower_assembler`.
   *
   * @param[in] edge_graph A valid Gudhi::Filtered_edges_vector.
   * @param[in] step_range A valid range of filtration values.
   * @exception std::invalid_argument In case step_range is empty.
   */
  template<typename SimplicialComplex,
      class InputStepRange = std::initializer_list<typename SimplicialComplex::Filtration_value>>
  void initialize_approximate_version(const Gudhi::Filtered_edges_vector<SimplicialComplex>& edge_graph,
                                      const InputStepRange& step_range) {
    GUDHI_CHECK(std::begin(step_range) != std::end(step_range),
                std::invalid_argument("At least one step_range is mandatory for initialize_approximate_version"));

    // Copy step range to be modified
    InputStepRange step_range_copy(step_range.begin(), step_range.end());
    // edge graph min length is required in the range
    step_range_copy.insert(step_range_copy.begin(), edge_graph.get_filtration_min());
    // If we want to go further edge graph max length, let's stop at edge graph max length
    if (step_range_copy.back() > edge_graph.get_filtration_max())
      step_range_copy.push_back(edge_graph.get_filtration_max());

    // Insert by default min and max of edge graph to ease the interface
    std::sort(step_range_copy.begin(), step_range_copy.end());
    // Remove duplicate values
    step_range_copy.erase(std::unique(step_range_copy.begin(), step_range_copy.end()), step_range_copy.end());
    // Remove all thresholds values that are < edge graph min length and > edge graph max length
    step_range_copy.erase(std::remove_if(step_range_copy.begin(), step_range_copy.end(), [edge_graph](const double& x) {
                                           return (x < edge_graph.get_filtration_min() ||
                                               x > edge_graph.get_filtration_max());
                                         }), step_range_copy.end());

#ifdef GUDHI_USE_TBB

    tbb::task_scheduler_init init_parallel;

    // Create the pipeline
    tbb::pipeline pipeline;

    using Filtered_edges_vector = Gudhi::Filtered_edges_vector<SimplicialComplex>;
    // Create strong collapse parallel stage and add it to the pipeline
    Strong_collapse_parallel_filter<Filtered_edges_vector, InputStepRange> collapse_parallel(number_of_points_,
                                                                                             edge_graph,
                                                                                             step_range_copy);
    pipeline.add_filter(collapse_parallel);

    // Create tower assembler parallel stage and add it to the pipeline
    Tower_assembler_parallel_filter<Filtered_edges_vector, InputStepRange> tower_parallel(number_of_points_,
                                                                                          edge_graph,
                                                                                          step_range_copy,
                                                                                          &tower_assembler_);
    pipeline.add_filter(tower_parallel);

    // Run the pipeline
    // Need more than one token in flight per thread to keep all threads
    // busy; 2-4 works
    pipeline.run(init_parallel.default_num_threads() * 4);

#else  // GUDHI_USE_TBB

    Flag_complex_sparse_matrix matrix_before_collapse(number_of_points_);

    for (auto threshold : step_range_copy) {
#ifdef DEBUG_TRACES
      std::cout << "Flag_complex_strong_collapse::initialize_approximate_version - threshold=" << threshold << std::endl;
#endif  // DEBUG_TRACES
      Flag_complex_sparse_matrix collapsed_matrix(number_of_points_,
                                                  edge_graph.sub_filter_edges_by_filtration(threshold));

      collapsed_matrix.strong_collapse();
      tower_assembler_.build_tower_for_two_complexes(matrix_before_collapse, collapsed_matrix, threshold);
      matrix_before_collapse = collapsed_matrix;
    }

#endif  // GUDHI_USE_TBB

  }

  /** \brief Strong collapse exact version initialization.
   *
   * It constructs in loop Flag complexes for every filtration values in edge_graph, collapse them,
   * and finally fill a Tower assembler, aka. `Flag_complex_tower_assembler`.
   *
   * @param[in] edge_graph A valid Gudhi::Filtered_edges_vector.
   */
  template<typename SimplicialComplex>
  void initialize_exact_version(const Gudhi::Filtered_edges_vector<SimplicialComplex>& edge_graph) {
#ifdef GUDHI_USE_TBB

    tbb::task_scheduler_init init_parallel;

    // Create the pipeline
    tbb::pipeline pipeline;

    using Filtered_edges_vector = Gudhi::Filtered_edges_vector<SimplicialComplex>;
    // Create strong collapse parallel stage and add it to the pipeline
    Strong_collapse_parallel_filter<Filtered_edges_vector> collapse_parallel(number_of_points_, edge_graph);
    pipeline.add_filter(collapse_parallel);

    // Create tower assembler parallel stage and add it to the pipeline
    Tower_assembler_parallel_filter<Filtered_edges_vector> tower_parallel(number_of_points_, edge_graph,
                                                                          &tower_assembler_);
    pipeline.add_filter(tower_parallel);

    // Run the pipeline
    // Need more than one token in flight per thread to keep all threads
    // busy; 2-4 works
    pipeline.run(init_parallel.default_num_threads() * 4);

#else  // GUDHI_USE_TBB

    Flag_complex_sparse_matrix matrix_before_collapse(number_of_points_);

    for (std::size_t index = 0; index < edge_graph.size(); index++) {
      Flag_complex_sparse_matrix collapsed_matrix(number_of_points_,
                                                  edge_graph.sub_filter_edges_by_index(index));

      collapsed_matrix.strong_collapse();
      tower_assembler_.build_tower_for_two_complexes(matrix_before_collapse, collapsed_matrix,
                                                     edge_graph.get_filtration_at(index));
      matrix_before_collapse = collapsed_matrix;
    }

#endif  // GUDHI_USE_TBB

  }

  /** \brief Returns the distance matrix constructed from the Tower assembler, aka. `Flag_complex_tower_assembler`.
   *
   * @return Distance matrix.
   */
  Distance_matrix get_distance_matrix() {
    return tower_assembler_.distance_matrix();
  }

};

}  // namespace strong_collapse

}  // namespace Gudhi

#endif  // FLAG_COMPLEX_STRONG_COLLAPSE_H_
