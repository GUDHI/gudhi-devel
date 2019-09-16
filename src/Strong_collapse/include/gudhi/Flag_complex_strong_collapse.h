/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Siddharth Pritam
 *
 *    Copyright (C) 2019 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
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
 * Strong collapse is constructing in loop Flag complexes, stored in a `Flag_complex_sparse_matrix` and collapse them.
 * Finally, one can get the distance matrix of the collapsed complex.
 *
 * \remark `Flag_complex_strong_collapse` is parallelized using \ref tbb.
 * Performances are quite improved when it is installed.
 */
template<class SimplicialComplex>
class Flag_complex_strong_collapse {
 private:
  Gudhi::strong_collapse::Flag_complex_tower_assembler<SimplicialComplex> tower_assembler_;
 public:
  /** \brief Flag_complex_strong_collapse approximate version constructor.
   *
   *
   * It constructs in loop Flag complexes but only for the given filtration values given and collapse them.
   *
   * @param[in] number_of_points Number of points.
   * @param[in] edge_graph A valid Gudhi::Filtered_edges_vector.
   * @param[in] step_range A valid range of filtration values.
   * @exception std::invalid_argument In case step_range is empty.
   */
  template<class InputStepRange = std::initializer_list<typename SimplicialComplex::Filtration_value>>
  Flag_complex_strong_collapse(std::size_t number_of_points,
                               const Gudhi::Filtered_edges_vector<SimplicialComplex>& edge_graph,
                               const InputStepRange& step_range)
      : tower_assembler_(number_of_points) {
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

    // Create strong collapse parallel stage and add it to the pipeline
    Strong_collapse_parallel_filter<SimplicialComplex, InputStepRange> collapse_parallel(number_of_points,
                                                                                         edge_graph,
                                                                                         step_range_copy);
    pipeline.add_filter(collapse_parallel);

    // Create tower assembler parallel stage and add it to the pipeline
    Tower_assembler_parallel_filter<SimplicialComplex, InputStepRange> tower_parallel(number_of_points,
                                                                                      edge_graph,
                                                                                      step_range_copy,
                                                                                      &tower_assembler_);
    pipeline.add_filter(tower_parallel);

    // Run the pipeline
    // Need more than one token in flight per thread to keep all threads
    // busy; 2-4 works
    pipeline.run(init_parallel.default_num_threads() * 4);
#else  // GUDHI_USE_TBB
    Flag_complex_sparse_matrix<SimplicialComplex> matrix_before_collapse(number_of_points);

    for (auto threshold : step_range_copy) {
#ifdef DEBUG_TRACES
      std::cout << "Flag_complex_strong_collapse::initialize_approximate_version - threshold=" << threshold << std::endl;
#endif  // DEBUG_TRACES
      Flag_complex_sparse_matrix<SimplicialComplex> collapsed_matrix(number_of_points,
                                                  edge_graph.sub_filter_edges_by_filtration(threshold));

      collapsed_matrix.strong_collapse();
      tower_assembler_.build_tower_for_two_complexes(matrix_before_collapse, collapsed_matrix, threshold);
      matrix_before_collapse = collapsed_matrix;
    }
#endif  // GUDHI_USE_TBB
  }

  /** \brief Flag_complex_strong_collapse exact version constructor.
   *
   * It constructs in loop Flag complexes for every filtration values in edge_graph and collapse them.
   *
   *
   * @param[in] number_of_points Number of points.
   * @param[in] edge_graph A valid Gudhi::Filtered_edges_vector.
   */
  Flag_complex_strong_collapse(std::size_t number_of_points,
                               const Gudhi::Filtered_edges_vector<SimplicialComplex>& edge_graph)
      : tower_assembler_(number_of_points) {
#ifdef GUDHI_USE_TBB
    tbb::task_scheduler_init init_parallel;

    // Create the pipeline
    tbb::pipeline pipeline;

    // Create strong collapse parallel stage and add it to the pipeline
    Strong_collapse_parallel_filter<SimplicialComplex> collapse_parallel(number_of_points, edge_graph);
    pipeline.add_filter(collapse_parallel);

    // Create tower assembler parallel stage and add it to the pipeline
    Tower_assembler_parallel_filter<SimplicialComplex> tower_parallel(number_of_points, edge_graph,
                                                                          &tower_assembler_);
    pipeline.add_filter(tower_parallel);

    // Run the pipeline
    // Need more than one token in flight per thread to keep all threads
    // busy; 2-4 works
    pipeline.run(init_parallel.default_num_threads() * 4);
#else  // GUDHI_USE_TBB
    Flag_complex_sparse_matrix<SimplicialComplex> matrix_before_collapse(number_of_points);

    for (std::size_t index = 0; index < edge_graph.size(); index++) {
      Flag_complex_sparse_matrix<SimplicialComplex> collapsed_matrix(number_of_points,
                                                  edge_graph.sub_filter_edges_by_index(index));

      collapsed_matrix.strong_collapse();
      tower_assembler_.build_tower_for_two_complexes(matrix_before_collapse, collapsed_matrix,
                                                     edge_graph.get_filtration_at(index));
      matrix_before_collapse = collapsed_matrix;
    }
#endif  // GUDHI_USE_TBB
  }

  /** \brief Returns the distance matrix constructed after strong collapse computation.
   *
   * @return Distance matrix.
   */
  typename Flag_complex_sparse_matrix<SimplicialComplex>::Distance_matrix get_distance_matrix() {
    return tower_assembler_.distance_matrix();
  }

};

}  // namespace strong_collapse

}  // namespace Gudhi

#endif  // FLAG_COMPLEX_STRONG_COLLAPSE_H_
