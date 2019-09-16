/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2019 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef STRONG_COLLAPSE_FLAG_COMPLEX_STRONG_COLLAPSE_TBB_FILTERS_H_
#define STRONG_COLLAPSE_FLAG_COMPLEX_STRONG_COLLAPSE_TBB_FILTERS_H_

// Not available if TBB is not used
#ifdef GUDHI_USE_TBB

#include <gudhi/Strong_collapse/Flag_complex_sparse_matrix.h>
#include <gudhi/Strong_collapse/Flag_complex_tower_assembler.h>
#include <gudhi/graph_simplicial_complex.h>

#include <tbb/pipeline.h>  // for tbb::pipeline
#include <tbb/concurrent_vector.h>  // for tbb::concurrent_vector

#include <cstddef>  // for std::size_t
#include <mutex>  // for std::mutex
#include <iomanip>  // for std::setw

namespace Gudhi {

namespace strong_collapse {

template<typename SimplicialComplex, class InputStepRange = std::initializer_list<typename SimplicialComplex::Filtration_value>>
class Strong_collapse_parallel_filter: public tbb::filter {
public:
  Strong_collapse_parallel_filter(std::size_t number_of_points,
                                  const Filtered_edges_vector<SimplicialComplex>& edge_graph) :
      filter(parallel),
      number_of_points_(number_of_points),
      edge_graph_(edge_graph),
      index_(0) { }

  Strong_collapse_parallel_filter(std::size_t number_of_points,
                                  const Filtered_edges_vector<SimplicialComplex>& edge_graph,
                                  const InputStepRange& step_range) :
      filter(parallel),
      number_of_points_(number_of_points),
      edge_graph_(edge_graph),
      step_range_(step_range),
      index_(0) { }

private:
  using Indexed_sparse_matrix_ptr = std::pair<std::size_t, Flag_complex_sparse_matrix<SimplicialComplex>*>;
  std::size_t number_of_points_;
  Filtered_edges_vector<SimplicialComplex> edge_graph_;
  InputStepRange step_range_;
  std::size_t index_;
  std::mutex index_mutex;

  void *operator()(void *) {
    index_mutex.lock();
    std::size_t index_backup = index_;
    index_++;
    index_mutex.unlock();
    Flag_complex_sparse_matrix<SimplicialComplex> *collapsed_matrix_ptr = nullptr;

    if (step_range_.size() > 0) {
      // Approximate version
      if (index_backup >= step_range_.size()) {
        // Parallel end condition
        return nullptr;
      }
      double threshold = *(step_range_.begin() + index_backup);
      collapsed_matrix_ptr = new Flag_complex_sparse_matrix<SimplicialComplex>(number_of_points_,
                                                            edge_graph_.sub_filter_edges_by_filtration(threshold));
    } else {
      // Exact version
      if (index_backup > edge_graph_.size()) {
        // Parallel end condition
        return nullptr;
      }
      collapsed_matrix_ptr = new Flag_complex_sparse_matrix<SimplicialComplex>(number_of_points_,
                                                            edge_graph_.sub_filter_edges_by_index(index_backup));
    }

    collapsed_matrix_ptr->strong_collapse();

    Indexed_sparse_matrix_ptr* indexed_matrix_ptr = new Indexed_sparse_matrix_ptr(index_backup, collapsed_matrix_ptr);
    return indexed_matrix_ptr;
  }
};

template<typename SimplicialComplex, class InputStepRange = std::initializer_list<double>>
class Tower_assembler_parallel_filter: public tbb::filter {
public:
  Tower_assembler_parallel_filter(std::size_t number_of_points,
                                  const Filtered_edges_vector<SimplicialComplex>& edge_graph,
                                  Flag_complex_tower_assembler<SimplicialComplex>* tower_assembler_ptr) :
      filter(serial_in_order),
      edge_graph_(edge_graph),
      tower_assembler_ptr_(tower_assembler_ptr),
      parallel_matrix_ptr_(edge_graph_.size(), nullptr),
      index_(0) {
    matrix_before_collapse_ptr_ = new Flag_complex_sparse_matrix<SimplicialComplex>(number_of_points);
  }

  Tower_assembler_parallel_filter(std::size_t number_of_points,
                                  const Filtered_edges_vector<SimplicialComplex>& edge_graph,
                                  const InputStepRange& step_range,
                                  Flag_complex_tower_assembler<SimplicialComplex>* tower_assembler_ptr) :
      filter(parallel),
      edge_graph_(edge_graph),
      step_range_(step_range),
      tower_assembler_ptr_(tower_assembler_ptr),
      parallel_matrix_ptr_(edge_graph_.size(), nullptr),
      index_(0) {
    matrix_before_collapse_ptr_ = new Flag_complex_sparse_matrix<SimplicialComplex>(number_of_points);
  }

  ~Tower_assembler_parallel_filter() {
    delete matrix_before_collapse_ptr_;
  }

  // Forbid copy/move constructor/assignment operator
  Tower_assembler_parallel_filter(const Tower_assembler_parallel_filter& other) = delete;
  Tower_assembler_parallel_filter& operator= (const Tower_assembler_parallel_filter& other) = delete;
  Tower_assembler_parallel_filter (Tower_assembler_parallel_filter&& other) = delete;
  Tower_assembler_parallel_filter& operator= (Tower_assembler_parallel_filter&& other) = delete;

private:
  using Indexed_sparse_matrix_ptr = std::pair<std::size_t, Flag_complex_sparse_matrix<SimplicialComplex>*>;
  Filtered_edges_vector<SimplicialComplex>  edge_graph_;
  InputStepRange step_range_;
  Flag_complex_tower_assembler<SimplicialComplex>* tower_assembler_ptr_;
  tbb::concurrent_vector<Flag_complex_sparse_matrix<SimplicialComplex>*> parallel_matrix_ptr_;
  std::size_t index_;
  Flag_complex_sparse_matrix<SimplicialComplex>* matrix_before_collapse_ptr_;
  std::mutex index_mutex;

  void *operator()(void * item) {
    // Cast as sender format
    Indexed_sparse_matrix_ptr* indexed_sparse_matrix_ptr = static_cast<Indexed_sparse_matrix_ptr*>(item);

    // Save the Flag_complex_sparse_matrix pointer value
    parallel_matrix_ptr_[indexed_sparse_matrix_ptr->first] = indexed_sparse_matrix_ptr->second;

    index_mutex.lock();
    // Flag_complex_sparse_matrix needs to be treated in order of index
    if (index_ == indexed_sparse_matrix_ptr->first) {
      // Also compute all the Flag_complex_sparse_matrix that arrived before
      while ((parallel_matrix_ptr_[index_] != nullptr) && index_ < parallel_matrix_ptr_.size()) {
        if (step_range_.size() > 0) {
          // Approximate version
          auto threshold = *(step_range_.begin() + index_);

          tower_assembler_ptr_->build_tower_for_two_complexes(*matrix_before_collapse_ptr_,
                                                              *(parallel_matrix_ptr_[index_]),
                                                              threshold);
        } else {
          // Exact version
          tower_assembler_ptr_->build_tower_for_two_complexes(*matrix_before_collapse_ptr_,
                                                              *(parallel_matrix_ptr_[index_]),
                                                              edge_graph_.get_filtration_at(index_));
        }
        delete matrix_before_collapse_ptr_;
        matrix_before_collapse_ptr_ = parallel_matrix_ptr_[index_];
        ++index_;
      }
    }
    index_mutex.unlock();
    delete indexed_sparse_matrix_ptr;
    return nullptr;
  }
};

}  // namespace strong_collapse

}  // namespace Gudhi

#endif  // GUDHI_USE_TBB

#endif  // STRONG_COLLAPSE_FLAG_COMPLEX_STRONG_COLLAPSE_TBB_FILTERS_H_