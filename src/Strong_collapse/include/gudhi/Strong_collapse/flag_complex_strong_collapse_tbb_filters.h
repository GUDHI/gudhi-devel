/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
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

#ifndef STRONG_COLLAPSE_FLAG_COMPLEX_STRONG_COLLAPSE_TBB_FILTERS_H_
#define STRONG_COLLAPSE_FLAG_COMPLEX_STRONG_COLLAPSE_TBB_FILTERS_H_

// Not available if TBB is not used
#ifdef GUDHI_USE_TBB

#include <gudhi/Strong_collapse/Flag_complex_sparse_matrix.h>
#include <gudhi/Strong_collapse/Flag_complex_tower_assembler.h>

#include <tbb/pipeline.h>  // for tbb::pipeline
#include <tbb/concurrent_vector.h>  // for tbb::concurrent_vector

#include <cstddef>  // for std::size_t
#include <mutex>  // for std::mutex
#include <iomanip>  // for std::setw

namespace Gudhi {

namespace strong_collapse {

using Indexed_sparse_matrix_ptr = std::pair<std::size_t, Flag_complex_sparse_matrix*>;

template<typename FilteredEdgesVector, class InputStepRange = std::initializer_list<double>>
class Strong_collapse_parallel_filter: public tbb::filter {
public:
  Strong_collapse_parallel_filter(std::size_t number_of_points, const FilteredEdgesVector& edge_graph) :
      filter(parallel),
      number_of_points_(number_of_points),
      edge_graph_(edge_graph),
      index_(0) { }

  Strong_collapse_parallel_filter(std::size_t number_of_points, const FilteredEdgesVector& edge_graph,
                                  const InputStepRange& step_range) :
      filter(parallel),
      number_of_points_(number_of_points),
      edge_graph_(edge_graph),
      step_range_(step_range),
      index_(0) { }

private:
  std::size_t number_of_points_;
  FilteredEdgesVector edge_graph_;
  InputStepRange step_range_;
  std::size_t index_;
  std::mutex index_mutex;

  void *operator()(void *) {
    index_mutex.lock();
    std::size_t index_backup = index_;
    index_++;
    index_mutex.unlock();
    Flag_complex_sparse_matrix *collapsed_matrix_ptr = nullptr;

    if (step_range_.size() > 0) {
      // Approximate version
      if (index_backup >= step_range_.size()) {
        // Parallel end condition
        return nullptr;
      }
      double threshold = *(step_range_.begin() + index_backup);
      collapsed_matrix_ptr = new Flag_complex_sparse_matrix(number_of_points_,
                                                            edge_graph_.sub_filter_edges_by_filtration(threshold));
    } else {
      // Exact version
      if (index_backup > edge_graph_.size()) {
        // Parallel end condition
        return nullptr;
      }
      collapsed_matrix_ptr = new Flag_complex_sparse_matrix(number_of_points_,
                                                            edge_graph_.sub_filter_edges_by_index(index_backup));
    }

    collapsed_matrix_ptr->strong_collapse();

    Indexed_sparse_matrix_ptr* indexed_matrix_ptr = new Indexed_sparse_matrix_ptr(index_backup, collapsed_matrix_ptr);
    return indexed_matrix_ptr;
  }
};

template<typename FilteredEdgesVector, class InputStepRange = std::initializer_list<double>>
class Tower_assembler_parallel_filter: public tbb::filter {
public:
  Tower_assembler_parallel_filter(std::size_t number_of_points, const FilteredEdgesVector& edge_graph,
                                        Flag_complex_tower_assembler* tower_assembler_ptr) :
      filter(serial_in_order),
      edge_graph_(edge_graph),
      tower_assembler_ptr_(tower_assembler_ptr),
      parallel_matrix_ptr_(edge_graph_.size(), nullptr),
      index_(0) {
    matrix_before_collapse_ptr_ = new Flag_complex_sparse_matrix(number_of_points);
  }

  Tower_assembler_parallel_filter(std::size_t number_of_points, const FilteredEdgesVector& edge_graph,
                                  const InputStepRange& step_range, Flag_complex_tower_assembler* tower_assembler_ptr) :
      filter(parallel),
      edge_graph_(edge_graph),
      step_range_(step_range),
      tower_assembler_ptr_(tower_assembler_ptr),
      parallel_matrix_ptr_(edge_graph_.size(), nullptr),
      index_(0) {
    matrix_before_collapse_ptr_ = new Flag_complex_sparse_matrix(number_of_points);
  }

  ~Tower_assembler_parallel_filter() {
    delete matrix_before_collapse_ptr_;
  }

private:
  FilteredEdgesVector edge_graph_;
  InputStepRange step_range_;
  Flag_complex_tower_assembler* tower_assembler_ptr_;
  tbb::concurrent_vector<Flag_complex_sparse_matrix*> parallel_matrix_ptr_;
  std::size_t index_;
  Flag_complex_sparse_matrix* matrix_before_collapse_ptr_;
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
          double threshold = *(step_range_.begin() + index_);

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