/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Loiseaux
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - 2025/04 Hannah Schreiber: Reorganization + documentation.
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file Projective_cover_kernel.h
 * @author David Loiseaux
 * @brief Contains the @ref Gudhi::multi_persistence::Projective_cover_kernel class.
 */

#ifndef MP_PROJECTIVE_COVER_KERNEL_H_INCLUDED
#define MP_PROJECTIVE_COVER_KERNEL_H_INCLUDED

#include <stdexcept>  //std::invalid_argument, std::runtime_error
#include <utility>    //std::move
#include <vector>
#include <set>

#include <gudhi/Debug_utils.h>
#include <gudhi/Multi_parameter_filtered_complex.h>
#include <gudhi/Matrix.h>

namespace Gudhi {
namespace multi_persistence {

/**
 * @class Projective_cover_kernel Projective_cover_kernel.h gudhi/Projective_cover_kernel.h
 * @ingroup multi_persistence
 *
 * @brief TODO (what it is + mention that it only works for 2 parameters)
 *
 * @tparam MultiFiltrationValue Filtration value class respecting the @ref MultiFiltrationValue concept.
 * @tparam columnType Column type to use for the matrix used internally. Default value:
 * @ref Gudhi::persistence_matrix::Column_types::NAIVE_VECTOR "NAIVE_VECTOR".
 */
template <class MultiFiltrationValue,
          Gudhi::persistence_matrix::Column_types columnType = Gudhi::persistence_matrix::Column_types::NAIVE_VECTOR>
class Projective_cover_kernel
{
 private:
  /**
   * @brief Options for matrix type.
   */
  struct Matrix_options : Gudhi::persistence_matrix::Default_options<columnType, true> {
    using Index = std::uint32_t;
  };

 public:
  using Filtration_value = MultiFiltrationValue;                                   /**< Filtration value type. */
  using Complex = Multi_parameter_filtered_complex<Filtration_value>;              /**< Complex data type. */
  using Index = typename Complex::Index;                                           /**< Index type. */
  using Dimension = typename Complex::Dimension;                                   /**< Dimension type. */
  using Filtration_value_container = typename Complex::Filtration_value_container; /**< Filt. value container type. */
  using Boundary_container = typename Complex::Boundary_container;                 /**< Boundary container type. */
  using Dimension_container = typename Complex::Dimension_container;               /**< Dimension container type. */
  using Matrix = Gudhi::persistence_matrix::Matrix<Matrix_options>;                /**< Matrix type. */

  // TODO: this only works for 2 parameter modules. Optimize w.r.t. this.
  /**
   * @brief Constructor.
   *
   * @param complex Complex containing the boundaries, dimensions and 2-parameter filtration values necessary,
   * ordered by dimension and then co-lexicographically with respect to the filtration values.
   * @param dim Dimension for which to compute the kernel.
   */
  Projective_cover_kernel(const Complex &complex, Dimension dim)
  {
    using namespace Gudhi::multi_filtration;

    if (complex.get_number_of_parameters() != 2) throw std::invalid_argument("Only available for 2-parameter modules.");
    if (!complex.is_ordered_by_dimension()) throw std::invalid_argument("Complex has to be ordered by dimension.");

    const auto &boundaries = complex.get_boundaries();
    const auto &filtValues = complex.get_filtration_values();
    const auto &dimensions = complex.get_dimensions();

    GUDHI_CHECK_code(
      Index i = 0;
      for (; i < complex.get_number_of_cycle_generators() - 1; ++i) {
        GUDHI_CHECK(filtValues[i].num_generators() == 1,
                    std::invalid_argument("Only available for 1-critical modules."));
        GUDHI_CHECK(dimensions[i] <= dimensions[i + 1],
                    std::invalid_argument("Cells have to be ordered by dimension."));
        if (dimensions[i] == dimensions[i + 1])
          GUDHI_CHECK(is_less_or_equal_than_lexicographically<true>(filtValues[i], filtValues[i + 1]),
                      std::invalid_argument("Cells with same dimension have to be ordered co-lexicographically."));
      }
      GUDHI_CHECK(filtValues[i].num_generators() == 1,
                  std::invalid_argument("Only available for 1-critical modules."));
    )

    Index startDim1, startDim2, end;
    _get_complex_indices(complex, dim, startDim1, startDim2, end);

    Index nberDim = startDim2 - startDim1;
    Index nberDim1 = end - startDim2;
    Index nberGen = nberDim + nberDim1;
    Index shift = startDim1;

    Matrix M(nberGen);
    Matrix N(nberGen);  // slave
    std::vector<bool> isReduced(nberGen);
    // TODO : pivot caches are small : maybe use a flat container instead ?
    std::vector<std::set<int>> pivotCache(nberGen);  // this[pivot] = cols of given pivot

    auto get_pivot = [&M, shift](Index i) { return _get_pivot(M, i, shift); };

    _initialize_matrices(startDim1, end, boundaries, M, N);
    _initialize_containers(nberDim, nberGen, shift, filtValues, dimensions, get_pivot, pivotCache, isReduced);
    SmallQueue lexicoIt = _initialize_lex_queue(nberDim, nberGen, shift, filtValues, pivotCache, get_pivot);

    while (!lexicoIt.empty()) {
      Filtration_value gridValue = lexicoIt.pop();
      for (int i : lexicoIt.get_current_cols()) {
        while (_reduce_column(complex, i, isReduced, pivotCache, M, N, gridValue, get_pivot));
        _update_after_new_pivot(i, lexicoIt, pivotCache, filtValues, gridValue, get_pivot);
      }
    }
  }

  /**
   * @brief Returns the kernel generators. (? TODO)
   */
  Boundary_container build_generators()
  {
    Index start = 0;
    while (dimensions_[start] == dimensions_[0]) ++start;
    return Boundary_container(boundaries_.begin() + start, boundaries_.end());
  }

  /**
   * @brief Returns a complex representing the projective cover kernel.
   */
  Complex create_complex() { return Complex(boundaries_, dimensions_, filtrationValues_); }

 private:
  /**
   * @brief Lexicographically ordered queue. Each filtration value is represented once and contains the list of indices
   * of the columns with this filtration value.
   */
  struct SmallQueue {
    SmallQueue() = default;

    struct MFWrapper {
      MFWrapper(const Filtration_value &g) : g(g) {};

      MFWrapper(const Filtration_value &g, int col) : g(g) { someCols.insert(col); }

      MFWrapper(const Filtration_value &g, std::initializer_list<int> cols) : g(g), someCols(cols.begin(), cols.end())
      {}

      void insert(int col) const { someCols.insert(col); }

      bool operator<(const MFWrapper &other) const { return is_strict_less_than_lexicographically(g, other.g); }

     public:
      Filtration_value g;
      mutable std::set<int> someCols;
    };

    void insert(const Filtration_value &g, int col)
    {
      auto it = queue.find(g);
      if (it != queue.end()) {
        it->insert(col);
      } else {
        queue.emplace(g, col);
      }
    };

    void insert(const Filtration_value &g, const std::initializer_list<int> &cols)
    {
      auto it = queue.find(g);
      if (it != queue.end()) {
        for (int c : cols) it->insert(c);
      } else {
        queue.emplace(g, cols);
      }
    };

    [[nodiscard]] bool empty() const { return queue.empty(); }

    Filtration_value pop()
    {
      if (queue.empty()) throw std::runtime_error("Queue is empty");

      auto out = std::move(*queue.begin());
      queue.erase(queue.begin());
      std::swap(lastCols, out.someCols);
      return out.g;
    }

    const auto &get_current_cols() const { return lastCols; }

   private:
    std::set<MFWrapper> queue;
    std::set<int> lastCols;
  };

  Boundary_container boundaries_;               /** Boundary container. */
  Dimension_container dimensions_;              /** Dimension container. */
  Filtration_value_container filtrationValues_; /** Filtration value container. */

  static void _get_complex_indices(const Complex &complex,
                                   Dimension dim,
                                   Index &startDim1,
                                   Index &startDim2,
                                   Index &end)
  {
    const auto &dims = complex.get_dimensions();
    startDim1 = 0;
    startDim2 = 0;
    end = 0;

    Index i = 0;
    auto size = complex.get_number_of_cycle_generators();
    for (; i < size && dims[i] < dim; ++i);
    if (i == size) throw std::invalid_argument("Given dimension has no generators.");
    startDim1 = i;
    for (; i < size && dims[i] == dim; ++i);
    if (i == size || dims[i] > dim + 1) throw std::invalid_argument("Given dimension has no generators.");
    startDim2 = i;
    for (; i < size && dims[i] == dim + 1; ++i);
    end = i;
  }

  static int _get_pivot(Matrix &M, Index i, Index shift)
  {
    const auto &col = M.get_column(i);
    return col.size() > 0 ? (*col.rbegin()).get_row_index() - shift : -1;
  }

  static void _initialize_matrices(Index start, Index end, const Boundary_container &boundaries, Matrix &M, Matrix &N)
  {
    for (Index i = start; i < end; i++) {
      M.insert_boundary(boundaries[i]);
      N.insert_boundary({i - start});
    }
  }

  template <typename F>
  static SmallQueue _initialize_lex_queue(Index nberDim,
                                          Index nberGen,
                                          Index shift,
                                          const Filtration_value_container &filtValues,
                                          const std::vector<std::set<int>> &pivotCache,
                                          F &&get_pivot)
  {
    SmallQueue lexicoIt;
    for (Index i = nberDim; i < nberGen; ++i) {
      int pivot = std::forward<F>(get_pivot)(i);
      if (pivot >= 0) {
        lexicoIt.insert(filtValues[i + shift], i);
        auto it = pivotCache[pivot].find(i);
        GUDHI_CHECK(it != pivotCache[pivot].end(), std::runtime_error("Column not registered in pivot cache."));
        it++;
        for (; it != pivotCache[pivot].end(); ++it) {
          int colIdx = *it;
          GUDHI_CHECK(static_cast<Index>(colIdx) > i, std::runtime_error("Column not registered in the right order."));
          auto prev = filtValues[colIdx + shift];
          prev.push_to_least_common_upper_bound(filtValues[i + shift]);
          lexicoIt.insert(std::move(prev), colIdx);
        }
      }
    }
    return lexicoIt;
  }

  template <typename F>
  static void _update_after_new_pivot(Index i,
                                      SmallQueue &lexicoIt,
                                      std::vector<std::set<int>> &pivotCache,
                                      const Filtration_value_container &filtValues,
                                      const Filtration_value &gridValue,
                                      F &&get_pivot)
  {
    int pivot = std::forward<F>(get_pivot)(i);

    if (pivot < 0) return;

    auto [it, wasThere] = pivotCache[pivot].insert(i);
    it++;
    for (; it != pivotCache[pivot].end(); ++it) {
      int colIdx = *it;
      GUDHI_CHECK(static_cast<Index>(colIdx) > i,
                  std::runtime_error("(update) Column not registered in the right order."));
      auto prev = filtValues[colIdx];
      if (!(prev >= filtValues[i])) {
        prev.push_to_least_common_upper_bound(filtValues[i]);
        if (is_strict_less_than_lexicographically(gridValue, prev)) {
          lexicoIt.insert(prev, colIdx);
        }
      }
    }
  }

  template <typename F>
  void _initialize_containers(Index nberDim,
                              Index nberGen,
                              Index shift,
                              const Filtration_value_container &filtValues,
                              const Dimension_container &dimensions,
                              F &&get_pivot,
                              std::vector<std::set<int>> &pivotCache,
                              std::vector<bool> &isReduced)
  {
    Index size = (nberGen - nberDim) * 2;
    boundaries_.reserve(size);
    filtrationValues_.reserve(size);
    dimensions_.reserve(size);

    for (Index i = nberDim; i < nberGen; ++i) {
      boundaries_.emplace_back();
      filtrationValues_.push_back(filtValues[i + shift]);
      dimensions_.push_back(dimensions[i + shift]);
      int pivot = std::forward<F>(get_pivot)(i);
      if (pivot < 0) {
        isReduced[i] = true;
      } else {
        auto &currentCache = pivotCache[pivot];
        currentCache.emplace_hint(currentCache.cend(), i);  // j is increasing
      }
    }
  }

  template <typename F>
  bool _reduce_column(const Complex &complex,
                      Index i,
                      std::vector<bool> &isReduced,
                      std::vector<std::set<int>> &pivotCache,
                      Matrix &M,
                      Matrix &N,
                      const Filtration_value &gridValue,
                      F &&get_pivot)
  {
    int pivot = std::forward<F>(get_pivot)(i);
    if (pivot < 0) {
      if (!isReduced[i]) {
        const auto &col = N.get_column(i);
        boundaries_.emplace_back(col.begin(), col.end());
        filtrationValues_.emplace_back(gridValue);
        dimensions_.emplace_back(complex.get_dimensions()[i] + 1);
        isReduced[i] = true;
      }
      return false;
    }
    // WARN : we lazy update variables linked with col i...
    if (pivotCache[pivot].size() == 0) {
      return false;
    }
    const auto &filtValues = complex.get_filtration_values();
    for (Index k : pivotCache[pivot]) {
      if (k >= i) {  // cannot reduce more here. this is a (local) pivot.
        return false;
      }
      if (filtValues[k] <= gridValue) {
        M.add_to(k, i);
        N.add_to(k, i);
        pivotCache[pivot].erase(i);
        // WARN : we update the pivot cache after the update loop
        GUDHI_CHECK(std::forward<F>(get_pivot)(i) < pivot, std::runtime_error("Column addition failed."));
        return true;  // pivot has changed
      }
    }
    return false;  // for loop exhausted (i may not be there because of lazy)
  }
};

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // MP_PROJECTIVE_COVER_KERNEL_H_INCLUDED
