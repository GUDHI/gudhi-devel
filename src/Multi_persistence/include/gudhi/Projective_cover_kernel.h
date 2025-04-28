/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Loiseaux
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - 2025/04 Hannah Schreiber: Reorganization.
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file Projective_cover_kernel.h
 * @author David Loiseaux
 * @brief Contains the @ref Gudhi::multi_persistence::Projective_cover_kernel class.
 */

#ifndef MP_PROJECTIVE_COVER_KERNEL_H_INCLUDED
#define MP_PROJECTIVE_COVER_KERNEL_H_INCLUDED

// #include <utility>  //std::move
// #include <vector>

#include <gudhi/Debug_utils.h>
#include <gudhi/Multi_parameter_filtered_complex.h>
#include <gudhi/Matrix.h>

#include <gudhi/Multi_parameter_filtration.h>

namespace Gudhi {
namespace multi_persistence {

// template <class Multi_filtration_value>
template <Gudhi::persistence_matrix::Column_types columnType = Gudhi::persistence_matrix::Column_types::NAIVE_VECTOR>
class Projective_cover_kernel
{
 private:
  struct Matrix_options : Gudhi::persistence_matrix::Default_options<columnType, true> {
    using Index = std::uint32_t;
  };

 public:
  // using Filtration_value = Multi_filtration_value;
  using Filtration_value = Gudhi::multi_filtration::Multi_parameter_filtration<double>;
  using Complex = Multi_parameter_filtered_complex<Filtration_value>;
  using Index = typename Complex::Index;
  using Dimension = typename Complex::Dimension;
  using Matrix = Gudhi::persistence_matrix::Matrix<Matrix_options>;

  // TODO : this only works for 2 parameter modules. Optimize w.r.t. this.
  // filtration values are assumed to be dim + co-lexicographically sorted
  Projective_cover_kernel(const Complex& complex, Dimension dim) : complex_(&complex)
  {
    if (complex.get_number_of_parameters() != 2) throw std::invalid_argument("Only available for 2-parameter modules.");
    if (complex.filtration_is_one_critical()) throw std::invalid_argument("Only available for 1-critical modules.");

    const auto& boundaries = complex.get_boundaries();
    const auto& filtValues = complex.get_filtration_values();

    Index startDim1, startDim2, end;
    _get_complex_indices(complex, dim, startDim1, startDim2, end);

    Index nberDim = startDim2 - startDim1;
    Index nberDim1 = end - startDim2;
    Index nberGen = nberDim + nberDim1;

    SmallQueue lexico_it;
    Matrix M(nberGen);
    Matrix N(nberGen);  // slave
    for (int i = startDim1; i < end; i++) {
      M.insert_boundary(boundaries[i]);
      N.insert_boundary({i});
    }

    boundaries_.reserve(nberDim1);
    filtrationValues_.reserve(nberDim1);
    dimensions_.reserve(nberDim1);

    auto get_pivot = [&M, startDim1](Index i) { return _get_pivot(M, i, startDim1); };

    std::vector<bool> isReduced(nberGen);
    auto pivotCache = _initialize_pivot_cache(isReduced, nberDim, nberGen, get_pivot);
  }

 private:
  struct SmallQueue {
    SmallQueue() {};

    struct MFWrapper {
      MFWrapper(const Filtration_value &g) : g(g) {};

      MFWrapper(const Filtration_value &g, int col) : g(g) { some_cols.insert(col); }

      MFWrapper(const Filtration_value &g, std::initializer_list<int> cols) : g(g), some_cols(cols.begin(), cols.end())
      {}

      inline void insert(int col) const { some_cols.insert(col); }

      inline bool operator<(const MFWrapper &other) const { return is_strict_less_than_lexicographically(g, other.g); }

     public:
      Filtration_value g;
      mutable std::set<int> some_cols;
    };

    inline void insert(const Filtration_value &g, int col)
    {
      auto it = queue.find(g);
      if (it != queue.end()) {
        it->insert(col);
      } else {
        queue.emplace(g, col);
      }
    };

    inline void insert(const Filtration_value &g, const std::initializer_list<int> &cols)
    {
      auto it = queue.find(g);
      if (it != queue.end()) {
        for (int c : cols) it->insert(c);
      } else {
        queue.emplace(g, cols);
      }
    };

    inline bool empty() const { return queue.empty(); }

    inline Filtration_value pop()
    {
      if (queue.empty()) throw std::runtime_error("Queue is empty");

      auto out = std::move(*queue.begin());
      queue.erase(queue.begin());
      std::swap(last_cols, out.some_cols);
      return out.g;
    }

    const auto &get_current_cols() const { return last_cols; }

   private:
    std::set<MFWrapper> queue;
    std::set<int> last_cols;
  };

  const Complex* complex_;
  std::vector<std::vector<Index> > boundaries_;
  std::vector<Filtration_value> filtrationValues_;
  std::vector<Dimension> dimensions_;

  static void _get_complex_indices(const Complex& complex,
                                   Dimension dim,
                                   Index& startDim1,
                                   Index& startDim2,
                                   Index& end)
  {
    startDim1 = 0;
    startDim2 = 0;
    end = 0;

    Index i = 0;
    auto size = complex.get_number_of_cycle_generators();
    for (; i < size && complex.get_dimension(i) < dim; ++i);
    if (i == size || complex.get_dimension(i) > dim) throw std::invalid_argument("Given dimension has no generators.");
    startDim1 = i;
    for (; i < size && complex.get_dimension(i) == dim; ++i);
    if (i == size || complex.get_dimension(i) > dim + 1)
      throw std::invalid_argument("Given dimension + 1 has no generators.");
    startDim2 = i;
    for (; i < size && complex.get_dimension(i) == dim + 1; ++i);
    end = i;
  }

  static int _get_pivot(const Matrix &M, Index i, Index shift)
  {
    const auto &col = M.get_column(i);
    return col.size() > 0 ? (*col.rbegin()).get_row_index() - shift : -1;
  }

  template <typename F>
  static std::vector<std::set<int>> _initialize_pivot_cache(std::vector<bool> &isReduced,
                                                            Index nberDim,
                                                            Index nberGen,
                                                            F &&get_pivot)
  {
    // TODO : pivot caches are small : maybe use a flat container instead ?
    std::vector<std::set<int>> pivot_cache(nberGen);  // this[pivot] = cols of given pivot (<=nd)
    for (Index j = nberDim; j < nberGen; ++j) {
      int pivot = get_pivot(j);
      if (pivot < 0) {
        isReduced[j] = true;
      } else {
        auto &currentCache = pivot_cache[pivot];
        currentCache.emplace_hint(currentCache.cend(), j);  // j is increasing
      }
    }
  }
};

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // MP_PROJECTIVE_COVER_KERNEL_H_INCLUDED
