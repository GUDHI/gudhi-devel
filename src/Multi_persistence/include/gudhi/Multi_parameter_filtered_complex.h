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
 * @file Multi_parameter_filtered_complex.h
 * @author David Loiseaux
 * @brief Contains the @ref Gudhi::multi_persistence::Multi_parameter_filtered_complex class.
 */

#ifndef MP_FILTERED_COMPLEX_H_INCLUDED
#define MP_FILTERED_COMPLEX_H_INCLUDED

#include <cstdint>  //std::uint32_t
#include <algorithm>
#include <numeric>
#include <ostream>
#include <utility>
#include <vector>

#include <gudhi/Debug_utils.h>
#include <gudhi/Multi_parameter_filtration.h>  //for lex order

namespace Gudhi {
namespace multi_persistence {

// TODO: better name
/**
 * @class Multi_parameter_filtered_complex Multi_parameter_filtered_complex.h gudhi/Multi_parameter_filtered_complex.h
 * @ingroup multi_persistence
 *
 * @brief Class storing the boundaries, the dimensions and the filtration values of all cells composing a complex.
 *
 * @tparam MultiFiltrationValue Filtration value class respecting the @ref MultiFiltrationValue concept.
 */
template <class MultiFiltrationValue>
class Multi_parameter_filtered_complex
{
 public:
  using Index = std::uint32_t;                        /**< Complex index type. */
  using Filtration_value = MultiFiltrationValue;      /**< Filtration value type. */
  using T = typename Filtration_value::value_type;    /**< Numerical type of an element in a filtration value. */
  using Filtration_value_container = std::vector<Filtration_value>; /**< Filtration value container type. */
  using Boundary = std::vector<Index>; /**< Cell boundary type, represented by the complex indices of its faces. */
  using Boundary_container = std::vector<Boundary>;   /**< Boundary container type. */
  using Dimension = int;                              /**< Dimension type. */
  using Dimension_container = std::vector<Dimension>; /**< Dimension container type. */

  /**
   * @brief Default constructor. Constructs an empty complex.
   */
  Multi_parameter_filtered_complex() : filtrationValues_(), maxDimension_(-1), isOrderedByDimension_(true) {}

  /**
   * @brief Constructs the complex by copying all three given containers into the class.
   *
   * @param boundaries Container of boundaries. A boundary has to be described by the indices of its faces in this
   * container. E.g., if a vertex \f$ v \f$ is stored at index \f$ i \f$ and another vertex at index \f$ j \f$, then
   * `boundaries[i]` and `boundaries[j]` are both empty and if the edge \f$ (v,u) \f$ is at index \f$ k \f$, then
   * `boundaries[k]` is equal to `{i, j}`. All boundaries are expected to be ordered by increasing index value.
   * @param dimensions Dimension container. The value at index \f$ i \f$ has to correspond to the dimension of the
   * cell at index \f$ i \f$ in `boundaries`.
   * @param filtrationValues Filtration value container. The value at index \f$ i \f$ has to correspond to the
   * filtration value of the cell at index \f$ i \f$ in `boundaries`.
   */
  Multi_parameter_filtered_complex(const Boundary_container& boundaries,
                                   const Dimension_container& dimensions,
                                   const Filtration_value_container& filtrationValues)
      : boundaries_(boundaries),
        dimensions_(dimensions),
        filtrationValues_(filtrationValues),
        maxDimension_(-1),
        isOrderedByDimension_(false)
  {
    _initialize_dimension_utils();
  }

  /**
   * @brief Constructs the complex by moving all three given containers to the class.
   *
   * @param boundaries Container of boundaries. A boundary has to be described by the indices of its faces in this
   * container. E.g., if a vertex \f$ v \f$ is stored at index \f$ i \f$ and another vertex at index \f$ j \f$, then
   * `boundaries[i]` and `boundaries[j]` are both empty and if the edge \f$ (v,u) \f$ is at index \f$ k \f$, then
   * `boundaries[k]` is equal to `{i, j}`. All boundaries are expected to be ordered by increasing index value.
   * @param dimensions Dimension container. The value at index \f$ i \f$ has to correspond to the dimension of the
   * cell at index \f$ i \f$ in `boundaries`.
   * @param filtrationValues Filtration value container. The value at index \f$ i \f$ has to correspond to the
   * filtration value of the cell at index \f$ i \f$ in `boundaries`.
   */
  Multi_parameter_filtered_complex(Boundary_container&& boundaries,
                                   Dimension_container&& dimensions,
                                   Filtration_value_container&& filtrationValues)
      : boundaries_(std::move(boundaries)),
        dimensions_(std::move(dimensions)),
        filtrationValues_(std::move(filtrationValues)),
        maxDimension_(0),
        isOrderedByDimension_(false)
  {
    _initialize_dimension_utils();
  }

  /**
   * @brief Returns the number of cells in the complex.
   */
  [[nodiscard]] Index get_number_of_cycle_generators() const { return boundaries_.size(); }

  /**
   * @brief Returns the number of parameters in the filtration.
   */
  [[nodiscard]] Index get_number_of_parameters() const
  {
    if (filtrationValues_.empty()) return 0;
    return filtrationValues_[0].num_parameters();
  }

  /**
   * @brief Returns true if and only if the boundaries are ordered by dimension. That is, if an index increases,
   * the represented cell at the new index can only have same or higher dimension than the cell at the index before.
   */
  [[nodiscard]] bool is_ordered_by_dimension() const { return isOrderedByDimension_; }

  /**
   * @brief Returns a const reference to the filtration value container.
   */
  const Filtration_value_container& get_filtration_values() const { return filtrationValues_; }

  /**
   * @brief Returns a reference to the filtration value container.
   * @warning The container is not const such that the user can easily modify/update a filtration value. But do not
   * modify the size of the container, its indices have still to correspond to the indices in the other containers.
   */
  Filtration_value_container& get_filtration_values() { return filtrationValues_; }

  /**
   * @brief Returns a const reference to the dimension container.
   */
  [[nodiscard]] const Dimension_container& get_dimensions() const { return dimensions_; }

  /**
   * @brief Returns a const reference to the boundary container.
   */
  [[nodiscard]] const Boundary_container& get_boundaries() const { return boundaries_; }

  /**
   * @brief Returns the maximal dimension of a cell in the complex.
   */
  [[nodiscard]] Dimension get_max_dimension() const { return maxDimension_; }

  /**
   * @brief Sorts the container internally such that the cells are ordered first by dimension and then
   * co-lexicographically by filtration values. If two cells have same dimension and same filtration value, they are
   * considered equal (i.e., they relative position from each other does not matter).
   * Note that the indices of the cells changes therefore.
   */
  void sort_by_dimension_co_lexicographically()
  {
    using namespace Gudhi::multi_filtration;

    sort([&](Index i, Index j) -> bool {
      if (dimensions_[i] == dimensions_[j]) {
        return is_strict_less_than_lexicographically<true>(filtrationValues_[i], filtrationValues_[j]);
      }
      return dimensions_[i] < dimensions_[j];
    });
  }

  /**
   * @brief Sorts the internal containers using the given comparaison method.
   * Note that the indices of the cells changes therefore.
   *
   * @tparam Comp Method type with signature (Index, Index)->bool.
   * @param comparaison Method taking two complex indices (those before the sort) as input and returns true if and
   * only if the cell at the first index is supposed to be placed before the cell at the second index.
   */
  template <typename Comp>
  void sort(Comp&& comparaison)
  {
    // TODO: test if it is not faster to just reconstruct everything instead of swapping
    // Note: perm and inv have to be build in any case
    // if we reconstruct, we additionally build three containers of vector of Index, of Index
    // and of Filtration_value, which will be swapped respectively with boundaries_, dimensions_
    // and filtrationValues_
    // in this version (swapping), we additionally build two containers of Index instead
    // so should theoretically be better, but not so sure if we replace the containers with
    // completely flat containers one day, i.e. with no cheap swap method
    std::vector<Index> perm(boundaries_.size());
    std::iota(perm.begin(), perm.end(), 0);
    std::vector<Index> pos = perm;
    std::vector<Index> invPos = perm;
    std::sort(perm.begin(), perm.end(), std::forward<Comp>(comparaison));
    std::vector<Index> invPerm(boundaries_.size());
    for (Index i = 0; i < perm.size(); ++i) invPerm[perm[i]] = i;

    Dimension lastDim = -1;
    isOrderedByDimension_ = true;

    for (Index curr = 0; curr < perm.size(); ++curr) {
      Index p = perm[curr];
      Index i = pos[p];
      if (i != curr) {
        GUDHI_CHECK(curr < i, "Something is wrong");
        std::swap(boundaries_[curr], boundaries_[i]);
        std::swap(dimensions_[curr], dimensions_[i]);
        swap(filtrationValues_[curr], filtrationValues_[i]);
        std::swap(pos[invPos[curr]], pos[p]);
        std::swap(invPos[curr], invPos[pos[invPos[curr]]]);
      }
      for (Index& b : boundaries_[curr]) b = invPerm[b];
      std::sort(boundaries_[curr].begin(), boundaries_[curr].end());
      if (lastDim > dimensions_[curr]) isOrderedByDimension_ = false;
      lastDim = dimensions_[curr];
    }
  }

  /**
   * @brief Removes completely from the complex all cells of dimension strictly higher than given.
   *
   * @warning If @ref is_ordered_by_dimension does not return true, the complex is sorted by dimension before pruning.
   * So, the indexing changes afterwards.
   *
   * @param maxDim Maximal dimension to keep.
   * @return Number of remaining cells in the complex.
   */
  Index prune_above_dimension(int maxDim)
  {
    if (!isOrderedByDimension_) sort_by_dimension_co_lexicographically();
    Index i = 0;
    while (i < dimensions_.size() && dimensions_[i] < maxDim + 1) ++i;
    boundaries_.resize(i);
    dimensions_.resize(i);
    filtrationValues_.resize(i);
    maxDimension_ = dimensions_.empty() ? -1 : dimensions_.back();
    return i;
  }

  /**
   * @brief Projects all filtration values into the given grid. If @p coordinate is false, the entries are set to
   * the nearest upper bound value with the same parameter in the grid. Otherwise, the entries are set to the indices
   * of those nearest upper bound values.
   * An index \f$ i \f$ of the grid corresponds to the same parameter as the index \f$ i \f$ in a generator of the
   * filtration value. The internal vectors correspond to the possible values of the parameters, ordered by increasing
   * value, forming therefore all together a 2D grid.
   *
   * @param grid Vector of vector with size at least number of filtration parameters.
   * @param coordinate If true, the values are set to the coordinates of the projection in the grid. If false,
   * the values are set to the values at the coordinates of the projection.
   */
  void coarsen_on_grid(const std::vector<std::vector<T> >& grid, bool coordinate = true)
  {
    for (auto gen = 0U; gen < filtrationValues_.size(); ++gen) {
      filtrationValues_[gen].project_onto_grid(grid, coordinate);
    }
  }

  /**
   * @brief Builds a new complex by reordering the cells in the given complex with the given permutation map.
   */
  friend Multi_parameter_filtered_complex build_permuted_complex(const Multi_parameter_filtered_complex& complex,
                                                                 const std::vector<Index>& permutation)
  {
    if (permutation.size() != complex.get_number_of_cycle_generators())
      throw std::invalid_argument("Invalid permutation size.");

    std::vector<Index> inv(permutation.size());
    for (Index i = 0; i < permutation.size(); ++i) inv[permutation[i]] = i;

    Boundary_container newBoundaries;
    newBoundaries.reserve(permutation.size());
    Dimension_container newDimensions;
    newDimensions.reserve(permutation.size());
    Filtration_value_container newFiltrationValues;
    newBoundaries.reserve(permutation.size());

    for (Index i : permutation) {
      Boundary boundary(complex.boundaries_[i]);
      for (Index& b : boundary) b = inv[b];
      std::sort(boundary.begin(), boundary.end());
      newBoundaries.emplace_back(std::move(boundary));
      newDimensions.push_back(complex.dimensions_[i]);
      newFiltrationValues.emplace_back(complex.filtrationValues_[i]);
    }

    return Multi_parameter_filtered_complex(
        std::move(newBoundaries), std::move(newDimensions), std::move(newFiltrationValues));
  }

  /**
   * @brief Builds a new complex by reordering the cells in the given complex the same way than
   * @ref sort_by_dimension_co_lexicographically. Returns a pair with the new complex as first element and the
   * permutation map used as second element.
   */
  friend std::pair<Multi_parameter_filtered_complex, std::vector<Index> > build_permuted_complex(
      const Multi_parameter_filtered_complex& complex)
  {
    using namespace Gudhi::multi_filtration;

    std::vector<Index> perm(complex.get_number_of_cycle_generators());
    std::iota(perm.begin(), perm.end(), 0);
    std::sort(perm.begin(), perm.end(), [&](Index i, Index j) -> bool {
      if (complex.dimensions_[i] == complex.dimensions_[j]) {
        return is_strict_less_than_lexicographically<true>(complex.filtrationValues_[i], complex.filtrationValues_[j]);
      }
      return complex.dimensions_[i] < complex.dimensions_[j];
    });
    auto out = build_permuted_complex(complex, perm);
    return std::make_pair(std::move(out), std::move(perm));
  }

  /**
   * @brief Builds a new complex from the given one by projecting its filtration values on a grid.
   * See @ref coarsen_on_grid with the paramater `coordinate` at true.
   */
  friend auto build_complex_coarsen_on_grid(const Multi_parameter_filtered_complex& complex,
                                            const std::vector<std::vector<T> >& grid)
  {
    using namespace Gudhi::multi_filtration;
    using Return_filtration_value = decltype(std::declval<Filtration_value>().template as_type<std::int32_t>());
    using Return_complex = Multi_parameter_filtered_complex<Return_filtration_value>;

    typename Return_complex::Filtration_value_container coords(complex.get_number_of_cycle_generators());
    for (Index gen = 0U; gen < coords.size(); ++gen) {
      coords[gen] = compute_coordinates_in_grid<std::int32_t>(complex.filtrationValues_[gen], grid);
    }
    return Return_complex(complex.boundaries_, complex.dimensions_, coords);
  }

  // /**
  //  * @brief Compares two boundaries and returns true if and only if the size of the first is strictly smaller than
  //  the
  //  * second, or, if the two sizes are the same, the first is lexicographically strictly smaller than the second.
  //  * The boundaries are assumed to be ordered by increasing values.
  //  */
  // static bool boundary_is_strictly_smaller_than(const Boundary& b1, const Boundary& b2) {
  //   // we want faces to be smaller than proper cofaces
  //   if (b1.size() < b2.size()) return true;
  //   if (b1.size() > b2.size()) return false;

  //   // lexico for others
  //   for (Index i = 0; i < b2.size(); ++i){
  //     if (b1[i] < b2[i]) return true;
  //     if (b1[i] > b2[i]) return false;
  //   }

  //   // equal
  //   return false;
  // }

  /**
   * @brief Outstream operator.
   */
  friend std::ostream& operator<<(std::ostream& stream, const Multi_parameter_filtered_complex& complex)
  {
    stream << "Boundary:\n";
    stream << "{\n";
    for (Index i = 0; i < complex.boundaries_.size(); ++i) {
      const auto& boundary = complex.boundaries_[i];
      stream << i << ": {";
      for (auto b : boundary) stream << b << ", ";
      if (!boundary.empty()) stream << "\b" << "\b ";
      stream << "},\n";
    }
    stream << "}\n";

    stream << "Dimensions: (max " << complex.get_max_dimension() << ")\n";
    stream << "{";
    for (auto d : complex.dimensions_) stream << d << ", ";
    if (!complex.dimensions_.empty()) {
      stream << "\b" << "\b";
    }
    stream << "}\n";

    stream << "Filtration values:\n";
    stream << "{\n";
    for (auto f : complex.filtrationValues_) stream << f << "\n";
    stream << "}\n";

    return stream;
  }

 private:
  Boundary_container boundaries_;               /**< Boundary container. */
  Dimension_container dimensions_;              /**< Dimension container. */
  Filtration_value_container filtrationValues_; /**< Filtration value container. */
  Dimension maxDimension_;                      /**< Maximal dimension of a cell. */
  bool isOrderedByDimension_;                   /**< True if and only if the containers are ordered by dimension. */

  /**
   * @brief Initializes maxDimension_ and isOrderedByDimension_
   */
  void _initialize_dimension_utils()
  {
    isOrderedByDimension_ = true;
    for (Index i = 0; i < dimensions_.size() - 1; ++i) {
      maxDimension_ = std::max(dimensions_[i], maxDimension_);
      if (dimensions_[i] > dimensions_[i + 1]) isOrderedByDimension_ = false;
    }
    maxDimension_ = std::max(dimensions_.back(), maxDimension_);
  }
};

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // MP_FILTERED_COMPLEX_H_INCLUDED
