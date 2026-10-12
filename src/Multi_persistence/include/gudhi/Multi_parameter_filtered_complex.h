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

#include <cstddef>
#include <cstdint>  //std::uint32_t
#include <algorithm>
#include <numeric>
#include <ostream>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#ifdef GUDHI_USE_TBB
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/parallel_sort.h>
#endif

#include <gudhi/Debug_utils.h>
#include <gudhi/serialization_utils.h>
#include <gudhi/Multi_parameter_filtration_value.h>  //for lex order
#include <gudhi/Multi_filtration/multi_filtration_conversions.h>

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
 * @tparam I Index type. Default value: std::uint32_t.
 * @tparam D Dimension type. Default value: int.
 */
template <class MultiFiltrationValue, typename I = std::uint32_t, typename D = int>
class Multi_parameter_filtered_complex {
 public:
  using Index = I;                                    /**< Complex index type. */
  using Filtration_value = MultiFiltrationValue;      /**< Filtration value type. */
  using T = typename Filtration_value::value_type;    /**< Numerical type of an element in a filtration value. */
  using Filtration_value_container = std::vector<Filtration_value>; /**< Filtration value container type. */
  using Boundary = std::vector<Index>;      /**< Cell boundary type, represented by the complex indices of its faces. */
  using Boundary_container = std::vector<Boundary>;   /**< Boundary container type. */
  using Dimension = D;                                /**< Dimension type. */
  using Dimension_container = std::vector<Dimension>; /**< Dimension container type. */

  static constexpr Dimension nullDimension = -1;      /**< Dimension default null value. */

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
  Multi_parameter_filtered_complex(const Boundary_container& boundaries, const Dimension_container& dimensions,
                                   const Filtration_value_container& filtrationValues)
      : boundaries_(boundaries),
        dimensions_(dimensions),
        filtrationValues_(filtrationValues),
        maxDimension_(nullDimension),
        isOrderedByDimension_(false) {
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
  Multi_parameter_filtered_complex(Boundary_container&& boundaries, Dimension_container&& dimensions,
                                   Filtration_value_container&& filtrationValues)
      : boundaries_(std::move(boundaries)),
        dimensions_(std::move(dimensions)),
        filtrationValues_(std::move(filtrationValues)),
        maxDimension_(nullDimension),
        isOrderedByDimension_(false) {
    _initialize_dimension_utils();
  }

  /**
   * @brief Constructs the complex from the given containers.
   *
   * @tparam BoundaryRange Boundary range with a size() and operator[] method. The boundary returned by operator[]
   * has to have a begin() and end() method returning `LegacyForwardIterator` iterators dereferencing into a type
   * convertible into @ref MultiFiltrationValue::value_type.
   * @tparam DimensionRange Dimension range with a begin() and end() method returning `LegacyForwardIterator` iterators
   * dereferencing into a type convertible into `D`.
   * @tparam FiltrationRange Filtration value range with a size() and operator[] method: the value
   * returned can either be a range of a type convertible into @ref MultiFiltrationValue::value_type, or a range of
   * ranges of a type convertible into @ref MultiFiltrationValue::value_type. The middle range also have to have the
   * same methods than the outer range. The very inner range, if existing, has to have a begin() and end() method.
   * @param boundaries Container of boundaries. A boundary has to be described by the indices of its faces in this
   * container. E.g., if a vertex \f$ v \f$ is stored at index \f$ i \f$ and another vertex at index \f$ j \f$, then
   * `boundaries[i]` and `boundaries[j]` are both empty and if the edge \f$ (v,u) \f$ is at index \f$ k \f$, then
   * `boundaries[k]` is equal to `{i, j}`. All boundaries are expected to be ordered by increasing index value.
   * @param dimensions Dimension container. The value at index \f$ i \f$ has to correspond to the dimension of the
   * cell at index \f$ i \f$ in `boundaries`.
   * @param filtrationValues Filtration value container. The value at index \f$ i \f$ has to correspond to the
   * filtration value of the cell at index \f$ i \f$ in `boundaries`.
   */
  template <class BoundaryRange, class DimensionRange, class FiltrationRange>
  Multi_parameter_filtered_complex(const BoundaryRange& boundaries, const DimensionRange& dimensions,
                                   const FiltrationRange& filtrationValues)
      : boundaries_(boundaries.size()),
        dimensions_(dimensions.begin(), dimensions.end()),
        filtrationValues_(),
        maxDimension_(nullDimension),
        isOrderedByDimension_(false) {
    GUDHI_CHECK(boundaries.size() == dimensions.size(),
                std::invalid_argument("There should be as many boundaries as dimensions."));
    GUDHI_CHECK(boundaries.size() == filtrationValues.size(),
                std::invalid_argument("There should be as many boundaries as filtration values."));

    for (std::size_t i = 0; i < boundaries.size(); ++i) {
      const auto& b = boundaries[i];
      boundaries_[i] = Boundary(b.begin(), b.end());
    }
    if (filtrationValues.size() != 0) {
      int numParam = 0;
      const auto& fil = filtrationValues[0];
      if (fil.size() != 0) {
        if constexpr (std::is_arithmetic_v<std::remove_cv_t<std::remove_reference_t<decltype(fil[0])> > >) {
          // 1-critical
          numParam = fil.size();
          filtrationValues_.reserve(filtrationValues.size());
          for (std::size_t i = 0; i < filtrationValues.size(); ++i) {
            const auto& f = filtrationValues[i];  // no for (auto& f: ...) because of Numpy_2D_span
            filtrationValues_.push_back({f.begin(), f.end()});
          }
        } else {
          // k-critical: assuming fil[0] is a range of same sized ranges
          if (fil[0].size() != 0) {
            numParam = fil[0].size();
          }
          filtrationValues_ = Filtration_value_container(filtrationValues.size(), Filtration_value::inf(numParam));
          for (std::size_t i = 0; i < filtrationValues.size(); ++i) {
            const auto& f = filtrationValues[i];
            for (std::size_t g = 0; g < f.size(); ++g) {
              const auto& gen = f[g];
              // TODO: or add_guaranteed_generator?
              filtrationValues_[i].add_generator(gen.begin(), gen.end());
            }
          }
        }
      }
    }

    _initialize_dimension_utils();
  }

  /**
   * @brief Copy constructor.
   */
  Multi_parameter_filtered_complex(const Multi_parameter_filtered_complex& complex) = default;

  /**
   * @brief Copy constructor.
   */
  template <class OtherFiltrationValue, typename OI, typename OD>
  Multi_parameter_filtered_complex(const Multi_parameter_filtered_complex<OtherFiltrationValue, OI, OD>& complex)
      : boundaries_(complex.get_boundaries().begin(), complex.get_boundaries().end()),
        dimensions_(complex.get_dimensions().begin(), complex.get_dimensions().end()),
        filtrationValues_(complex.get_filtration_values().size()),
        maxDimension_(complex.get_max_dimension()),
        isOrderedByDimension_(complex.is_ordered_by_dimension()) {
    const auto& fils = complex.get_filtration_values();
    for (Index i = 0; i < filtrationValues_.size(); ++i) {
      filtrationValues_[i] = fils[i].template as_type<MultiFiltrationValue>();
    }
  }

  /**
   * @brief Move constructor.
   */
  Multi_parameter_filtered_complex(Multi_parameter_filtered_complex&& complex) noexcept = default;

  /**
   * @brief Destructor.
   */
  ~Multi_parameter_filtered_complex() = default;

  /**
   * @brief Assign operator.
   */
  Multi_parameter_filtered_complex& operator=(const Multi_parameter_filtered_complex& other) = default;

  /**
   * @brief Assign operator.
   */
  template <class OtherFiltrationValue, typename OI, typename OD>
  Multi_parameter_filtered_complex& operator=(
      const Multi_parameter_filtered_complex<OtherFiltrationValue, OI, OD>& other) {
    boundaries_ = Boundary_container(other.get_boundaries().begin(), other.get_boundaries().end());
    dimensions_ = Dimension_container(other.get_dimensions().begin(), other.get_dimensions().end());
    const auto& fils = other.get_filtration_values();
    filtrationValues_ = Filtration_value_container(fils.size());
    for (Index i = 0; i < filtrationValues_.size(); ++i) {
      filtrationValues_[i] = fils[i].template as_type<MultiFiltrationValue>();
    }
    maxDimension_ = other.get_max_dimension();
    isOrderedByDimension_ = other.is_ordered_by_dimension();

    return *this;
  }

  /**
   * @brief Move assign operator.
   */
  Multi_parameter_filtered_complex& operator=(Multi_parameter_filtered_complex&& other) noexcept = default;

  [[nodiscard]] bool is_empty() const { return boundaries_.empty(); }

  /**
   * @brief Returns the number of cells in the complex.
   */
  [[nodiscard]] Index get_number_of_cycle_generators() const { return boundaries_.size(); }

  /**
   * @brief Returns the number of parameters in the filtration.
   */
  [[nodiscard]] Index num_parameters() const {
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
  void sort_by_dimension_co_lexicographically() {
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
  void sort(Comp&& comparaison) {
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
        GUDHI_CHECK(curr < i, std::runtime_error("Got curr " + std::to_string(curr) + " >= i " + std::to_string(i)));
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
  Index prune_above_dimension(int maxDim) {
    if (!isOrderedByDimension_) sort_by_dimension_co_lexicographically();
    Index i = 0;
    while (i < dimensions_.size() && dimensions_[i] < maxDim + 1) ++i;
    boundaries_.resize(i);
    dimensions_.resize(i);
    filtrationValues_.resize(i);
    maxDimension_ = dimensions_.empty() ? nullDimension : dimensions_.back();
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
   * @tparam OneDimArray A range of values \f$ U \f$ convertible into `T`. Has to implement
   * a begin, end and operator[] method and a `value_type` definition equal to \f$ U \f$.
   * @param grid Vector of @p OneDimArray with size at least number of filtration parameters. Each array has to be
   * ordered by increasing value.
   * @param coordinate If true, the values are set to the coordinates of the projection in the grid. If false,
   * the values are set to the values at the coordinates of the projection.
   */
  template <typename OneDimArray>
  void coarsen_on_grid(const std::vector<OneDimArray>& grid, bool coordinate = true) {
#ifdef GUDHI_USE_TBB
    tbb::parallel_for(Index(0), Index(filtrationValues_.size()), [&](Index gen) {
      // TODO : preallocate for tbb
      // preallocate what?
      filtrationValues_[gen].project_onto_grid(grid, coordinate);
    });
#else
    for (auto& fil : filtrationValues_) {
      fil.project_onto_grid(grid, coordinate);
    }
#endif
  }

  void make_filtration_non_decreasing() {
    auto order = [&]() -> bool {
      bool modified = false;
      for (std::size_t i = 1; i < boundaries_.size(); ++i) {
        for (auto b : boundaries_[i]) {
          modified |= intersect_lifetimes(filtrationValues_[i], filtrationValues_[b]);
        }
      }
      return modified;
    };

    if (isOrderedByDimension_) {
      order();
      return;
    }

    while (order());  // should always stop
  }

  /**
   * @brief Builds a new complex by reordering the cells in the given complex with the given permutation map.
   */
  friend Multi_parameter_filtered_complex build_permuted_complex(const Multi_parameter_filtered_complex& complex,
                                                                 const std::vector<Index>& permutation) {
    if (permutation.size() > complex.get_number_of_cycle_generators())
      throw std::invalid_argument("Invalid permutation size.");

    const Index nullIndex = -1;
    std::vector<Index> inv(complex.get_number_of_cycle_generators(), nullIndex);
    for (Index i = 0; i < permutation.size(); ++i) inv[permutation[i]] = i;

    Boundary_container newBoundaries;
    newBoundaries.reserve(permutation.size());
    Dimension_container newDimensions;
    newDimensions.reserve(permutation.size());
    Filtration_value_container newFiltrationValues;
    newBoundaries.reserve(permutation.size());

    for (Index i : permutation) {
      Boundary boundary;
      for (Index b : complex.boundaries_[i]) {
        if (inv[b] != nullIndex) boundary.push_back(inv[b]);
      }
      std::sort(boundary.begin(), boundary.end());
      newBoundaries.emplace_back(std::move(boundary));
      newDimensions.push_back(complex.dimensions_[i]);
      newFiltrationValues.emplace_back(complex.filtrationValues_[i]);
    }

    return Multi_parameter_filtered_complex(std::move(newBoundaries), std::move(newDimensions),
                                            std::move(newFiltrationValues));
  }

  /**
   * @brief Builds a new complex by reordering the cells in the given complex the same way than
   * @ref sort_by_dimension_co_lexicographically. Returns a pair with the new complex as first element and the
   * permutation map used as second element.
   */
  friend std::pair<Multi_parameter_filtered_complex, std::vector<Index> > build_permuted_complex(
      const Multi_parameter_filtered_complex& complex) {
    using namespace Gudhi::multi_filtration;

    std::vector<Index> perm(complex.get_number_of_cycle_generators());
    std::iota(perm.begin(), perm.end(), 0);
#ifdef GUDHI_USE_TBB
    tbb::parallel_sort(perm.begin(), perm.end(), [&](Index i, Index j) -> bool {
#else
    std::sort(perm.begin(), perm.end(), [&](Index i, Index j) -> bool {
#endif
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
                                            const std::vector<std::vector<T> >& grid) {
    using namespace Gudhi::multi_filtration;
    using Return_filtration_value = decltype(std::declval<Filtration_value>().template as_type<std::int32_t>());
    using Return_complex = Multi_parameter_filtered_complex<Return_filtration_value, I, D>;

    typename Return_complex::Filtration_value_container coords(complex.get_number_of_cycle_generators());
#ifdef GUDHI_USE_TBB
    tbb::parallel_for(Index(0), Index(coords.size()), [&](Index gen) {
      coords[gen] = compute_coordinates_in_grid<std::int32_t>(complex.filtrationValues_[gen], grid);
    });
#else
    for (Index gen = 0U; gen < coords.size(); ++gen) {
      coords[gen] = compute_coordinates_in_grid<std::int32_t>(complex.filtrationValues_[gen], grid);
    }
#endif
    // all r-value to prevent copy of coords (the other are copied anyway)
    auto b = complex.boundaries_;
    auto d = complex.dimensions_;
    return Return_complex(std::move(b), std::move(d), std::move(coords));
  }

  /**
   * @brief Equality operator.
   */
  friend bool operator==(const Multi_parameter_filtered_complex& a, const Multi_parameter_filtered_complex& b) {
    if (&a == &b) return true;
    if (a.maxDimension_ != b.maxDimension_) return false;
    // TODO: test up to permutation instead ?
    return a.isOrderedByDimension_ == b.isOrderedByDimension_ && a.boundaries_ == b.boundaries_ &&
           a.dimensions_ == b.dimensions_ && a.filtrationValues_ == b.filtrationValues_;
  }

  /**
   * @brief Serialize given value into the buffer at given pointer.
   *
   * @param value Value to serialize.
   * @param start Pointer to the start of the space in the buffer where to store the serialization.
   * @return End position of the serialization in the buffer.
   */
  friend char* serialize_value_to_char_buffer(const Multi_parameter_filtered_complex& value, char* start) {
    char* curr = start;
    curr = serialize_value_to_char_buffer(value.boundaries_, curr);
    curr = serialize_value_to_char_buffer(value.dimensions_, curr);
    curr = serialize_value_to_char_buffer(value.filtrationValues_, curr);
    curr = serialize_value_to_char_buffer(value.maxDimension_, curr);
    curr = serialize_value_to_char_buffer(value.isOrderedByDimension_, curr);
    return curr;
  }

  /**
   * @brief Deserialize the value from a buffer at given pointer and stores it in given value.
   *
   * @param value Value to fill with the deserialized summand.
   * @param start Pointer to the start of the space in the buffer where the serialization is stored.
   * @return End position of the serialization in the buffer.
   */
  friend const char* deserialize_value_from_char_buffer(Multi_parameter_filtered_complex& value, const char* start) {
    const char* curr = start;
    curr = deserialize_value_from_char_buffer(value.boundaries_, curr);
    curr = deserialize_value_from_char_buffer(value.dimensions_, curr);
    curr = deserialize_value_from_char_buffer(value.filtrationValues_, curr);
    curr = deserialize_value_from_char_buffer(value.maxDimension_, curr);
    curr = deserialize_value_from_char_buffer(value.isOrderedByDimension_, curr);
    return curr;
  }

  /**
   * @brief Returns the serialization size of the given summand.
   */
  friend std::size_t get_serialization_size_of(const Multi_parameter_filtered_complex& value) {
    std::size_t size = get_serialization_size_of(value.boundaries_);
    size += get_serialization_size_of(value.dimensions_);
    size += get_serialization_size_of(value.filtrationValues_);
    size += get_serialization_size_of(value.maxDimension_);
    size += get_serialization_size_of(value.isOrderedByDimension_);
    return size;
  }

  /**
   * @brief Outstream operator.
   */
  friend std::ostream& operator<<(std::ostream& stream, const Multi_parameter_filtered_complex& complex) {
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
  void _initialize_dimension_utils() {
    isOrderedByDimension_ = true;
    if (dimensions_.empty()) return;

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
