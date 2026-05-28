/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Loiseaux
 *
 *    Copyright (C) 2021 Inria
 *
 *    Modification(s):
 *      - 2026/02 Hannah Schreiber: reorganization + small optimizations + documentation
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file Module.h
 * @author David Loiseaux
 * @brief Contains the @ref Gudhi::multi_persistence::Module class.
 */

#ifndef MP_MODULE_H_INCLUDED
#define MP_MODULE_H_INCLUDED

#include <cstddef>
#include <ostream>    //std::ostream
#include <algorithm>  // std::max
#include <stdexcept>
#include <utility>
#include <array>
#include <vector>

#ifdef GUDHI_USE_TBB
#include <oneapi/tbb/parallel_for.h>
#endif
#include <boost/range/any_range.hpp>
#include <boost/range/adaptor/type_erased.hpp>
#include <boost/range/adaptor/filtered.hpp>

#include <gudhi/Debug_utils.h>
#include <gudhi/simple_mdspan.h>
#include <gudhi/Multi_persistence/Box.h>
#include <gudhi/Multi_persistence/Summand.h>

namespace Gudhi {
namespace multi_persistence {

/**
 * @class Module Module.h gudhi/Multi_persistence/Module.h
 * @ingroup multi_persistence
 *
 * @brief Class representing a multi-parameter persistence module as a set of summands represented by their birth
 * and death corners. The module is associated to a box to which the module gets restricted.
 *
 * @tparam T Value type of a parameter in a filtration value.
 */
template <typename T>
class Module {
 public:
  using value_type = T;                                     /**< Value type of filtration value parameter. */
  using Dimension = int;                                    /** Dimension type. */
  using Summand_t = Summand<value_type, Dimension>;         /**< Summand type. */

 private:
  using Module_t = std::vector<Summand_t>;                  /**< Set of summands defining the module. */

 public:
  using iterator = typename Module_t::iterator;             /**< Iterator type. */
  using const_iterator = typename Module_t::const_iterator; /**< Const iterator type. */
  using Index = typename Module_t::size_type;               /**< Summand container indexation type. */
  using Bar = std::array<value_type, 2>;                    /**< Bar type. */
  // using Summand_of_dimension_range =
  //     boost::any_range<Summand_t, boost::forward_traversal_tag, const Summand_t &, std::ptrdiff_t>;

  static constexpr T T_inf = Summand_t::T_inf;     /**< Infinity. */
  static constexpr T T_m_inf = Summand_t::T_m_inf; /**< Minus infinity. */

  /**
   * @brief Default constructor. Builds empty module.
   */
  Module() : module_(), box_(), maxDim_(Summand_t::template get_null_value<Dimension>()) {}

  /**
   * @brief Builds an empty module with given box.
   * 
   * @param box Box restricting the module.
   */
  Module(const Box<value_type> &box) : module_(), box_(box), maxDim_(Summand_t::template get_null_value<Dimension>()) {}

  /**
   * @brief Returns the begin iterator. Iterators are LegacyRandomAccessIterator.
   */
  iterator begin() { return module_.begin(); }

  /**
   * @brief Returns the end iterator.
   */
  iterator end() { return module_.end(); }

  /**
   * @brief Returns the begin const iterator. Iterators are LegacyRandomAccessIterator.
   */
  const_iterator begin() const { return module_.cbegin(); }

  /**
   * @brief Returns the end const iterator.
   */
  const_iterator end() const { return module_.cend(); }

  // TODO: nanobind binding for multipers
  /**
   * @brief Returns a range over all summands of the module with given dimension.
   * 
   * @param dimension Summand dimension.
   */
  /* Summand_of_dimension_range */ auto get_summand_of_dimension_range(Dimension dimension) const {
    auto has_dimension = [&](const Summand_t &sum) { return sum.get_dimension() == dimension; };
    // cython does not handle the adaptor properly without explicitly using boost::adaptors::type_erased
    return module_ |
           boost::adaptors::filtered(has_dimension) /*  |
  boost::adaptors::type_erased<Summand_t, boost::forward_traversal_tag, const Summand_t &, std::ptrdiff_t>() */
        ;
  }

  /**
   * @brief Returns the restricting box.
   */
  const Box<value_type> &get_box() const { return box_; }

  /**
   * @brief Sets the box.
   * 
   * @param box Box restricting the module.
   */
  void set_box(const Box<value_type> &box) { box_ = box; }

  /**
   * @brief Returns a reference to the summand at given index.
   * 
   * @param index Index in summand container.
   *
   * @note If you have to modify the dimension of the summand, you also have to update the maximal
   * dimension of the module with @ref set_max_dimension. Otherwise, there is no guarantee for the
   * correctness of @ref get_max_dimension which is used internally!
   */
  Summand_t &get_summand(Index index) { return module_[index]; }

  /**
   * @brief Returns a const reference to the summand at given index.
   * 
   * @param index Index in summand container.
   */
  const Summand_t &get_summand(Index index) const { return module_[index]; }

  /**
   * @brief Adds a summand to the module at given index.
   * 
   * @param i Index to where the summands has to be added. If the index is already existing in the container,
   * the given summand will replace the summand existing at that index.
   * @param summand Summand to add.
   * @param dimension If specified, sets the dimension of the summand to the given value (which therefore has to be
   * positive). Otherwise, keeps the value as already defined in the summand.
   */
  void add_summand(Index i, const Summand_t &summand,
                   Dimension dimension = Summand_t::template get_null_value<Dimension>()) {
    if (module_.size() <= i) resize(i + 1, summand.get_number_of_parameters());
    module_[i] = summand;
    if (dimension != Summand_t::template get_null_value<Dimension>()) {
      GUDHI_CHECK(dimension >= 0, std::invalid_argument("Summand dimension has to be positive."));
      module_[i].set_dimension(dimension);
    }
    maxDim_ = std::max(maxDim_, module_[i].get_dimension());
  }

  /**
   * @brief Adds a summand at the end of the summand container in the module.
   * 
   * @param summand Summand to add.
   * @param dimension If specified, sets the dimension of the summand to the given value (which therefore has to be
   * positive). Otherwise, keeps the value as already defined in the summand.
   */
  void add_summand(const Summand_t &summand, Dimension dimension = Summand_t::template get_null_value<Dimension>()) {
    add_summand(module_.size(), summand, dimension);
  }

  /**
   * @brief Returns the maximal dimension of the summands in the module.
   */
  [[nodiscard]] Dimension get_max_dimension() const { return maxDim_; }

  /**
   * @brief Sets the maximal dimension of the summands in the module. Has to be used if the dimension of the summands
   * where modified after being added to the module.
   */
  void set_max_dimension(Dimension maxDim) { maxDim_ = maxDim; }

  /**
   * @brief Returns the number of summands in the module.
   */
  [[nodiscard]] Index size() const { return module_.size(); }

  /**
   * @brief Resizes the summand container. If the container contained less than @p size summands, allocates as
   * many trivial summands necessary to fill the count. If the container contained more than @p size summands,
   * the container is truncated to the first @p size summands.
   * 
   * @param size New number of summands.
   * @param numberOfParameters Number of parameters of the new summands.
   */
  void resize(Index size, int numberOfParameters) { module_.resize(size, Summand_t(numberOfParameters)); }

  /**
   * @brief Returns the bounding box of the module independently of the restricting box.
   */
  Box<value_type> compute_bounds() const {
    Dimension numParam = box_.get_lower_corner().size();
    typename Box<value_type>::Point_t lower_bound(numParam, T_inf);
    typename Box<value_type>::Point_t upper_bound(numParam, T_m_inf);
    for (const auto &summand : module_) {
      auto summandBounds = summand.compute_bounds();
      const auto &[m, M] = summandBounds.get_bounding_corners();
      GUDHI_CHECK(m.size() == numParam && M.size() == numParam,
                  std::logic_error("Stored box does not have the same number of parameters than a summand."));
      for (auto parameter = 0; parameter < numParam; parameter++) {
        lower_bound[parameter] = std::min(m[parameter], lower_bound[parameter]);
        upper_bound[parameter] = std::min(M[parameter], upper_bound[parameter]);
      }
    }
    return Box(lower_bound, upper_bound);
  }

  /**
   * @brief Computes the intersection of the (positive slope) line with the summands. That corresponds to the barcode
   * of the 1-dimensional filtration along this line in the module.
   * 
   * @param l Line with positive slope.
   * @param dimension If specified, only the barcode in that dimension is stored in the output. Otherwise,
   * all dimensions are returned.
   * @return A vector of vector of @ref Bar "": the first axis corresponds to a dimension and the second axis
   * to the bars in that dimension. If the argument @p dimension was given, the first axis will still contain all
   * dimensions, but all except one sub-vector will be empty.
   */
  std::vector<std::vector<Bar>> get_barcode_from_line(
      const Line<value_type> &l, Dimension dimension = Summand_t::template get_null_value<Dimension>()) const {
    std::vector<std::vector<Bar>> barcode(get_max_dimension() + 1);
    for (Dimension i = 0; i < get_max_dimension(); ++i) {
      if (dimension == Summand_t::template get_null_value<Dimension>() || i == dimension) {
        barcode[i].reserve(size());
      }
    }
    for (const auto &summand : module_) {
      if (dimension == Summand_t::template get_null_value<Dimension>() || summand.get_dimension() == dimension) {
        barcode[summand.get_dimension()].push_back(summand.get_bar(l));
      }
    }
    return barcode;
  }

  /**
   * @brief Removes from the module all summands with all corners at plus infinity.
   */
  void clean() {
    module_.erase(
        std::remove_if(module_.begin(), module_.end(), [](const Summand_t &s) { return s.get_upset().is_plus_inf(); }),
        module_.end());
    maxDim_ = Summand_t::template get_null_value<Dimension>();
    for (const auto &sum : module_) maxDim_ = std::max(maxDim_, sum.get_dimension());
  }

  /**
   * @brief Identifies all corners of the summands which are close with respect to the given precision value.
   *
   * Two birth/death corners in a same summand are considered close with respect to @p precision if the maximal
   * coordinate difference between both is smaller than @p precision. The new identified corner takes all maximal
   * coordinates of the two, covering therefore both.
   *
   * Note that the merge is done in ordered pairs, so if \f$ b_1 \f$ is close enough to \f$ b_2 \f$ and \f$ b_2 \f$
   * close enough to \f$ b_3 \f$, once \f$ b_1 \f$ and \f$ b_2 \f$ merge, the new point is eventually not merged with
   * \f$ b_3 \f$.
   * 
   * @param precision Distance threshold.
   */
  void fill(value_type precision) {
    // TODO: parallelize
    for (Summand_t &sum : module_) {
      sum.identify_births(precision);
      sum.identify_deaths(precision);
    }
  }

  /**
   * @brief Rescales the corners of all summands for each parameter, that is, for each corner in the summand,
   * the \f$ p^{th} \f$ coordinate is multiplied by `rescaleFactors[p]`.
   * 
   * @tparam RandomAccessValueRange Range with a size() and operator[] method.
   * @param rescaleFactors Rescale factors. There must be at least as many elements than parameters in the summands.
   * @param dimension If specified, only the summands of this dimension are rescaled.
   */
  template <class RandomAccessValueRange>
  void rescale(const RandomAccessValueRange &rescaleFactors,
               Dimension dimension = Summand_t::template get_null_value<Dimension>()) {
    for (auto &summand : module_) {
      if (dimension == Summand_t::template get_null_value<Dimension>() || summand.get_dimension() == dimension)
        summand.rescale(rescaleFactors);
    }
  }

  /**
   * @brief Translates each parameter of the summands in the module, that is, for each corner in the summand,
   * the \f$ p^{th} \f$ coordinate is summed with `translation[p]`.
   *
   * @tparam RandomAccessValueRange Range with a size() and operator[] method.
   * @param translation Translation factors. There must be at least as many elements than parameters in the summand.
   * @param dimension If specified, only the summands of this dimension are translated.
   */
  template <class RandomAccessValueRange>
  void translate(const RandomAccessValueRange &translation,
                 Dimension dimension = Summand_t::template get_null_value<Dimension>()) {
    for (auto &summand : module_) {
      if (dimension == Summand_t::template get_null_value<Dimension>() || summand.get_dimension() == dimension)
        summand.translate(translation);
    }
  }

  /**
   * @brief Snaps the corner coordinates of all summands in the module to their closest integer and interprets them as
   * indices in a grid to replaces them with the values at those indices in the given grid. For example, if \f$ c \f$
   * is a birth or death corner in a summand, its new value at parameter \f$ p \f$ will be
   * \f$ c[p] = grid[p][snap(c[p])] \f$. If \f$ snap(c[p]) \f$ is out of bound, the value is replaced with infinity.
   * 
   * @tparam GridRange Range with size() and operator[] method. The operator[] method must return a type with the
   * same methods and a value type convertible to `T`.
   * @param grid 2-dimensional range with first axis at least as long as number of parameters.
   */
  template <class GridRange>
  void evaluate_in_grid(const GridRange &grid) {
#ifdef GUDHI_USE_TBB
    tbb::parallel_for(std::size_t(0), module_.size(), [&](std::size_t i) { module_[i].evaluate_in_grid(grid); });
#else
    for (auto &summand : module_) {
      summand.evaluate_in_grid(grid);
    }
#endif
  }

  /**
   * @brief Equality operator.
   */
  friend bool operator==(const Module &a, const Module &b) {
    if (a.get_max_dimension() != b.get_max_dimension()) return false;
    if (a.box_ != b.box_) return false;
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
      if (a.get_summand(i) != b.get_summand(i)) return false;
    }
    return true;
  }

  /**
   * @brief Outstream operator.
   */
  friend std::ostream &operator<<(std::ostream &stream, const Module &m) {
    stream << "Module of max dim " << m.maxDim_ << ":\n";
    for (const auto &s : m) {
      stream << s << "\n";
    }
    return stream;
  }

  /**
   * @brief Swap operator.
   */
  friend void swap(Module &mod1, Module &mod2) noexcept {
    mod1.module_.swap(mod2.module_);
    swap(mod1.box_, mod2.box_);
    std::swap(mod1.maxDim_, mod2.maxDim_);
  }

 private:
  Module_t module_;     /**< Summand container. */
  Box<value_type> box_; /**< Restricting box. */
  Dimension maxDim_;    /**< Maximal dimension of the module. */
};

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // MP_MODULE_H_INCLUDED
