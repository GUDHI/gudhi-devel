/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Loiseaux, Hannah Schreiber
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file multi_persistence_landscapes.h
 * @author David Loiseaux, Hannah Schreiber
 * @brief Contains the methods @ref Gudhi::multi_persistence::build_complex_from_scc_file,
 * @ref Gudhi::multi_persistence::write_complex_to_scc_file, @ref Gudhi::multi_persistence::build_slicer_from_scc_file,
 * @ref Gudhi::multi_persistence::build_complex_from_bitmap and @ref Gudhi::multi_persistence::build_slicer_from_bitmap.
 */

#ifndef MP_LANDSCAPES_HELPERS_H_
#define MP_LANDSCAPES_HELPERS_H_

#include <algorithm>
#include <array>
#include <cstddef>
#include <utility>
#include <vector>

#ifdef GUDHI_USE_TBB
#include <oneapi/tbb/enumerable_thread_specific.h>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/task_arena.h>
#endif

#include <gudhi/Debug_utils.h>
#include <gudhi/Multi_persistence/Box.h>
#include <gudhi/Multi_persistence/Module.h>
#include <gudhi/Multi_persistence/Summand.h>
#include <gudhi/Multi_persistence/summand_helpers.h>
#include <gudhi/Multi_persistence/utils.h>

namespace Gudhi {
namespace multi_persistence {

namespace detail {

/**
 * @ingroup multi_persistence
 * @private
 */
template <typename T>
inline void _update_top_values(T value, std::vector<T>& top) {
  if (top.empty() || !(value > top.back())) return;

  // if top.size() is very large, it could be even worth it to use std::pop_heap/std::push_heap instead
  std::size_t pos = top.size() - 1;
  while (pos > 0 && top[pos - 1] < value) {
    top[pos] = top[pos - 1];
    --pos;
  }
  top[pos] = value;
}

/**
 * @ingroup multi_persistence
 * @private
 */
template <typename T>
inline void _get_top_values(const std::vector<std::array<T, 2>>& bars, double t, std::vector<double>& top) {
  std::fill(top.begin(), top.end(), 0.0);
  for (const auto& bar : bars) {
    const double value = std::max(0.0, std::min(t - static_cast<double>(bar[0]), static_cast<double>(bar[1]) - t));
    _update_top_values(value, top);
  }
}

/**
 * @ingroup multi_persistence
 * @private
 */
template <class Slicer, typename T>
inline std::vector<std::array<T, 2>> _compute_barcode_on_line(Slicer& slicer, T xBasePoint, T yBasePoint, T xDirection,
                                                              T yDirection, int degree, bool initialize,
                                                              bool ignoreInf) {
  slicer.push_to(Line<T>({xBasePoint, yBasePoint}, {xDirection, yDirection}));

  if (initialize) {
    slicer.initialize_persistence_computation(ignoreInf);
  } else {
    // when non-vine, does same then initialize_persistence_computation
    slicer.update_persistence_computation(ignoreInf);
  }
  auto barcode = slicer.template get_flat_barcode<true, T, false>();
  GUDHI_CHECK(degree >= 0 && static_cast<std::size_t>(degree) < barcode.size(),
              std::out_of_range("Landscape degree is outside barcode degree range."));
  // if the slicer contains more than just dimension degree and degree + 1, this feels quite wasteful...
  return barcode[degree];
}

/**
 * @ingroup multi_persistence
 * @private
 */
template <typename T>
inline std::vector<std::pair<typename Module<T>::Summand_t const*, Box<T>>> _build_module_landscape_summand_cache(
    const Module<T>& module, typename Module<T>::Dimension dimension) {
  using Summand = typename Module<T>::Summand_t;

  std::vector<std::pair<Summand const*, Box<T>>> summands;
  summands.reserve(module.size());
  for (const auto& summand : module) {
    if (summand.get_number_of_death_corners() != 0 && summand.get_number_of_birth_corners() != 0 &&
        summand.get_dimension() == dimension)
      summands.emplace_back(&summand, summand.compute_bounds());
  }
  return summands;
}

/**
 * @ingroup multi_persistence
 * @private
 */
template <typename T, class RandomAccessValueRange>
inline bool _could_have_positive_module_landscape(
    const std::pair<typename Module<T>::Summand_t const*, Box<T>>& summand, const RandomAccessValueRange& x) {
  if (x.size() != summand.second.get_number_of_coordinates()) return true;
  const auto& lower = summand.second.get_lower_corner();
  const auto& upper = summand.second.get_upper_corner();
  for (std::size_t parameter = 0; parameter < lower.size(); ++parameter) {
    if (x[parameter] <= lower[parameter] || x[parameter] >= upper[parameter]) return false;
  }
  return true;
}

/**
 * @ingroup multi_persistence
 * @private
 */
template <class RandomAccessValueRange>
inline std::size_t _module_landscape_top_size(const RandomAccessValueRange& ks) {
  auto maxKValue = *std::max_element(ks.begin(), ks.end());
  GUDHI_CHECK(maxKValue >= 0, std::invalid_argument("Landscape ks must be positive"));
  GUDHI_CHECK(static_cast<std::size_t>(maxKValue) < std::numeric_limits<std::size_t>::max(),
              std::length_error("Landscape index is too large."));
  return static_cast<std::size_t>(maxKValue) + 1;
}

/**
 * @ingroup multi_persistence
 * @private
 */
template <typename T, typename D, class RandomAccessValueRange>
inline auto _compute_summand_landscape_value(const Summand<T, D>& sum, const RandomAccessValueRange& x) {
  using signedT = maybe_make_signed_t<T>;
  signedT landscapeValue = 0;
  // TODO: if the types of Births and Deaths in Summand changes (to become a template for example)
  // the input to _get_summand_diagonal has to get adapted to it, as it makes use of
  // Underlying_container::value_type working like a vector
  const auto &births = sum.get_upset().get_underlying_container();
  const auto &deaths = sum.get_downset().get_underlying_container();
  for (const auto &birth : births) {
    for (const auto& death : deaths) {
      signedT value = std::min(_get_summand_diagonal<signedT>(birth, x), _get_summand_diagonal<signedT>(x, death));
      landscapeValue = std::max(landscapeValue, value);
    }
  }
  return landscapeValue;
}

/**
 * @ingroup multi_persistence
 * @private
 */
template <typename T, class RandomAccessValueRange1, class RandomAccessValueRange2>
inline void _set_module_landscape_pixel(
    std::vector<maybe_make_signed_t<T>>& images, std::size_t pixel, std::size_t planeSize,
    const std::vector<std::pair<typename Module<T>::Summand_t const*, Box<T>>>& summands,
    const RandomAccessValueRange1& x, const RandomAccessValueRange2& ks, std::vector<maybe_make_signed_t<T>>& top) {
  using SignedT = maybe_make_signed_t<T>;

  std::fill(top.begin(), top.end(), SignedT(0));
  for (const auto& summand : summands) {
    if (_could_have_positive_module_landscape(summand, x))
      _update_top_values(_compute_summand_landscape_value(*summand.first, x), top);
  }
  for (std::size_t index = 0; index < ks.size(); ++index) {
    const auto k = static_cast<std::size_t>(ks[index]);
    images[index * planeSize + pixel] = k < top.size() ? top[k] : SignedT(0);
  }
}

}  // namespace detail

/**
 * @ingroup multi_persistence
 *
 * @brief TODO
 *
 * @tparam Slicer Either @ref Slicer or @ref Thread_safe_slicer class with any valid template combination.
 * @tparam GridAxisRange Range with size() and operator[] method returning a value convertible into
 * @ref Slicer::value_type.
 * @tparam DirectionRange Range with size() and operator[] method returning a value convertible into
 * @ref Slicer::value_type.
 * @tparam IndexRange Integer range with begin(), end(), size() and operator[] method returning a value convertible
 * into `std::size_t`.
 * @param mainSlicer
 * @param xGrid
 * @param yGrid
 * @param direction
 * @param xStride
 * @param yStride
 * @param dt
 * @param degree
 * @param ks
 * @param ignoreInf
 * @param n_jobs If TBB is linked, allows to specify the number of threads that should be used for parallelization.
 */
template <class Slicer, class GridAxisRange, class DirectionRange, class IndexRange>
inline std::vector<double> compute_slicer_landscapes_on_grid(Slicer& mainSlicer, const GridAxisRange& xGrid,
                                                             const GridAxisRange& yGrid,
                                                             const DirectionRange& direction, std::size_t xStride,
                                                             std::size_t yStride, double dt, int degree,
                                                             const IndexRange& ks, bool ignoreInf,
                                                             [[maybe_unused]] int n_jobs = 0) {
  using T = typename Slicer::value_type;

  const std::size_t xSize = xGrid.size();
  const std::size_t ySize = yGrid.size();

  GUDHI_CHECK(xGrid.size() > 0 && yGrid.size() > 0, std::invalid_argument("Grid axis have to be non empty."));
  GUDHI_CHECK(direction.size() == 2, std::invalid_argument("Direction has to be 2-dimensional."));
  GUDHI_CHECK(xStride > 0 && yStride > 0, std::invalid_argument("Grid strides have to be strictly positive."));
  GUDHI_CHECK(std::isfinite(dt) && dt > 0, std::invalid_argument("Grid step has to be finite and strictly positive."));
  GUDHI_CHECK(xGrid.size() < std::numeric_limits<std::size_t>::max() / yGrid.size(),
              std::invalid_argument("Grid is too large: " + std::to_string(xGrid.size()) + " * " +
                                    std::to_string(yGrid.size()) +
                                    " >= " + std::to_string(std::numeric_limits<std::size_t>::max()) + "."));

  std::vector<double> out(ks.size() * xGrid.size() * yGrid.size(), 0.0);

  if (ks.size() == 0) return out;

  Gudhi::Simple_mdspan view(out.data(), ks.size(), xGrid.size(), yGrid.size());
  const std::size_t numberOfLines =
      (std::min(xStride, xSize) * ySize) + (xStride < xSize ? (xSize - xStride) * std::min(yStride, ySize) : 0);
  auto maxKValue = *std::max_element(ks.begin(), ks.end());

  GUDHI_CHECK(maxKValue >= 0, std::invalid_argument("Landscape ks must be positive"));

  auto get_grid_line_start = [xSize, ySize, xStride, yStride](std::size_t lineIdx) -> std::array<std::size_t, 2> {
    const std::size_t block1_size = std::min(xStride, xSize) * ySize;
    if (lineIdx < block1_size) {
      const std::size_t q = lineIdx / ySize;
      return {q, lineIdx - (q * ySize)};
    }
    const std::size_t idx2 = lineIdx - block1_size;
    const std::size_t b = std::min(yStride, ySize);
    const std::size_t q = idx2 / b;
    return {xStride + q, idx2 - (q * b)};
  };

  auto retrieve_landscape_values = [&](std::size_t x, std::size_t y, const std::vector<std::array<T, 2>>& bars,
                                       double t, std::vector<double>& top) {
    detail::_get_top_values(bars, t, top);
    for (std::size_t k = 0; k < ks.size(); ++k) {
      view(k, x, y) = top[static_cast<std::size_t>(ks[k])];
    }
  };

  auto compute_values_on_line = [&](auto& slicer, std::size_t lineIdx, bool initialize, std::vector<double>& top) {
    const auto [i0, j0] = get_grid_line_start(lineIdx);
    std::vector<std::array<T, 2>> barcode = detail::_compute_barcode_on_line(
        slicer, xGrid[i0], yGrid[j0], direction[0], direction[1], degree, initialize, ignoreInf);
    const std::size_t length = std::min(((xSize - 1 - i0) / xStride) + 1, ((ySize - 1 - j0) / yStride) + 1);
    for (std::size_t step = 0; step < length; ++step) {
      retrieve_landscape_values(i0 + (step * xStride), j0 + (step * yStride), barcode, dt * static_cast<double>(step),
                                top);
    }
  };

  auto compute_value_in_line_range = [&](auto& slicer, std::size_t begin, std::size_t end) {
    bool initialize = true;
    std::vector<double> top(static_cast<std::size_t>(maxKValue) + 1, 0.0);
    for (std::size_t lineIdx = begin; lineIdx < end; ++lineIdx) {
      compute_values_on_line(slicer, lineIdx, initialize, top);
      initialize = false;
    }
  };

#if defined(GUDHI_USE_TBB)
  if (n_jobs == 1) {
    compute_value_in_line_range(mainSlicer, 0, numberOfLines);
  } else {
    const std::size_t targetChunks = n_jobs > 0 ? static_cast<std::size_t>(n_jobs) * 4 : std::size_t(64);
    const std::size_t grainSize = std::max<std::size_t>(1, (numberOfLines + targetChunks - 1) / targetChunks);
    auto run = [&] {
      tbb::parallel_for(tbb::blocked_range<std::size_t>(0, numberOfLines, grainSize), [&](const auto& range) {
        tbb::this_task_arena::isolate([&] {
          auto s = mainSlicer.weak_copy();
          compute_value_in_line_range(s, range.begin(), range.end());
        });
      });
    };
    if (n_jobs > 0) {
      tbb::task_arena arena(n_jobs);
      arena.execute(run);
    } else {
      run();
    }
  }
#else
  compute_value_in_line_range(mainSlicer, 0, numberOfLines);
#endif

  return out;
}

/**
 * @ingroup multi_persistence
 *
 * @brief Computes a set of landscape images for each given `k` (corresponding to the \f$ k^{th} \f$ landscape
 * function).
 * 
 * @tparam T Value type of a parameter in a filtration value.
 * @tparam U Template argument of @ref Box. Has to be either T or std::make_signed_t<T>.
 * @tparam RandomAccessValueRange1 Range of integers with a begin(), end(), size() and operator[] method.
 * @tparam RandomAccessValueRange2 Range of integers with a size() and operator[] method.
 * @param module Module.
 * @param dimension Dimension of the summands to be used.
 * @param ks Range of positive \f$ k \f$'s to compute the \f$ k^{th} \f$ landscape function of the module.
 * @param box Box in which to restrict the landscapes.
 * @param resolution Image resolution. Should have size 2.
 * @param n_jobs If TBB is linked, allows to specify the number of threads that should be used for parallelization.
 * @return A continuous vector of landscape values which should be interpreted as a c-ordered 3-dimensional array with
 * a first axis corresponding the the \f$ k \f$'s, a second axis corresponding to the image axis with the first
 * resolution and a third axis corresponding to the image axis with the second resolution.
 */
template <typename T, typename U, class RandomAccessValueRange1, class RandomAccessValueRange2>
inline std::vector<maybe_make_signed_t<T>> compute_set_of_module_landscapes(
    const Module<T>& module, typename Module<T>::Dimension dimension, const RandomAccessValueRange1& ks,
    const Box<U>& box, const RandomAccessValueRange2& resolution, int n_jobs = 0) {
  static_assert(std::is_same_v<U, T> || std::is_same_v<U, maybe_make_signed_t<T>>,
                "Box template parameter is not compatible with Summand value type.");

  GUDHI_CHECK(resolution.size() >= 2, std::invalid_argument("Not enough resolution values."));

  using SignedT = maybe_make_signed_t<T>;

  const auto nx = static_cast<std::size_t>(resolution[0]);
  const auto ny = static_cast<std::size_t>(resolution[1]);
  const auto planeSize = nx * ny;
  std::vector<SignedT> images(ks.size() * planeSize);
  if (ks.size() == 0 || planeSize == 0) return images;

  const auto summands = detail::_build_module_landscape_summand_cache(module, dimension);
  const auto topSize = std::min(detail::_module_landscape_top_size(ks), summands.size());
  if (topSize == 0) return images;

  const U step_x = (box.get_upper_corner()[0] - box.get_lower_corner()[0]) / static_cast<U>(nx);
  const U step_y = (box.get_upper_corner()[1] - box.get_lower_corner()[1]) / static_cast<U>(ny);
  auto set_pixel = [&](std::size_t pixel, std::vector<SignedT>& top) {
    const auto row = pixel / ny;
    const auto column = pixel - row * ny;
    const std::array<U, 2> x{
        box.get_lower_corner()[0] + step_x * static_cast<U>(row),
        box.get_lower_corner()[1] + step_y * static_cast<U>(column),
    };
    detail::_set_module_landscape_pixel(images, pixel, planeSize, summands, x, ks, top);
  };

#ifdef GUDHI_USE_TBB
  tbb::enumerable_thread_specific<std::vector<SignedT>> top([&] { return std::vector<SignedT>(topSize); });
  oneapi::tbb::task_arena arena(n_jobs);
  arena.execute(
      [&] { tbb::parallel_for(std::size_t(0), planeSize, [&](std::size_t pixel) { set_pixel(pixel, top.local()); }); });
#else
  std::vector<SignedT> top(topSize);
  for (std::size_t pixel = 0; pixel < planeSize; ++pixel) set_pixel(pixel, top);
#endif

  return images;
}

/**
 * @ingroup multi_persistence
 *
 * @brief Computes a set of landscape images for each given `k` (corresponding to the \f$ k^{th} \f$ landscape
 * function).
 * 
 * @tparam T Value type of a parameter in a filtration value.
 * @tparam RandomAccessValueRange Range of integers with a begin(), end(), size() and operator[] method.
 * @tparam RandomAccessArray Range of arithmetic values with a size() and operator[] method.
 * @param module Module.
 * @param dimension Dimension of the summands to be used.
 * @param ks Range of positive \f$ k \f$'s to compute the \f$ k^{th} \f$ landscape function of the module.
 * @param grid Grid partitioning the image. Should have size 2 and the sub-arrays partition an axis of the image each.
 * @param n_jobs If TBB is linked, allows to specify the number of threads that should be used for parallelization.
 * @return A continuous vector of landscape values which should be interpreted as a c-ordered 3-dimensional array with
 * a first axis corresponding the the \f$ k \f$'s, a second axis corresponding to the image axis with the first
 * grid resolution and a third axis corresponding to the image axis with the second grid resolution.
 */
template <typename T, class RandomAccessValueRange, class RandomAccessArray>
inline std::vector<maybe_make_signed_t<T>> compute_set_of_module_landscapes(const Module<T>& module,
                                                                            typename Module<T>::Dimension dimension,
                                                                            const RandomAccessValueRange& ks,
                                                                            const std::vector<RandomAccessArray>& grid,
                                                                            int n_jobs = 0) {
  GUDHI_CHECK(grid.size() >= 2, std::invalid_argument("First axis of the grid has not enough values."));

  using SignedT = maybe_make_signed_t<T>;
  using GridT = std::decay_t<decltype(grid[0][0])>;

  const auto nx = grid[0].size();
  const auto ny = grid[1].size();
  const auto planeSize = nx * ny;
  std::vector<SignedT> images(ks.size() * planeSize);
  if (ks.size() == 0 || planeSize == 0) return images;

  const auto summands = detail::_build_module_landscape_summand_cache(module, dimension);
  const auto topSize = std::min(detail::_module_landscape_top_size(ks), summands.size());
  if (topSize == 0) return images;

  auto set_pixel = [&](std::size_t pixel, std::vector<SignedT>& top) {
    const auto row = pixel / ny;
    const auto column = pixel - row * ny;
    const std::array<GridT, 2> x{grid[0][row], grid[1][column]};
    detail::_set_module_landscape_pixel(images, pixel, planeSize, summands, x, ks, top);
  };

#ifdef GUDHI_USE_TBB
  tbb::enumerable_thread_specific<std::vector<SignedT>> top([&] { return std::vector<SignedT>(topSize); });
  oneapi::tbb::task_arena arena(n_jobs);
  arena.execute(
      [&] { tbb::parallel_for(std::size_t(0), planeSize, [&](std::size_t pixel) { set_pixel(pixel, top.local()); }); });
#else
  std::vector<SignedT> top(topSize);
  for (std::size_t pixel = 0; pixel < planeSize; ++pixel) set_pixel(pixel, top);
#endif

  return images;
}

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // MP_LANDSCAPES_HELPERS_H_