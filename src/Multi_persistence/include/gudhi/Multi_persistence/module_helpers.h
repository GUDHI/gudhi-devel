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
 * @file module_helpers.h
 * @author David Loiseaux
 * @brief Contains the helper methods @ref Gudhi::multi_persistence::.
 */

#ifndef MP_MODULE_HELPERS_H_
#define MP_MODULE_HELPERS_H_

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <stdexcept>
#include <type_traits>
#include <vector>

#ifdef GUDHI_USE_TBB
#include <oneapi/tbb/task_arena.h>
#include <oneapi/tbb/parallel_for.h>
#endif

#include <gudhi/Debug_utils.h>
#include <gudhi/simple_mdspan.h>
#include <gudhi/Multi_persistence/Box.h>
#include <gudhi/Multi_persistence/Module.h>
#include <gudhi/Multi_persistence/summand_helpers.h>

namespace Gudhi {
namespace multi_persistence {

/**
 * @ingroup multi_persistence
 *
 * @brief Constructs a new identical module but with permuted summand indices.
 *
 * @tparam T First template parameter of @ref Module.
 * @tparam RandomAccessValueRange Range of intergers with a size() and operator[] method.
 * @param module Module to permute.
 * @param permutation Permutation map of the summand indices, such that the \f$ i^{th} \f$ summand in the new module
 * corresponds to the summand at index \f$ perm[i] \f$ in the original module.
 */
template <typename T, class RandomAccessValueRange>
Module<T> build_permuted_module(const Module<T> &module, const RandomAccessValueRange &permutation) {
  GUDHI_CHECK(permutation.size() <= module.size(),
              std::invalid_argument("Permutation size is greater than the module size."));

  Module<T> out;

  if (module.size() == 0) return out;

  out.resize(module.size(), 1);
  for (std::size_t i = 0; i < permutation.size(); ++i) {
    out.get_summand(i) = module.get_summand(permutation[i]);
  }
  out.set_max_dimension(module.get_max_dimension());
  return out;
}

/**
 * @ingroup multi_persistence
 *
 * @brief Builds a copy of the given module but only keeping the summands of given dimension.
 *
 * @tparam T First template parameter of @ref Module.
 * @param module Module.
 * @param dim Dimension.
 */
template <typename T>
Module<T> build_module_of_dimension(const Module<T> &module, int dim) {
  if (dim < 0 || module.size() == 0) return {};

  Module<T> out;

  // My guess is that iterating twice over the summands is cheaper then having potentially several
  // reallocations because of the push_backs. But I could be wrong? Is it worth benchmarking?
  std::size_t size = 0;
  auto r = module.get_summands_of_dimension_range(dim);
  for (auto it = r.begin(); it != r.end(); ++it) ++size;

  out.resize(size, 1);
  std::size_t i = 0;
  for (const auto &sum : module.get_summands_of_dimension_range(dim)) {
    out.get_summand(i) = sum;
    ++i;
  }
  out.set_max_dimension(dim);
  return out;
}

/**
 * @ingroup multi_persistence
 *
 * @brief Builds a copy of the given module but only keeping the summands of given dimensions.
 *
 * @tparam T First template parameter of @ref Module.
 * @tparam ContinuousRandomAccessRange Continuous (in memory) dimension range with random access iterators.
 * Must have a begin() and end() method.
 * @param module Module.
 * @param dims Range of dimensions corresponding to the dimensions of the summands to copy.
 */
template <typename T, class ContinuousRandomAccessRange,
          class = std::enable_if_t<!std::is_arithmetic_v<ContinuousRandomAccessRange>>>
Module<T> build_module_of_dimension(const Module<T> &module, const ContinuousRandomAccessRange &dims) {
  if (module.size() == 0) return {};

  Module<T> out;

  // Assuming dims will not have more than 10 values, this should be even faster then a std::set::find
  auto contains = [&dims](auto v) { return std::find(dims.begin(), dims.end(), v) != dims.end(); };

  // My guess is that iterating twice over the summands is cheaper then having potentially several
  // reallocations because of the push_backs. But I could be wrong? Is it worth benchmarking?
  std::size_t size = 0;
  for (auto it = module.begin(); it != module.end(); ++it) {
    if (contains(it->get_dimension())) ++size;
  }

  out.resize(size, 1);
  std::size_t i = 0;
  typename Module<T>::Dimension maxDim = -1;
  for (const auto &sum : module) {
    if (contains(sum.get_dimension())) {
      out.get_summand(i) = sum;
      ++i;
      if (sum.get_dimension() > maxDim) maxDim = sum.get_dimension();
    }
  }
  out.set_max_dimension(maxDim);
  return out;
}

/**
 * @ingroup multi_persistence
 *
 * @private
 */
template <typename T, typename U>
inline std::vector<maybe_make_signed_t<T>> _get_module_landscape_values(const Module<T> &mod, const std::vector<U> &x,
                                                                        typename Module<T>::Dimension dimension) {
  using signedT = maybe_make_signed_t<T>;

  std::vector<signedT> values;
  values.reserve(mod.size());
  for (std::size_t i = 0; i < mod.size(); i++) {
    const auto &summand = mod.get_summand(i);
    if (summand.get_dimension() == dimension) values.push_back(compute_summand_landscape_value(summand, x));
  }
  std::sort(values.begin(), values.end(), [](signedT x, signedT y) { return x > y; });
  return values;
}

// TODO: extend in higher resolution dimension
/**
 * @ingroup multi_persistence
 *
 * @brief Computes a set of landscape images for each given `k` (corresponding to the \f$ k^{th} \f$ landscape
 * function).
 *
 * @tparam T Value type of a parameter in a filtration value.
 * @tparam U Template argument of @ref Box. Has to be either T or std::make_signed_t<T>.
 * @tparam RandomAccessValueRange1 Range of unsigned integers with a size() and operator[] method.
 * @tparam RandomAccessValueRange2 Range of unsigned integers with a size() and operator[] method.
 * @param mod Module.
 * @param dimension Dimension of the summands to be used.
 * @param ks Range of \f$ k \f$'s to compute the \f$ k^{th} \f$ landscape function of the module.
 * @param box Box in which to restrict the landscapes.
 * @param resolution Image resolution. Should have size 2.
 * @param n_jobs If TBB is linked, allows to specify the number of threads that should be used for parallelization.
 * @return A continuous vector of landscape values which should be interpreted as a c-ordered 3-dimensional array with
 * a first axis corresponding the the \f$ k \f$'s, a second axis corresponding to the image axis with the first
 * resolution and a third axis corresponding to the image axis with the second resolution.
 */
template <typename T, typename U, class RandomAccessValueRange1, class RandomAccessValueRange2>
inline std::vector<maybe_make_signed_t<T>> compute_set_of_module_landscapes(
    const Module<T> &mod, typename Module<T>::Dimension dimension, const RandomAccessValueRange1 &ks, const Box<U> &box,
    const RandomAccessValueRange2 &resolution, [[maybe_unused]] int n_jobs = 0) {
  static_assert(std::is_same_v<U, T> || std::is_same_v<U, maybe_make_signed_t<T>>,
                "Box template parameter is not compatible with Summand value type.");

  GUDHI_CHECK(resolution.size() >= 2, std::invalid_argument("Not enough resolution values."));

  using signedT = maybe_make_signed_t<T>;

  std::vector<signedT> images(ks.size() * resolution[0] * resolution[1]);
  Simple_mdspan view(images.data(), ks.size(), resolution[0], resolution[1]);
  U stepX = (box.get_upper_corner()[0] - box.get_lower_corner()[0]) / static_cast<U>(resolution[0]);
  U stepY = (box.get_upper_corner()[1] - box.get_lower_corner()[1]) / static_cast<U>(resolution[1]);

  auto get_image_values = [&](unsigned int i) {
    return [&, i](unsigned int j) {
      std::vector<signedT> landscapes =
          _get_module_landscape_values<T, U>(mod,
                                             {box.get_lower_corner()[0] + (stepX * static_cast<U>(i)),
                                              box.get_lower_corner()[1] + (stepY * static_cast<U>(j))},
                                             dimension);
      for (std::size_t k_idx = 0; k_idx < ks.size(); ++k_idx) {
        unsigned int k = ks[k_idx];
        view(k_idx, i, j) = k < landscapes.size() ? landscapes[k] : 0;
      }
    };
  };

#ifdef GUDHI_USE_TBB
  using ResT = typename RandomAccessValueRange2::value_type;
  oneapi::tbb::task_arena arena(n_jobs);
  arena.execute([&] {
    tbb::parallel_for(ResT(0), resolution[0],
                      [&](unsigned int i) { tbb::parallel_for(ResT(0), resolution[1], get_image_values(i)); });
  });
#else
  for (unsigned int i = 0; i < resolution[0]; ++i) {
    auto get_image_values_at = get_image_values(i);
    for (unsigned int j = 0; j < resolution[1]; ++j) {
      get_image_values_at(j);
    }
  }
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
 * @tparam RandomAccessValueRange Range of unsigned integers with a size() and operator[] method.
 * @tparam RandomAccessArray Range of arithmetic values with a size() and operator[] method.
 * @param mod Module.
 * @param dimension Dimension of the summands to be used.
 * @param ks Range of \f$ k \f$'s to compute the \f$ k^{th} \f$ landscape function of the module.
 * @param grid Grid partitioning the image. Should have size 2 and the sub-arrays partition an axis of the image each.
 * @param n_jobs If TBB is linked, allows to specify the number of threads that should be used for parallelization.
 * @return A continuous vector of landscape values which should be interpreted as a c-ordered 3-dimensional array with
 * a first axis corresponding the the \f$ k \f$'s, a second axis corresponding to the image axis with the first
 * grid resolution and a third axis corresponding to the image axis with the second grid resolution.
 */
template <typename T, class RandomAccessValueRange, class RandomAccessArray>
inline std::vector<maybe_make_signed_t<T>> compute_set_of_module_landscapes(const Module<T> &mod,
                                                                            typename Module<T>::Dimension dimension,
                                                                            const RandomAccessValueRange &ks,
                                                                            const std::vector<RandomAccessArray> &grid,
                                                                            [[maybe_unused]] int n_jobs = 0) {
  GUDHI_CHECK(grid.size() >= 2, std::invalid_argument("First axis of the grid has not enough values."));

  if (grid[0].size() == 0 || grid[1].size() == 0) return {};

  using signedT = maybe_make_signed_t<T>;
  using gT = std::decay_t<decltype(grid[0][0])>;

  std::vector<signedT> images(ks.size() * grid[0].size() * grid[1].size());
  Simple_mdspan view(images.data(), ks.size(), grid[0].size(), grid[1].size());

  auto get_image_values = [&](std::size_t i) {
    return [&, i](std::size_t j) {
      std::vector<signedT> landscapes = _get_module_landscape_values<T, gT>(mod, {grid[0][i], grid[1][j]}, dimension);
      for (std::size_t k_idx = 0; k_idx < ks.size(); ++k_idx) {
        unsigned int k = ks[k_idx];
        view(k_idx, i, j) = k < landscapes.size() ? landscapes[k] : 0;
      }
    };
  };

#ifdef GUDHI_USE_TBB
  oneapi::tbb::task_arena arena(n_jobs);
  arena.execute([&] {
    tbb::parallel_for(std::size_t(0), grid[0].size(),
                      [&](std::size_t i) { tbb::parallel_for(std::size_t(0), grid[1].size(), get_image_values(i)); });
  });
#else
  for (std::size_t i = 0; i < grid[0].size(); ++i) {
    auto get_image_values_at = get_image_values(i);
    for (std::size_t j = 0; j < grid[1].size(); ++j) {
      get_image_values_at(j);
    }
  }
#endif

  return images;
}

/**
 * @ingroup multi_persistence
 *
 * @brief Computes the distance of all given points to all summands in the module.
 * TODO: proper definition of the distance.
 *
 * @tparam T Value type of a parameter in a filtration value.
 * @tparam RandomAccessPointRange Range with size() and operator[] method. The operator[] method must return a
 * type with the same methods and a value type convertible to `T`.
 * @param mod Module.
 * @param pts Range of points with number of coordinates corresponding to the number of parameters in the module.
 * @param negative If true, the distance is allowed to be signed.
 * @param n_jobs If TBB is linked, allows to specify the number of threads that should be used for parallelization.
 * @return A continuous vector of distance values which should be interpreted as a c-ordered 2-dimensional array with
 * a first axis corresponding to the points and a second axis to the summands.
 */
template <typename T, class RandomAccessPointRange>
inline std::vector<maybe_make_signed_t<T>> compute_module_distances_to(const Module<T> &mod,
                                                                       const RandomAccessPointRange &pts, bool negative,
                                                                       [[maybe_unused]] int n_jobs = 0) {
  std::vector<maybe_make_signed_t<T>> res(pts.size() * mod.size());
  Gudhi::Simple_mdspan data(res.data(), pts.size(), mod.size());

  auto get_distances_of_point = [&](std::size_t i) {
    for (std::size_t j = 0; j < data.extent(1); ++j) {
      data(i, j) = compute_summand_distance_to(mod.get_summand(j), pts[i], negative);
    }
  };

#ifdef GUDHI_USE_TBB
  oneapi::tbb::task_arena arena(n_jobs);  // limits the number of threads
  arena.execute([&] { tbb::parallel_for(std::size_t(0), pts.size(), get_distances_of_point); });
#else
  for (std::size_t i = 0; i < pts.size(); ++i) {
    get_distances_of_point(i);
  }
#endif

  return res;
}

/**
 * @ingroup multi_persistence
 *
 * @brief For a birth and death corner in a summand of the module, let the diagonal between those two be
 * \f$ min\{death[p] - birth[p]\} \f$ for all parameters \f$ p \f$. This method returns for all summands in the module
 * the maximal diagonal of all birth-death pairs in the intersection between the summand and the box.
 *
 * @tparam T Value type of a parameter in a filtration value.
 * @tparam U Template argument of @ref Box. Has to be either T or std::make_signed_t<T>.
 * @param mod Module.
 * @param box Box to intersect with. The box is ignored if trivial.
 */
template <typename T, typename U>
inline std::vector<maybe_make_signed_t<T>> compute_module_interleavings(const Module<T> &mod, const Box<U> &box) {
  static_assert(std::is_same_v<U, T> || std::is_same_v<U, maybe_make_signed_t<T>>,
                "Box template parameter is not compatible with Module value type.");

  std::vector<maybe_make_signed_t<T>> interleavings(mod.size());

  auto get_interleaving = [&](std::size_t i) {
    interleavings[i] = compute_summand_interleaving(mod.get_summand(i), box);
  };

#ifdef GUDHI_USE_TBB
  tbb::parallel_for(std::size_t(0), interleavings.size(), get_interleaving);
#else
  for (std::size_t i = 0; i < interleavings.size(); ++i) {
    get_interleaving(i);
  }
#endif

  return interleavings;
}

/**
 * @ingroup multi_persistence
 *
 * @private
 */
template <typename T, class RandomAccessValueRange>
inline double _get_module_pixel_value(typename Module<T>::const_iterator start, typename Module<T>::const_iterator end,
                                      const RandomAccessValueRange &x, double delta, double p, bool normalize,
                                      double moduleWeight, const std::vector<double> &interleavings) {
  double value = 0;

  if (p == 0) {
    for (auto it = start; it != end; it++) {
      value += compute_summand_local_weight(*it, x, delta);
    }
    if (normalize) value /= moduleWeight;
    return value;
  }

  if (p != Module<double>::T_inf) {
    std::size_t c = 0;
    for (auto it = start; it != end; it++) {
      double summandWeight = interleavings[c++];
      double summandXWeight = compute_summand_local_weight(*it, x, delta);
      value += std::pow(summandWeight, p) * summandXWeight;
    }
    if (normalize) value /= moduleWeight;
    return value;
  }

  for (auto it = start; it != end; it++) {
    value = std::max(value, compute_summand_local_weight(*it, x, delta));
  }
  return value;
}

/**
 * @ingroup multi_persistence
 *
 * @private
 */
template <typename T, typename U, class RandomAccessPointRange, class OutputIt>
inline void _compute_module_pixels_of_degree(typename Module<T>::const_iterator start,
                                             typename Module<T>::const_iterator end, double delta, double p,
                                             bool normalize, const Box<U> &box,
                                             const RandomAccessPointRange &coordinates, int n_jobs, OutputIt dFirst) {
  std::vector<double> interleavings(std::distance(start, end));

  auto compute_module_weight = [&](auto &&op) -> double {
    double moduleWeight = 0;
    std::size_t c = 0;
    for (auto it = start; it != end; it++) {
      interleavings[c] = compute_summand_interleaving<T, double>(*it, box);
      moduleWeight += op(moduleWeight, interleavings[c]);
      ++c;
    }
    return moduleWeight;
  };

  auto op = [p](double mw, double inter) -> double {
    if (p == 0) {
      return inter > 0;
    }
    if (p != Module<double>::T_inf) {
      // /!\ TODO: deal with inf summands (for the moment,  depends on the box...)
      if (inter > 0 && inter != Module<double>::T_inf) return std::pow(inter, p);
      return 0;
    }
    if (inter > 0 && inter != Module<double>::T_inf) return std::max(mw, inter);
    return 0;
  };

  double moduleWeight = compute_module_weight(op);
  std::cout << "moduleWeight: " << moduleWeight << "\n";
  if (moduleWeight == 0) return;

#ifdef GUDHI_USE_TBB
  oneapi::tbb::task_arena arena(n_jobs);  // limits the number of threads
  arena.execute([&] {
    tbb::parallel_for(std::size_t(0), coordinates.size(), [&](std::size_t i) {
      *(dFirst + i) =
          _get_module_pixel_value<T>(start, end, coordinates[i], delta, p, normalize, moduleWeight, interleavings);
    });
  });
#else
  for (std::size_t i = 0; i < coordinates.size(); ++i) {
    *(dFirst + i) =
        _get_module_pixel_value<T>(start, end, coordinates[i], delta, p, normalize, moduleWeight, interleavings);
  }
#endif
}

/**
 * @ingroup multi_persistence
 *
 * @brief Assumes that the summands in the module are ordered by increasing dimension. Computes the persistence images
 * of the module for the given dimensions.
 *
 * @tparam T Value type of a parameter in a filtration value.
 * @tparam RandomAccessPointRange Range with size() and operator[] method. The operator[] method must return a
 * type with the same methods and a value type convertible to `T`.
 * @tparam DimensionRange Integer range with a size() and begin() method.
 * @tparam U Template argument of @ref Box. Has to be either T or std::make_signed_t<T>. Default: T.
 * @param mod Module.
 * @param coordinates Image coordinates. One value will be computed for each of them for each requested dimension.
 * @param dimensions Range of dimensions to compute. Has to be ordered by increasing value.
 * @param box Box within to compute the module weight. Default: trivial box.
 * @param delta Radius. If positive, the weight computed is the interleaving distance to 0 of the
 * summand restricted to the current point and radius. If negative, the weight is the volume of the largest rectangle
 * spanned by a birth and a death corner of the summand intersected with the current point and radius.
 * @param p \f$ p \f$ coefficient of the \f$ p \f$-norm to use.
 * @param normalize Indicates if the values have to be normalized or not.
 * @param n_jobs If TBB is linked, allows to specify the number of threads that should be used for parallelization.
 * @return A continuous vector of pixel values which should be interpreted as a c-ordered 2-dimensional array with
 * a first axis corresponding to the dimensions and a second axis to the coordinates.
 */
template <typename T, class RandomAccessPointRange, class DimensionRange, typename U = T>
inline std::vector<double> compute_module_pixels(const Module<T> &mod, const RandomAccessPointRange &coordinates,
                                                 const DimensionRange &dimensions, const Box<U> &box = {},
                                                 double delta = 0.1, double p = 1, bool normalize = true,
                                                 int n_jobs = 0) {
  auto numDegrees = dimensions.size();
  std::vector<double> out(numDegrees * coordinates.size());
  Gudhi::Simple_mdspan view(out.data(), numDegrees, coordinates.size());

  auto start = mod.begin();
  auto end = mod.begin();
  auto dimIt = dimensions.begin();
  for (std::size_t degreeIdx = 0; degreeIdx < numDegrees; degreeIdx++) {
    auto d = *dimIt;
    ++dimIt;
    start = end;
    while (start != mod.end() && start->get_dimension() != d) start++;
    if (start == mod.end()) break;
    end = start;
    while (end != mod.end() && end->get_dimension() == d) end++;
    _compute_module_pixels_of_degree<T>(start, end, delta, p, normalize, box, coordinates, n_jobs, &view(degreeIdx, 0));
  }
  return out;
}

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // MP_MODULE_HELPERS_H_
