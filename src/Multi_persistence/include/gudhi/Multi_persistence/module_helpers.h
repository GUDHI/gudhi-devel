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

// RandomAccessValueRange: std::vector<int>
template <typename T, class RandomAccessValueRange>
Module<T> build_permuted_module(const Module<T> &module, const RandomAccessValueRange &permutation) {
  GUDHI_CHECK(permutation.size() <= module.size(),
              std::invalid_argument("Permutation size is greater than the module size."));

  Module<T> out(module.get_box());

  if (module.size() == 0) return out;

  out.resize(module.size(), 1);
  for (std::size_t i = 0; i < permutation.size(); ++i) {
    out.get_summand(i) = module.get_summand(permutation[i]);
  }
  out.set_max_dimension(module.get_max_dimension());
  return out;
}

/**
 * @private
 */
template <typename T>
inline std::vector<T> _get_module_landscape_values(const Module<T> &mod, const std::vector<T> &x,
                                                   typename Module<T>::Dimension dimension) {
  std::vector<T> values;
  values.reserve(mod.size());
  for (std::size_t i = 0; i < mod.size(); i++) {
    const auto &summand = mod.get_summand(i);
    if (summand.get_dimension() == dimension) values.push_back(compute_summand_landscape_value(summand, x));
  }
  std::sort(values.begin(), values.end(), [](T x, T y) { return x > y; });
  return values;
}

// // TODO: extend in higher resolution dimension
// template <typename T, class RandomAccessValueRange>
// inline std::vector<T> compute_module_landscape(const Module<T> &mod, typename Module<T>::Dimension dimension,
//                                                unsigned int k, const Box<T> &box,
//                                                const RandomAccessValueRange &resolution) {
//   GUDHI_CHECK(resolution.size() >= 2, std::invalid_argument("Not enough resolution values."));

//   std::vector<T> image;
//   image.resize(resolution[0] * resolution[1]);
//   Simple_mdspan view(image.data(), resolution[0], resolution[1]);
//   T stepX = (box.get_upper_corner()[0] - box.get_lower_corner()[0]) / resolution[0];
//   T stepY = (box.get_upper_corner()[1] - box.get_lower_corner()[1]) / resolution[1];

//   auto get_image_values = [&](unsigned int i) {
//     return [&, i](unsigned int j) {
//       auto landscape = _get_module_landscape_values(
//           mod, {box.get_lower_corner()[0] + (stepX * i), box.get_lower_corner()[1] + (stepY * j)}, dimension);
//       view(i, j) = k < landscape.size() ? landscape[k] : 0;
//     };
//   };

// #ifdef GUDHI_USE_TBB
//   tbb::parallel_for(0U, resolution[0],
//                     [&](unsigned int i) { tbb::parallel_for(0U, resolution[1], get_image_values(i)); });
// #else
//   for (unsigned int i = 0; i < resolution[0]; ++i) {
//     auto get_image_values_at = get_image_values(i);
//     for (unsigned int j = 0; j < resolution[1]; ++j) {
//       get_image_values_at(j);
//     }
//   }
// #endif

//   return image;
// }

template <typename T, class RandomAccessValueRange1, class RandomAccessValueRange2>
inline std::vector<T> compute_set_of_module_landscapes(const Module<T> &mod, typename Module<T>::Dimension dimension,
                                                       const RandomAccessValueRange1 &ks, const Box<T> &box,
                                                       const RandomAccessValueRange2 &resolution,
                                                       [[maybe_unused]] int n_jobs = 0) {
  GUDHI_CHECK(resolution.size() >= 2, std::invalid_argument("Not enough resolution values."));

  using ResT = typename RandomAccessValueRange2::value_type;

  std::vector<T> images(ks.size() * resolution[0] * resolution[1]);
  Simple_mdspan view(images.data(), ks.size(), resolution[0], resolution[1]);
  T stepX = (box.get_upper_corner()[0] - box.get_lower_corner()[0]) / resolution[0];
  T stepY = (box.get_upper_corner()[1] - box.get_lower_corner()[1]) / resolution[1];

  auto get_image_values = [&](unsigned int i) {
    return [&, i](unsigned int j) {
      std::vector<T> landscapes = _get_module_landscape_values(
          mod, {box.get_lower_corner()[0] + (stepX * i), box.get_lower_corner()[1] + (stepY * j)}, dimension);
      for (std::size_t k_idx = 0; k_idx < ks.size(); ++k_idx) {
        unsigned int k = ks[k_idx];
        view(k_idx, i, j) = k < landscapes.size() ? landscapes[k] : 0;
      }
    };
  };

#ifdef GUDHI_USE_TBB
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

template <typename T, class RandomAccessValueRange, class RandomAccessArray>
inline std::vector<T> compute_set_of_module_landscapes(const Module<T> &mod, typename Module<T>::Dimension dimension,
                                                       const RandomAccessValueRange &ks,
                                                       const std::vector<RandomAccessArray> &grid,
                                                       [[maybe_unused]] int n_jobs = 0) {
  GUDHI_CHECK(grid.size() >= 2, std::invalid_argument("First axis of the grid has not enough values."));

  std::vector<T> images(ks.size() * grid[0].size() * grid[1].size());
  Simple_mdspan view(images.data(), ks.size(), grid[0].size(), grid[1].size());

  auto get_image_values = [&](std::size_t i) {
    return [&, i](std::size_t j) {
      std::vector<T> landscapes = _get_module_landscape_values(mod, {grid[0][i], grid[1][j]}, dimension);
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

// template <typename T, class MultiFiltrationValue>
// inline std::vector<int> compute_module_euler_curve(const Module<T> &mod,
//                                                    const std::vector<MultiFiltrationValue> &points) {
//   std::vector<int> eulerCurve(points.size());

//   auto get_curve = [&](std::size_t i) {
//     for (const auto &sum : mod) {
//       if (sum.contains(points[i])) {
//         int sign = sum.get_dimension() % 2 ? -1 : 1;
//         eulerCurve[i] += sign;
//       }
//     }
//   };

// #ifdef GUDHI_USE_TBB
//   tbb::parallel_for(std::size_t(0), eulerCurve.size(), get_curve);
// #else
//   for (std::size_t i = 0; i < eulerCurve.size(); ++i) {
//     get_curve(i);
//   }
// #endif

//   return eulerCurve;
// }

template <typename T, class RandomAccessPointRange>
inline std::vector<T> compute_module_distances_to(const Module<T> &mod, const RandomAccessPointRange &pts,
                                                  bool negative, int n_jobs) {
  std::vector<T> res(pts.size() * mod.size());
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

// template <typename T, class RandomAccessPointRange>
// inline std::vector<typename Module<T>::Summand_t::Index> compute_module_lower_and_upper_generators_of(
//     const Module<T> &mod, const RandomAccessPointRange &pts, bool full, int n_jobs) {
//   unsigned int tupleSize = full ? 4 : 2;
//   std::vector<typename Module<T>::Summand_t::Index> out(pts.size() * mod.size() * tupleSize);
//   Gudhi::Simple_mdspan data(out.data(), pts.size(), mod.size(), tupleSize);

//   auto get_generators = [&](std::size_t i) {
//     return [&, i](std::size_t j) {
//       compute_summand_lower_and_upper_generator_of(mod.get_summand(j), pts[i], full, &data(i, j, 0));
//     };
//   };

// #ifdef GUDHI_USE_TBB
//   oneapi::tbb::task_arena arena(n_jobs);  // limits the number of threads
//   arena.execute([&] {
//     tbb::parallel_for(std::size_t(0), pts.size(),
//                       [&](std::size_t i) { tbb::parallel_for(std::size_t(0), mod.size(), get_generators(i)); });
//   });
// #else
//   for (std::size_t i = 0; i < pts.size(); ++i) {
//     auto get_generators_at = get_generators(i);
//     for (std::size_t j = 0; j < mod.size(); ++j) {
//       get_generators_at(j);
//     }
//   }
// #endif

//   return out;
// }

template <typename T>
inline std::vector<T> compute_module_interleavings(const Module<T> &mod, const Box<T> &box) {
  std::vector<T> interleavings(mod.size());

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
 * @private
 */
template <typename T, class RandomAccessValueRange>
inline T _get_module_pixel_value(typename Module<T>::const_iterator start, typename Module<T>::const_iterator end,
                                 const RandomAccessValueRange &x, T delta, T p, bool normalize, T moduleWeight,
                                 const std::vector<T> &interleavings) {
  T value = 0;

  if (p == 0) {
    for (auto it = start; it != end; it++) {
      value += compute_summand_local_weight(*it, x, delta);
    }
    if (normalize) value /= moduleWeight;
    return value;
  }

  if (p != Module<T>::T_inf) {
    std::size_t c = 0;
    for (auto it = start; it != end; it++) {
      T summandWeight = interleavings[c++];
      T summandXWeight = compute_summand_local_weight(*it, x, delta);
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
 * @private
 */
template <typename T, class RandomAccessPointRange, class OutputIt>
inline void _compute_module_pixels_of_degree(typename Module<T>::const_iterator start,
                                             typename Module<T>::const_iterator end, T delta, T p, bool normalize,
                                             const Box<T> &box, const RandomAccessPointRange &coordinates, int n_jobs,
                                             OutputIt dFirst) {
  std::vector<T> interleavings(std::distance(start, end));

  auto compute_module_weight = [&](auto &&op) -> T {
    T moduleWeight = 0;
    std::size_t c = 0;
    for (auto it = start; it != end; it++) {
      interleavings[c] = compute_summand_interleaving(*it, box);
      moduleWeight += op(moduleWeight, interleavings[c]);
      ++c;
    }
    return moduleWeight;
  };

  auto op = [&](T mw, T inter) -> T {
    if (p == 0) {
      return inter > 0;
    }
    if (p != Module<T>::T_inf) {
      // /!\ TODO: deal with inf summands (for the moment,  depends on the box...)
      if (inter > 0 && inter != Module<T>::T_inf) return std::pow(inter, p);
      return 0;
    }
    if (inter > 0 && inter != Module<T>::T_inf) return std::max(mw, inter);
    return 0;
  };

  T moduleWeight = compute_module_weight(op);
  if (moduleWeight == 0) return;

#ifdef GUDHI_USE_TBB
  oneapi::tbb::task_arena arena(n_jobs);  // limits the number of threads
  arena.execute([&] {
    tbb::parallel_for(std::size_t(0), coordinates.size(), [&](std::size_t i) {
      *(dFirst + i) =
          _get_module_pixel_value(start, end, coordinates[i], delta, p, normalize, moduleWeight, interleavings);
    });
  });
#else
  for (std::size_t i = 0; i < coordinates.size(); ++i) {
    *(dFirst + i) =
        _get_module_pixel_value(start, end, coordinates[i], delta, p, normalize, moduleWeight, interleavings);
  }
#endif
}

template <typename T, class RandomAccessPointRange, class DimensionRange>
inline std::vector<T> compute_module_pixels(const Module<T> &mod, const RandomAccessPointRange &coordinates,
                                            const DimensionRange &dimensions, const Box<T> &box = {}, T delta = 0.1,
                                            T p = 1, bool normalize = true, int n_jobs = 0) {
  auto numDegrees = dimensions.size();
  std::vector<T> out(numDegrees * coordinates.size());
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
    _compute_module_pixels_of_degree(start, end, delta, p, normalize, box, coordinates, n_jobs, &view(degreeIdx, 0));
  }
  return out;
}

// /**
//  * @private
//  */
// template <typename T, class RandomAccessArray>
// inline std::vector<int> _project_generator_into_grid(const Gudhi::multi_filtration::Multi_parameter_generator<T> &pt,
//                                                      const std::vector<RandomAccessArray> &grid) {
//   std::size_t num_parameters = grid.size();
//   std::vector<int> out(num_parameters);
//   if (pt.is_plus_inf() || pt.is_nan()) [[unlikely]] {
//     for (std::size_t i = 0; i < num_parameters; ++i) out[i] = grid[i].size() - 1;
//     return out;
//   }
//   if (pt.is_minus_inf()) [[unlikely]] {
//     for (std::size_t i = 0; i < num_parameters; ++i) out[i] = 0;
//     return out;
//   }
//   // pt has to be of size num_parameters now
//   for (std::size_t i = 0; i < num_parameters; ++i) {
//     if (pt[i] >= grid[i].back()) [[unlikely]]
//       out[i] = grid[i].size() - 1;
//     else if (pt[i] <= grid[i][0]) [[unlikely]] {
//       out[i] = 0;
//     } else {
//       auto temp = std::distance(grid[i].begin(), std::lower_bound(grid[i].begin(), grid[i].end(), pt[i]));
//       if (std::abs(grid[i][temp] - pt[i]) < std::abs(grid[i][temp - 1] - pt[i])) {
//         out[i] = temp;
//       } else {
//         out[i] = temp - 1;
//       }
//     }
//   }
//   return out;
// }

// template <typename T, class RandomAccessArray>
// inline std::vector<std::vector<std::vector<int>>> project_module_into_grid(const Module<T> &mod,
//                                                                            const std::vector<RandomAccessArray> &grid) {
//   std::vector<std::vector<std::vector<int>>> out(3);
//   auto &idx = out[0];
//   auto &births = out[1];
//   auto &deaths = out[2];

//   idx.resize(2);
//   idx[0].resize(mod.size());
//   idx[1].resize(mod.size());

//   // some heuristic: usually
//   births.reserve(2 * mod.size());
//   deaths.reserve(2 * mod.size());
//   for (std::size_t i = 0; i < mod.size(); ++i) {
//     auto &interval = mod.get_summand(i);
//     idx[0][i] = interval.get_upset().size();
//     for (const auto &pt : interval.get_upset()) {
//       births.push_back(_project_generator_into_grid(pt, grid));
//     }
//     idx[1][i] = interval.get_downset().size();
//     for (const auto &pt : interval.get_downset()) {
//       deaths.push_back(_project_generator_into_grid(pt, grid));
//     }
//   }
//   return out;
// }

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // MP_MODULE_HELPERS_H_
