/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Loiseaux
 *
 *    Copyright (C) 2021 Inria
 *
 *    Modification(s):
 *      - 2022/03 Hannah Schreiber: Integration of the new Vineyard_persistence class, renaming and cleanup.
 *      - 2022/05 Hannah Schreiber: Addition of Summand class and Module class.
 *      - 2026/06 Hannah Schreiber: Reorganizing and cleanup.
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file multiparameter_module_approximation.h
 * @author David Loiseaux
 * @brief Contains the @ref Gudhi::multi_persistence::multiparameter_module_approximation method.
 */

#ifndef MP_MULTIPARAM_MODULE_APPROX_H_
#define MP_MULTIPARAM_MODULE_APPROX_H_

#include <iostream>     // std::clog
#include <cstddef>      // std::size_t
#include <stdexcept>    // std::invalid_argument
#include <type_traits>  // std::is_signed_v
#include <utility>      // std::swap, std::move, std::forward
#include <cmath>        // std::ceil, std::fabs, std::max
#include <initializer_list>
#include <vector>
#include <array>

#ifdef GUDHI_USE_TBB
#include <oneapi/tbb/task_arena.h>
#include <oneapi/tbb/parallel_for.h>
#endif

#include <gudhi/Clock.h>
#include <gudhi/Debug_utils.h>
#include <gudhi/Multi_persistence/Line.h>
#include <gudhi/Multi_persistence/Box.h>
#include <gudhi/Multi_persistence/Point.h>
#include <gudhi/Multi_persistence/Module.h>
#include <gudhi/Slicer.h>

namespace Gudhi {
namespace multi_persistence {

/**
 * @ingroup multi_persistence
 *
 * @private
 * @brief Builds an approximation of the multi-parameter module from the given multi-parameter filtration.
 *
 * @tparam T Has to be signed.
 */
template <typename T>
class Multiparameter_module_approximator {
 public:
  template <class CoordinateRange = std::initializer_list<T>, class DirectionRange = std::initializer_list<T>>
  Multiparameter_module_approximator(const CoordinateRange& boxCorner1, const CoordinateRange& boxCorner2, T precision,
                                     bool thresholdToBox, const DirectionRange& direction = {}, bool verbose = false)
      : box_(),
        direction_(direction.begin(), direction.end()),
        gridSize_(boxCorner1.size()),
        signs_(boxCorner1.size()),
        precision_(precision),
        thresholdToBox_(thresholdToBox),
        verbose_(verbose) {
    static_assert(std::is_signed_v<T>, "Template parameter has to be a signed arithmetic type.");

    GUDHI_CHECK(boxCorner1.size() == boxCorner2.size(),
                std::invalid_argument("Box corners must have same numbers of coordinates."));
    GUDHI_CHECK(
        direction.size() == 0 || direction.size() == boxCorner1.size(),
        std::invalid_argument("The direction, if specified, and the box must have same numbers of coordinates."));

    const std::size_t numParameters = boxCorner1.size();
    int signsShifts = 0;
    std::size_t argMaxSignsShifts = -1;
    Coordinates bottomCorner(boxCorner1.begin(), boxCorner1.end());
    Coordinates topCorner(boxCorner2.begin(), boxCorner2.end());
    for (std::size_t i = 0; i < numParameters; i++) {
      auto& b = bottomCorner[i];
      auto& t = topCorner[i];
      gridSize_[i] = static_cast<int>(std::ceil((std::fabs(t - b) / precision_))) + 1;
      signs_[i] = (t > b);
      if (t < b) {
        std::swap(b, t);
        int localShift = gridSize_[i];
        if (direction_.size()) {
          localShift =
              direction_[i] > 0 ? static_cast<int>(std::ceil(static_cast<T>(gridSize_[i]) / direction_[i])) : 0;
        }
        if (localShift > signsShifts) {
          signsShifts = std::max(signsShifts, localShift);
          argMaxSignsShifts = i;
        }
      }
    }
    if (signsShifts > 0) {
      // this may be too much for large num_parameters
      for (std::size_t i = 0; i < numParameters; i++) gridSize_[i] += signsShifts;
      gridSize_[argMaxSignsShifts] = 1;
      _progress("Had to flatten/shift coordinate", argMaxSignsShifts, "by", signsShifts);
    }
    box_ = Box<T>(std::move(bottomCorner), std::move(topCorner));

    // TODO: check also that no corner is at infinity?
    GUDHI_CHECK(!box_.is_trivial(), std::invalid_argument("Given corners yield a trivial box."));

    _progress("Num parameters:", numParameters);
    _progress("Box:", box_);
    _progress("Grid size:", gridSize_);
    _progress("Signs:", signs_);
    _progress("Max error:", precision_);
    _progress("Threshold to box:", thresholdToBox_);
  }

  template <class MultiFiltrationValue, class PersistenceAlgorithm, class CoordinateRange = std::initializer_list<T>>
  Module<T> compute_approximation(Slicer<MultiFiltrationValue, PersistenceAlgorithm>& slicer,
                                  const CoordinateRange& basepoint, bool complete,
                                  [[maybe_unused]] int n_jobs = 0) const {
    static_assert(std::is_same_v<typename MultiFiltrationValue::value_type, T>,
                  "Slicer value type not compatible with class value type.");

    Module<T> out;
    const std::size_t numParameters = box_.get_number_of_coordinates();
    GUDHI_CHECK(
        basepoint.size() == numParameters,
        std::invalid_argument("Base point must have as many coordinates than number of parameters at initialization."));
    if (numParameters == 0) return out;

#ifdef GUDHI_USE_TBB
    // TODO: for now, only _approximate_module uses some parallelization, but it could be used a bit everywhere
    oneapi::tbb::task_arena arena(n_jobs);
    arena.execute([&]() {
      _initialize_module(out, slicer, basepoint);
      _approximate_module(out, slicer, basepoint);
      _clean_module(out, complete);
    });
#else
    _initialize_module(out, slicer, basepoint);
    _approximate_module(out, slicer, basepoint);
    _clean_module(out, complete);
#endif

    return out;
  }

 private:
  using Coordinates = typename Box<T>::Point_t;
  template <typename U>
  using Container = Point<U>;  // works like std::vector but printable
  using Point_t = typename Line<T>::Point_t;

  Box<T> box_;
  Coordinates direction_;
  Container<int> gridSize_;
  Container<bool> signs_;
  T precision_;
  bool thresholdToBox_;
  bool verbose_;

  template <bool sign>
  class LineIterator {
   public:
    using Point_t = typename Line<T>::Point_t;

    LineIterator(Point_t&& basepoint, const Point_t& direction, T precision, int num_iterations)
        : step_(precision), remainingIterations_(num_iterations), currentLine(std::move(basepoint), direction) {
      if constexpr (!sign) {
        step_ = -step_;  // T is always signed, so that is fine.
      }
    };

    const Line<T>& operator*() const { return currentLine; }

    LineIterator<sign>& next(std::size_t i) {
      if (is_finished()) return *this;
      auto& basepoint = currentLine.base_point();
      basepoint[i] += step_;
      --remainingIterations_;
      return *this;
    }

    [[nodiscard]] bool is_finished() const { return remainingIterations_ <= 0; }

   private:
    T step_;
    int remainingIterations_;
    Line<T> currentLine;
  };

  template <typename... Args>
  void _progress(const Args&... args) const {
    if (verbose_) {
      ((std::clog << args << " "), ...);
      std::clog << "\n";
    }
  }

  template <class S, class CoordinateRange>
  void _initialize_module(Module<T>& mod, S& slicer, const CoordinateRange& basepoint) const {
    Line<T> currentLine(basepoint.begin(), basepoint.end(), direction_.begin(), direction_.end());

    _progress("First line basepoint", currentLine.base_point());

    _progress("Initializing mma...");
    Gudhi::Clock timer("Initializing mma done");

    // fills the first barcode
    slicer.push_to(currentLine);
    slicer.initialize_persistence_computation(false);
    auto barcode = slicer.template get_flat_barcode<true>();
    auto numBars = 0;
    for (const auto& b : barcode) numBars += b.size();
    mod.resize(numBars, box_.get_number_of_coordinates());
    mod.set_max_dimension(barcode.size() - 1);
    std::size_t i = 0;
    for (unsigned int dim = 0; dim < barcode.size(); ++dim) {
      for ([[maybe_unused]] const auto& bar : barcode[dim]) {
        mod.get_summand(i).set_dimension(dim);
        ++i;
      }
    }
    _add_barcode_to_module(mod, currentLine, barcode);

    timer.end();
    _progress("Instantiated", numBars, "summands.");
    if (verbose_) timer.print();
  }

  template <class S, class CoordinateRange>
  void _approximate_module(Module<T>& mod, S& slicer, const CoordinateRange& basepoint) const {
    const std::size_t numParameters = box_.get_number_of_coordinates();
    _progress("Computing mma...");
    Gudhi::Clock timer("Computing mma done");

    for (std::size_t i = 1; i < numParameters; i++) {
      // the loop is on the faces of the lower box
      // should be parallelizable, up to a mutex on mod

      // skip faces with codim d_i=0
      if (!direction_.size() || direction_[i] > 0) {
        _progress("Face", i, " /", numParameters, "with grid size", gridSize_);
        _rec_approximate_module<0>(mod, Point_t(basepoint.begin(), basepoint.end()), i, numParameters - 1,
                                   slicer.weak_copy());
      }
    }

    // last one
    if (!direction_.size() || direction_[0] > 0) {
      _progress("Face", numParameters, " /", numParameters, "with grid size", gridSize_);
      _rec_approximate_module<1>(mod, Point_t(basepoint.begin(), basepoint.end()), 0, numParameters - 1,
                                 std::move(slicer));
    }

    timer.end();
    if (verbose_) timer.print();
  }

  void _clean_module(Module<T>& mod, bool complete) const {
    _progress("Cleaning output...");
    Gudhi::Clock timer("Cleaning output done");

    mod.clean();
    if (complete) {
      _progress("Completing output ...");
      for (std::size_t i = 0; i < box_.get_number_of_coordinates(); i++) mod.fill(precision_);
    }

    timer.end();
    if (verbose_) timer.print();
  }

  template <int axis, class Slicer>
  void _rec_approximate_module(Module<T>& module, Point_t&& basepoint, std::size_t gridSizeIdx, int paramToIterate,
                               Slicer&& currentPersistence) const {
    if (paramToIterate > axis && (static_cast<int>(gridSizeIdx) == paramToIterate || gridSize_[paramToIterate] == 0)) {
      // no need to copy basepoint, we just skip the dim here
      _rec_approximate_module<axis>(module, std::move(basepoint), gridSizeIdx, paramToIterate - 1,
                                    std::forward<Slicer>(currentPersistence));
      return;
    }

    if (signs_[paramToIterate]) {
      _rec_iterate_over_lines<axis, true>(module, std::move(basepoint), gridSizeIdx, paramToIterate,
                                          currentPersistence);
    } else {
      _rec_iterate_over_lines<axis, false>(module, std::move(basepoint), gridSizeIdx, paramToIterate,
                                           currentPersistence);
    }
  }

  template <int axis, bool sign, class Slicer>
  void _rec_iterate_over_lines(Module<T>& module, Point_t&& basepoint, std::size_t gridSizeIdx, int paramToIterate,
                               Slicer& currentPersistence) const {
    LineIterator<sign> line_iterator(std::move(basepoint), direction_, precision_,
                                     static_cast<int>(gridSizeIdx) == paramToIterate ? 0 : gridSize_[paramToIterate]);

    if (paramToIterate <= axis) {
      _add_vineyard_trajectory_to_module<axis, sign>(module, currentPersistence, line_iterator);
      return;
    }

    while (!line_iterator.is_finished()) {
      // TODO : multithread, but needs matrix to be thread safe + put mutex on module
      const auto& currPoint = (*line_iterator).base_point();
      _rec_approximate_module<axis>(module, Point_t(currPoint), gridSizeIdx, paramToIterate - 1,
                                    currentPersistence.weak_copy());
      line_iterator.next(paramToIterate);
    }
  }

  template <int axis, bool sign, class Slicer>
  void _add_vineyard_trajectory_to_module(Module<T>& module, Slicer& slicer, LineIterator<sign>& line_iterator) const {
    // Line iterator should be on the biggest axis
    while (!line_iterator.is_finished()) {
      const Line<T>& new_line = *(line_iterator.next(axis));
      slicer.push_to(new_line);
      slicer.update_persistence_computation();
      _add_barcode_to_module(module, new_line, slicer.template get_flat_barcode<true>());
    };
  }

  void _add_barcode_to_module(Module<T>& mod, const Line<T>& line,
                              const std::vector<std::vector<std::array<T, 2>>>& barcode) const {
#ifdef GUDHI_USE_TBB
    std::vector<std::size_t> shifts(barcode.size(), 0U);
    for (std::size_t i = 1U; i < barcode.size(); i++) {
      shifts[i] = shifts[i - 1] + barcode[i - 1].size();
    }
    tbb::parallel_for(size_t(0), barcode.size(), [&](size_t dim) {
      tbb::parallel_for(size_t(0), barcode[dim].size(), [&](size_t j) {
        mod.get_summand(shifts[dim] + j).add_bar(line, barcode[dim][j][0], barcode[dim][j][1], box_, thresholdToBox_);
      });
    });
#else
    typename Module<T>::Index count = 0;
    for (const auto& barDim : barcode) {
      for (const auto& bar : barDim) mod.get_summand(count++).add_bar(line, bar[0], bar[1], box_, thresholdToBox_);
    }
#endif
  }
};

/**
 * @ingroup multi_persistence
 *
 * @brief Computes the approximated module from the multi-parameter filtration stored in the given @ref Slicer
 * within the given box.
 *
 * The algorithm is based on @cite mma.
 *
 * @tparam MultiFiltrationValue First template parameter of @ref Slicer.
 * @tparam PersistenceAlgorithm Second template parameter of @ref Slicer.
 * @tparam CoordinateRange Range of an arithmetic value convertible to @ref MultiFiltrationValue::value_type.
 * Has to have a size(), begin() and end() method.
 * @tparam DirectionRange Range of an arithmetic value convertible to @ref MultiFiltrationValue::value_type.
 * Has to have a size(), begin() and end() method.
 * @param slicer @ref Slicer storing the multi-parameter filtration.
 * @param precision Approximation error.
 * @param boxCorner1 First corner of the box in which to compute the module. Has to have as many coordinates
 * than parameters in @p slicer.
 * @param boxCorner2 Second corner of the box, diagonally opposite to @p boxCorner1 such that the two corners
 * span the wanted box. Has to have as many coordinates than parameters in @p slicer.
 * @param direction Slope of the lines slicing the module to approximate. If empty, slope 1 is chosen.  Has to have
 * as many coordinates than parameters in @p slicer. Default: empty.
 * @param thresholdToBox If true, slice values outside the box are mapped to the border of the box. If false,
 * they are mapped to minus infinity (if below the box) or plus infinity (if above the box). Default: false.
 * @param complete If true, corners in the summands are identified if they are closer then @p precision. Default: true.
 * @param verbose If true, outputs information and timings in `std::clog`. Default: false.
 * @param n_jobs If TBB is linked, allows to specify the number of threads that should be used for parallelization.
 * Default: -1.
 * @return @ref Module with template parameter @ref MultiFiltrationValue::value_type.
 */
template <class MultiFiltrationValue, class PersistenceAlgorithm,
          class CoordinateRange = std::initializer_list<typename MultiFiltrationValue::value_type>,
          class DirectionRange = std::initializer_list<typename MultiFiltrationValue::value_type>>
Module<typename MultiFiltrationValue::value_type> multiparameter_module_approximation(
    Slicer<MultiFiltrationValue, PersistenceAlgorithm>& slicer, typename MultiFiltrationValue::value_type precision,
    const CoordinateRange& boxCorner1, const CoordinateRange& boxCorner2, const DirectionRange& direction = {},
    bool thresholdToBox = false, bool complete = true, bool verbose = false, [[maybe_unused]] int n_jobs = -1) {
  GUDHI_CHECK(slicer.get_number_of_parameters() == boxCorner1.size(),
              std::invalid_argument("Box corners need to have as many coordinates than parameters in the slicer."));

  Multiparameter_module_approximator mma(boxCorner1, boxCorner2, precision, thresholdToBox, direction, verbose);
  return mma.compute_approximation(slicer, boxCorner1, complete, n_jobs);
}

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // MP_MULTIPARAM_MODULE_APPROX_H_
