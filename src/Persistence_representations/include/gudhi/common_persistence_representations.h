/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - 2025/06 Hannah Schreiber: Various small bug fixes (missing `inline`s, `DEBUG_TRACES`s etc.)
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file common_persistence_representations.h
 * @author Pawel Dlotko
 * @brief This file contain an implementation of some common procedures used in the Persistence_representations module.
 */

#ifndef COMMON_PERSISTENCE_REPRESENTATIONS_H_
#define COMMON_PERSISTENCE_REPRESENTATIONS_H_

#include <stdexcept>  // std::invalid_argument
#include <utility>    // std::pair, std::make_pair
#include <sstream>    // std::ostringstream
#include <cmath>      // std::fabs
#include <vector>

#include <boost/math/constants/constants.hpp>

namespace Gudhi {
namespace Persistence_representations {

/**
 * In this module, we use the name Persistence_diagram for the representation of a diagram in a vector of pairs of two
 * double.
 *
 * @ingroup Persistence_representations
 */
using Persistence_diagram = std::vector<std::pair<double, double> >;

/**
 * Pi value.
 *
 * @ingroup Persistence_representations
 */
inline constexpr double pi = boost::math::constants::pi<double>();

// TODO: I don't think this is a good idea? Having modifiable global variables like that.
// But the description of `almost_equal` shows that it is the intended use.
// At the same time, it seems to be the only method using it. Can we not add it as an argument instead?
/**
 * Epsilon value for rounding tolerances. Initialized to 0.000005.
 * The value can be modified if another tolerance threshold is needed.
 *
 * @ingroup Persistence_representations
 */
inline double epsi = 0.000005;

/**
 * A procedure used to compare doubles. Typically given two doubles \f$ A \f$ and \f$ B \f$, comparing both with
 * \f$ A == B \f$ is not good idea. In this case, we use the procedure @ref almost_equal with the global variable
 * @ref epsi "". Setting up @ref epsi gives the user a tolerance on what should be consider equal.
 *
 * @ingroup Persistence_representations
 **/
inline bool almost_equal(double a, double b)
{
  if (std::fabs(a - b) < epsi) return true;
  return false;
}

// landscapes
/**
 * @private Extra functions needed in construction of barcodes.
 *
 * @ingroup Persistence_representations
 **/
inline constexpr double minus_length(const std::pair<double, double>& a) { return a.first - a.second; }

/**
 * @private Extra functions needed in construction of barcodes.
 *
 * @ingroup Persistence_representations
 */
inline constexpr double birth_plus_deaths(const std::pair<double, double>& a) { return a.first + a.second; }

// landscapes
/**
 * @private
 * This method provides a comparison of points that is used in construction of persistence landscapes. The ordering is
 * lexicographical for the first coordinate, and reverse-lexicographical for the second coordinate.
 *
 * @ingroup Persistence_representations
 **/
inline constexpr bool compare_points_sorting(const std::pair<double, double>& f, const std::pair<double, double>& s)
{
  if (f.first < s.first) {
    return true;
  } else {  // f.first >= s.first
    if (f.first > s.first) {
      return false;
    } else {  // f.first == s.first
      if (f.second > s.second) {
        return true;
      } else {
        return false;
      }
    }
  }
}

// TODO: the following methods are mostly used with `std::make_pair` as argument in `Persistence_landscape(_on_grid)".
// If called often, it would be more effective to have a (double, double, double, double) version instead.

// landscapes
/**
 * Given two points in \f$ R^2 \f$, the procedure compute the parameters \f$ A \f$ and \f$ B \f$ of the line
 * \f$ y = Ax + B \f$ that goes through those two points.
 *
 * @ingroup Persistence_representations
 **/
inline std::pair<double, double> compute_parameters_of_a_line(const std::pair<double, double>& p1,
                                                              const std::pair<double, double>& p2)
{
  const double a = (p2.second - p1.second) / (p2.first - p1.first);
  const double b = p1.second - a * p1.first;
  return std::make_pair(a, b);
}

// landscapes
/**
 * This procedure, given two points which lies on the opposite sides of the \f$ x \f$-axis, computes the value of
 * \f$ x \f$ for which the line connecting those two points crosses the \f$ x \f$-axis.
 **/
inline double find_zero_of_a_line_segment_between_those_two_points(const std::pair<double, double>& p1,
                                                                   const std::pair<double, double>& p2)
{
  if (p1.first == p2.first) return p1.first;
  if (p1.second * p2.second > 0) {
    std::ostringstream errMessage;
    errMessage << "In function find_zero_of_a_line_segment_between_those_two_points the arguments are: (" << p1.first
               << "," << p1.second << ") and (" << p2.first << "," << p2.second
               << "). There is no zero in line between those two points. Program terminated.";
    throw std::invalid_argument(errMessage.str());
  }
  // we assume here, that x \in [ p1.first, p2.first ] and p1 and p2 are points between which we will put the line
  // segment
  const auto [a, b] = compute_parameters_of_a_line(p1, p2);
  return -b / a;
}

// landscapes
/**
 * This procedure takes two points in \f$ R^2 \f$ and a double value \f$ x \f$. It computes the line parsing through
 * those two points and return the value of that linear function at \f$ x \f$.
 *
 * @ingroup Persistence_representations
 **/
inline double function_value(const std::pair<double, double>& p1, const std::pair<double, double>& p2, double x)
{
  // we assume here, that x \in [ p1.first, p2.first ] and p1 and p2 are points between which we will put the line
  // segment
  const auto [a, b] = compute_parameters_of_a_line(p1, p2);
  return (a * x + b);
}

}  // namespace Persistence_representations
}  // namespace Gudhi

#endif  // COMMON_PERSISTENCE_REPRESENTATIONS_H_
