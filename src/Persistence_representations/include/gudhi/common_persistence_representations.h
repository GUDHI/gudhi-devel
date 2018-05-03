/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016 INRIA (France)
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef COMMON_PERSISTENCE_REPRESENTATIONS_H_
#define COMMON_PERSISTENCE_REPRESENTATIONS_H_

#include <utility>
#include <string>
#include <cmath>
#include <boost/math/constants/constants.hpp>



namespace Gudhi {
namespace Persistence_representations {
// this file contain an implementation of some common procedures used in Persistence_representations.

static constexpr double pi = boost::math::constants::pi<double>();


/**
 * In this module, we use the name Persistence_diagram for the representation of a diagram in a vector of pairs of two double.
 */
using Persistence_diagram = std::vector<std::pair<double,double> >;

/**
 * In this module, we use the name Weight for the representation of a function taking a pair of two double and returning a double.
 */
using Weight = std::function<double (std::pair<double,double>) >;
// double epsi = std::numeric_limits<double>::epsilon();
double epsi = 0.000005;

/**
 *  A procedure used to compare doubles. Typically given two doubles A and B, comparing A == B is not good idea. In this
 *  case, we use the procedure almostEqual with the epsi defined at
 *  the top of the file. Setting up the epsi gives the user a tolerance on what should be consider equal.
**/
inline bool almost_equal(double a, double b) {
  if (std::fabs(a - b) < epsi) return true;
  return false;
}

// landscapes
/**
 * Extra functions needed in construction of barcodes.
**/
double minus_length(std::pair<double, double> a) { return a.first - a.second; }
double birth_plus_deaths(std::pair<double, double> a) { return a.first + a.second; }

// landscapes
/**
 * Given two points in R^2, the procedure compute the parameters A and B of the line y = Ax + B that crosses those two points.
**/
std::pair<double, double> compute_parameters_of_a_line(std::pair<double, double> p1, std::pair<double, double> p2) {
  double a = (p2.second - p1.second) / (p2.first - p1.first);
  double b = p1.second - a * p1.first;
  return std::make_pair(a, b);
}

// landscapes
/**
 * This procedure given two points which lies on the opposite sides of x axis, compute x for which the line connecting those two points crosses x axis.
**/
double find_zero_of_a_line_segment_between_those_two_points(std::pair<double, double> p1,
                                                            std::pair<double, double> p2) {
  if (p1.first == p2.first) return p1.first;
  if (p1.second * p2.second > 0) {
    std::ostringstream errMessage;
    errMessage << "In function find_zero_of_a_line_segment_between_those_two_points the arguments are: (" << p1.first
               << "," << p1.second << ") and (" << p2.first << "," << p2.second
               << "). There is no zero in line between those two points. Program terminated.";
    std::string errMessageStr = errMessage.str();
    const char* err = errMessageStr.c_str();
    throw(err);
  }
  // we assume here, that x \in [ p1.first, p2.first ] and p1 and p2 are points between which we will put the line
  // segment
  double a = (p2.second - p1.second) / (p2.first - p1.first);
  double b = p1.second - a * p1.first;
  return -b / a;
}

// landscapes
/**
 * This method provides a comparison of points that is used in construction of persistence landscapes. The ordering is
 * lexicographical for the first coordinate, and reverse-lexicographical for the second coordinate.
**/
bool compare_points_sorting(std::pair<double, double> f, std::pair<double, double> s) {
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

// landscapes
/**
 * This procedure takes two points in R^2 and a double value x. It computes the line parsing through those two points
 *and return the value of that linear function at x.
**/
double function_value(std::pair<double, double> p1, std::pair<double, double> p2, double x) {
  // we assume here, that x \in [ p1.first, p2.first ] and p1 and p2 are points between which we will put the line
  // segment
  double a = (p2.second - p1.second) / (p2.first - p1.first);
  double b = p1.second - a * p1.first;
  return (a * x + b);
}

}  // namespace Persistence_representations
}  // namespace Gudhi

#endif  // COMMON_PERSISTENCE_REPRESENTATIONS_H_
