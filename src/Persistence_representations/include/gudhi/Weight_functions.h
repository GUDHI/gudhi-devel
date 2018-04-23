/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Mathieu Carriere
 *
 *    Copyright (C) 2018  INRIA (France)
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

#ifndef WEIGHT_FUNCTIONS_H_
#define WEIGHT_FUNCTIONS_H_

// gudhi include
#include <gudhi/read_persistence_from_file.h>
#include <gudhi/common_persistence_representations.h>

// standard include
#include <cmath>
#include <iostream>
#include <vector>
#include <limits>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include <utility>
#include <functional>

namespace Gudhi {
namespace Persistence_representations {

/** \fn  static double pss_weight(std::pair<double,double> p)
 * \brief Persistence Scale Space kernel weight function.
 * \ingroup Persistence_representations
 *
 * @param[in] p point in 2D.
 */
static double pss_weight(std::pair<double,double> p)     {if(p.second > p.first)  return 1; else return -1;}

/** \fn static double linear_weight(std::pair<double,double> p)
 * \brief Linear weight function.
 * \ingroup Persistence_representations
 *
 * @param[in] p point in 2D.
 */
static double linear_weight(std::pair<double,double> p)  {return std::abs(p.second - p.first);}

/** \fn static double const_weight(std::pair<double,double> p)
 * \brief Constant weight function.
 * \ingroup Persistence_representations
 *
 * @param[in] p point in 2D.
 */
static double const_weight(std::pair<double,double> p)   {return 1;}

/** \fn static std::function<double (std::pair<double,double>) > arctan_weight(double C, double alpha)
 * \brief Returns the arctan weight function with parameters C and alpha.
 * \ingroup Persistence_representations
 *
 * @param[in] C     positive constant.
 * @param[in] alpha positive power.
 */
static std::function<double (std::pair<double,double>) > arctan_weight(double C, double alpha)  {return [=](std::pair<double,double> p){return C * atan(std::pow(std::abs(p.second - p.first), alpha));};}

}  // namespace Persistence_representations
}  // namespace Gudhi

#endif  // WEIGHT_FUNCTIONS_H_
