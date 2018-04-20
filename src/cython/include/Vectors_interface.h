/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Mathieu Carriere
 *
 *    Copyright (C) 2018 INRIA
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

#ifndef INCLUDE_VECTORS_INTERFACE_H_
#define INCLUDE_VECTORS_INTERFACE_H_

#include <gudhi/Landscape.h>
#include <gudhi/Persistence_image.h>
#include <gudhi/Persistence_weighted_gaussian.h>

#include <iostream>
#include <vector>
#include <utility>  // for std::pair

using Weight = std::function<double (std::pair<double,double>) >;

namespace Gudhi {

namespace persistence_diagram {

  std::vector<std::vector<double> > compute_ls(const std::vector<std::pair<double, double> >& diag, int nb_ls, double min_x, double max_x, int res_x) {
    Gudhi::Persistence_representations::Landscape L(diag, nb_ls, min_x, max_x, res_x);
    return L.vectorize();
  }

  std::vector<std::vector<double> > compute_pim(const std::vector<std::pair<double, double> >& diag, double min_x, double max_x, int res_x, double min_y, double max_y, int res_y, std::string weight, double sigma, double C, double p) {
    Weight weight_fn;
    if(weight.compare("linear") == 0)  weight_fn = Gudhi::Persistence_representations::Persistence_weighted_gaussian::linear_weight;
    if(weight.compare("arctan") == 0)  weight_fn = Gudhi::Persistence_representations::Persistence_weighted_gaussian::arctan_weight(C,p);
    if(weight.compare("const")  == 0)  weight_fn = Gudhi::Persistence_representations::Persistence_weighted_gaussian::const_weight;
    Gudhi::Persistence_representations::Persistence_image P(diag, min_x, max_x, res_x, min_y, max_y, res_y, weight_fn, sigma);
    return P.vectorize();
  }

}  // namespace persistence_diagram

}  // namespace Gudhi


#endif  // INCLUDE_VECTORS_INTERFACE_H_
