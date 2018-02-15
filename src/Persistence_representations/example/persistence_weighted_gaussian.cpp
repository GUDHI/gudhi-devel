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

#include <gudhi/Persistence_weighted_gaussian.h>

#include <iostream>
#include <vector>
#include <utility>

using PD = std::vector<std::pair<double,double> >;
using PWG = Gudhi::Persistence_representations::Persistence_weighted_gaussian;

int main(int argc, char** argv) {

  std::vector<std::pair<double, double> > persistence1;
  std::vector<std::pair<double, double> > persistence2;

  persistence1.push_back(std::make_pair(1, 2));
  persistence1.push_back(std::make_pair(6, 8));
  persistence1.push_back(std::make_pair(0, 4));
  persistence1.push_back(std::make_pair(3, 8));

  persistence2.push_back(std::make_pair(2, 9));
  persistence2.push_back(std::make_pair(1, 6));
  persistence2.push_back(std::make_pair(3, 5));
  persistence2.push_back(std::make_pair(6, 10));

  PWG PWG1(persistence1);
  PWG PWG2(persistence2);
  double sigma = 1;
  double tau = 1;
  int m = 1000;



  // Linear PWG

  std::cout << PWG1.compute_scalar_product  (PWG2, sigma, PWG::arctan_weight, m) << std::endl;
  std::cout << PWG1.compute_scalar_product  (PWG2, sigma, PWG::arctan_weight, -1) << std::endl;

  std::cout << PWG1.distance  (PWG2, sigma, PWG::arctan_weight, m) << std::endl;
  std::cout << PWG1.distance  (PWG2, sigma, PWG::arctan_weight, -1) << std::endl;







  // Gaussian PWG

  std::cout << std::exp( -PWG1.distance (PWG2, sigma, PWG::arctan_weight, m, 2) ) / (2*tau*tau) << std::endl;
  std::cout << std::exp( -PWG1.distance (PWG2, sigma, PWG::arctan_weight, -1, 2) ) / (2*tau*tau) << std::endl;







  // PSS

  PD pd1 = persistence1; int numpts = persistence1.size();    for(int i = 0; i < numpts; i++)  pd1.emplace_back(persistence1[i].second,persistence1[i].first);
  PD pd2 = persistence2;     numpts = persistence2.size();    for(int i = 0; i < numpts; i++)  pd2.emplace_back(persistence2[i].second,persistence2[i].first);

  PWG pwg1(pd1);
  PWG pwg2(pd2);

  std::cout << pwg1.compute_scalar_product  (pwg2, 2*std::sqrt(sigma), PWG::pss_weight, m) / (16*pi*sigma) << std::endl;
  std::cout << pwg1.compute_scalar_product  (pwg2, 2*std::sqrt(sigma), PWG::pss_weight, -1) / (16*pi*sigma) << std::endl;

  std::cout << pwg1.distance  (pwg2, 2*std::sqrt(sigma), PWG::pss_weight, m) / (16*pi*sigma) << std::endl;
  std::cout << pwg1.distance  (pwg2, 2*std::sqrt(sigma), PWG::pss_weight, -1) / (16*pi*sigma) << std::endl;


  return 0;
}
