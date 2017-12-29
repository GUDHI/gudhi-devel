/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Authors:       Mathieu Carrière
 *
 *    Copyright (C) 2017  INRIA
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

#define NUMPI 3.14159265359
#include <gudhi/kernel.h>

int main() {

  std::vector< std::pair<double, double> > v1, v2;

  double sigma = 2; double tau = 5;

  v1.emplace_back(2.7, 3.7);
  v1.emplace_back(9.6, 14.);
  v1.emplace_back(34.2, 34.974);

  v2.emplace_back(2.8, 4.45);
  v2.emplace_back(9.5, 14.1);

  std::cout << "SW exact = "    << Gudhi::kernel::sw             (v1, v2)          << std::endl;
  std::cout << "SW approx = "   << Gudhi::kernel::approx_sw      (v1, v2)          << std::endl;
  std::cout << "PSS exact = "   << Gudhi::kernel::pss            (v1,v2,sigma)     << std::endl;
  std::cout << "PSS approx = "  << Gudhi::kernel::approx_pss     (v1,v2,sigma)     << std::endl;
  std::cout << "PWG exact = "   << Gudhi::kernel::lpwg           (v1,v2,sigma)     << std::endl;
  std::cout << "PWG approx = "  << Gudhi::kernel::approx_lpwg    (v1,v2,sigma)     << std::endl;
  std::cout << "GPWG exact = "  << Gudhi::kernel::gpwg           (v1,v2,sigma,tau) << std::endl;
  std::cout << "GPWG approx = " << Gudhi::kernel::approx_gpwg    (v1,v2,sigma,tau) << std::endl;

}
