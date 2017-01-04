/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Authors:       Francois Godi, small modifications by Pawel Dlotko
 *
 *    Copyright (C) 2015  INRIA
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

#include <gudhi/Bottleneck.h>
#include <iostream>

int main() {

  std::vector< std::pair<double, double> > v1, v2;

  v1.emplace_back(2.7, 3.7);
  v1.emplace_back(9.6, 14.);
  v1.emplace_back(34.2, 34.974);
  v1.emplace_back(3., std::numeric_limits<double>::infinity());

  v2.emplace_back(2.8, 4.45);
  v2.emplace_back(9.5, 14.1);
  v2.emplace_back(3.2, std::numeric_limits<double>::infinity());


  double b = Gudhi::persistence_diagram::bottleneck_distance(v1, v2);

  std::cout << "Bottleneck distance = " << b << std::endl;

  b = Gudhi::persistence_diagram::bottleneck_distance(v1, v2, 0.1);

  std::cout << "Approx bottleneck distance = " << b << std::endl;

}
