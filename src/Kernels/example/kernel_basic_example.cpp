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

#include <gudhi/kernel.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>


void usage(int nbArgs, char *const progName) {
  std::cerr << "Error: Number of arguments (" << nbArgs << ") is not correct\n";
  std::cerr << "Usage: " << progName << " PD1 PD2 \n";
  std::cerr << "       i.e.: " << progName << " ../../../../data/persistence_diagram/PD1 ../../../../data/persistence_diagram/PD2 \n";
  exit(-1);  // ----- >>
}

int main(int argc, char **argv) {

  if (argc != 3) usage(argc, argv[0]);

  double sigma = 2; double tau = 5;

  std::string PDname1(argv[1]); std::string PDname2(argv[2]);
  std::vector< std::pair<double, double> > v1, v2; std::string line; double b,d;

  std::ifstream input1(PDname1);
  while(std::getline(input1,line)){
    std::stringstream stream(line); stream >> b; stream >> d; v1.push_back(std::pair<double,double>(b,d));
  }

  std::ifstream input2(PDname2);
  while(std::getline(input2,line)){
    std::stringstream stream(line); stream >> b; stream >> d; v2.push_back(std::pair<double,double>(b,d));
  }

  std::cout << "SWK exact = "    << Gudhi::kernel::swk             (v1,v2,sigma)     << std::endl;
  std::cout << "SWK approx = "   << Gudhi::kernel::approx_swk      (v1,v2,sigma)     << std::endl;
  std::cout << "PSSK exact = "   << Gudhi::kernel::pssk            (v1,v2,sigma)     << std::endl;
  std::cout << "PSSK approx = "  << Gudhi::kernel::approx_pssk     (v1,v2,sigma)     << std::endl;
  std::cout << "LPWGK exact = "  << Gudhi::kernel::lpwgk           (v1,v2,sigma)     << std::endl;
  std::cout << "LPWGK approx = " << Gudhi::kernel::approx_lpwgk    (v1,v2,sigma)     << std::endl;
  std::cout << "GPWGK exact = "  << Gudhi::kernel::gpwgk           (v1,v2,sigma,tau) << std::endl;
  std::cout << "GPWGK approx = " << Gudhi::kernel::approx_gpwgk    (v1,v2,sigma,tau) << std::endl;

}
