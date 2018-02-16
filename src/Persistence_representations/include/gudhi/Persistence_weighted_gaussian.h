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

#ifndef PERSISTENCE_WEIGHTED_GAUSSIAN_H_
#define PERSISTENCE_WEIGHTED_GAUSSIAN_H_

#ifdef GUDHI_USE_TBB
#include <tbb/parallel_for.h>
#endif

// gudhi include
#include <gudhi/read_persistence_from_file.h>

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
#include <boost/math/constants/constants.hpp>

double pi = boost::math::constants::pi<double>();
using PD = std::vector<std::pair<double,double> >;
using Weight = std::function<double (std::pair<double,double>) >;

namespace Gudhi {
namespace Persistence_representations {

class Persistence_weighted_gaussian{

 protected:
    PD diagram;
    Weight weight;
    double sigma;
    int approx;

 public:

  Persistence_weighted_gaussian(PD _diagram){diagram = _diagram; sigma = 1.0; approx = 1000; weight = arctan_weight;}
  Persistence_weighted_gaussian(PD _diagram, double _sigma, int _approx, Weight _weight){diagram = _diagram; sigma = _sigma; approx = _approx; weight = _weight;}
  PD get_diagram(){return this->diagram;}
  double get_sigma(){return this->sigma;}
  int get_approx(){return this->approx;}
  Weight get_weight(){return this->weight;}


  // **********************************
  // Utils.
  // **********************************


  static double pss_weight(std::pair<double,double> p){
    if(p.second > p.first)  return 1;
    else return -1;
  }

  static double arctan_weight(std::pair<double,double> p){
    return atan(p.second - p.first);
  }

  std::vector<std::pair<double,double> > Fourier_feat(PD diag, std::vector<std::pair<double,double> > z, Weight weight = arctan_weight){
    int md = diag.size(); std::vector<std::pair<double,double> > b; int mz = z.size();
    for(int i = 0; i < mz; i++){
      double d1 = 0; double d2 = 0; double zx = z[i].first; double zy = z[i].second;
      for(int j = 0; j < md; j++){
        double x = diag[j].first; double y = diag[j].second;
        d1 += weight(diag[j])*cos(x*zx + y*zy);
        d2 += weight(diag[j])*sin(x*zx + y*zy);
      }
      b.emplace_back(d1,d2);
    }
    return b;
  }

  std::vector<std::pair<double,double> > random_Fourier(double sigma, int m = 1000){
    std::normal_distribution<double> distrib(0,1); std::vector<std::pair<double,double> > z; std::random_device rd;
    for(int i = 0; i < m; i++){
      std::mt19937 e1(rd()); std::mt19937 e2(rd());
      double zx = distrib(e1); double zy = distrib(e2);
      z.emplace_back(zx/sigma,zy/sigma);
    }
    return z;
  }



  // **********************************
  // Scalar product + distance.
  // **********************************


  double compute_scalar_product(Persistence_weighted_gaussian second){

    PD diagram1 = this->diagram; PD diagram2 = second.diagram;

    if(this->approx == -1){
      int num_pts1 = diagram1.size(); int num_pts2 = diagram2.size(); double k = 0;
      for(int i = 0; i < num_pts1; i++)
        for(int j = 0; j < num_pts2; j++)
          k += this->weight(diagram1[i])*this->weight(diagram2[j])*exp(-((diagram1[i].first  - diagram2[j].first)  *  (diagram1[i].first  - diagram2[j].first) +
                                                                         (diagram1[i].second - diagram2[j].second) *  (diagram1[i].second - diagram2[j].second))
                                                                        /(2*this->sigma*this->sigma));
      return k;
    }
    else{
      std::vector<std::pair<double,double> > z =  random_Fourier(this->sigma, this->approx);
      std::vector<std::pair<double,double> > b1 = Fourier_feat(diagram1,z,this->weight);
      std::vector<std::pair<double,double> > b2 = Fourier_feat(diagram2,z,this->weight);
      double d = 0; for(int i = 0; i < this->approx; i++) d += b1[i].first*b2[i].first + b1[i].second*b2[i].second;
      return d/this->approx;
    }
  }

  double distance(Persistence_weighted_gaussian second, double power = 1) {
    if(this->sigma != second.get_sigma() || this->approx != second.get_approx()){
      std::cout << "Error: different representations!" << std::endl; return 0;
    }
    else return std::pow(this->compute_scalar_product(*this) + second.compute_scalar_product(second)-2*this->compute_scalar_product(second),  power/2.0);
  }


};

}  // namespace Persistence_weighted_gaussian
}  // namespace Gudhi

#endif  // PERSISTENCE_WEIGHTED_GAUSSIAN_H_
