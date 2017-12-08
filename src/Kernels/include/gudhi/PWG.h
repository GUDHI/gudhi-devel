/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Mathieu Carrière
 *
 *    Copyright (C) 2017  INRIA (France)
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

#ifndef PWG_H_
#define PWG_H_

#define NUMPI 3.14159265359

#include <stdlib.h>
#include <stdio.h>
#include <cstdlib>
#include <cstdio>
#include <iomanip>

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include <limits>
#include <assert.h>

#include <cmath>
#include <math.h>

#include <memory>
#include <stdexcept>
#include <omp.h>

#include <gmp.h>
#include <gmpxx.h>
#include <random>
#include <chrono>
#include <ctime>

using PD = std::vector<std::pair<double,double> >;

namespace Gudhi {
namespace persistence_weighted_gaussian {

double compute_exact_linear_pwg(PD PD1, PD PD2, double sigma, double C, int p){

  int num_pts1 = PD1.size();
  int num_pts2 = PD2.size();

  double k = 0;
  for(int i = 0; i < num_pts1; i++){
    for(int j = 0; j < num_pts2; j++){
      k += atan(C*pow(PD1[i].second-PD1[i].first,p))*atan(C*pow(PD2[j].second-PD2[j].first,p))*\
            exp(  -(    pow(PD1[i].first-PD2[j].first,2) + pow(PD1[i].second-PD2[j].second,2)     )/(2*pow(sigma,2))  );
    }
  }

  return k;

}

double compute_exact_gaussian_pwg(PD PD1, PD PD2, double sigma, double C, int p, double tau){

  int num_pts1 = PD1.size();
  int num_pts2 = PD2.size();

  double k1 = 0;
  for(int i = 0; i < num_pts1; i++){
    for(int j = 0; j < num_pts1; j++){
      k1 += atan(C*pow(PD1[i].second-PD1[i].first,p))*atan(C*pow(PD1[j].second-PD1[j].first,p))*\
            exp(  -(    pow(PD1[i].first-PD1[j].first,2) + pow(PD1[i].second-PD1[j].second,2)     )/(2*pow(sigma,2))  );
    }
  }

  double k2 = 0;
  for(int i = 0; i < num_pts2; i++){
    for(int j = 0; j < num_pts2; j++){
      k2 += atan(C*pow(PD2[i].second-PD2[i].first,p))*atan(C*pow(PD2[j].second-PD2[j].first,p))*\
            exp(  -(    pow(PD2[i].first-PD2[j].first,2) + pow(PD2[i].second-PD2[j].second,2)     )/(2*pow(sigma,2))  );
    }
  }

  double k3 = compute_exact_linear_pwg(PD1,PD2,sigma,C,p);
  return exp( - (k1+k2-2*k3) / (2*pow(tau,2))  );

}

double compute_exact_gaussian_RKHSdist(PD PD1, PD PD2, double sigma, double C, int p){

  int num_pts1 = PD1.size();
  int num_pts2 = PD2.size();

  double k1 = 0;
  for(int i = 0; i < num_pts1; i++){
    for(int j = 0; j < num_pts1; j++){
      k1 += atan(C*pow(PD1[i].second-PD1[i].first,p))*atan(C*pow(PD1[j].second-PD1[j].first,p))*\
            exp(  -(    pow(PD1[i].first-PD1[j].first,2) + pow(PD1[i].second-PD1[j].second,2)     )/(2*pow(sigma,2))  );
    }
  }

  double k2 = 0;
  for(int i = 0; i < num_pts2; i++){
    for(int j = 0; j < num_pts2; j++){
      k2 += atan(C*pow(PD2[i].second-PD2[i].first,p))*atan(C*pow(PD2[j].second-PD2[j].first,p))*\
            exp(  -(    pow(PD2[i].first-PD2[j].first,2) + pow(PD2[i].second-PD2[j].second,2)     )/(2*pow(sigma,2))  );
    }
  }

  double k3 = compute_exact_linear_pwg(PD1,PD2,sigma,C,p);
  return std::sqrt(k1+k2-2*k3);

}

double compute_approximate_linear_pwg_from_Fourier_features(const std::vector<std::pair<double,double> >& B1, \
                                                            const std::vector<std::pair<double,double> >& B2){
  double d = 0; int M  = B1.size();
  for(int i = 0; i < M; i++) d += B1[i].first*B2[i].first + B1[i].second*B2[i].second;
  return (1.0/M)*d;
}

double compute_approximate_gaussian_pwg_from_Fourier_features(const std::vector<std::pair<double,double> >& B1, \
                                                              const std::vector<std::pair<double,double> >& B2, double tau){
  int M  = B1.size();
  double d3 = compute_approximate_linear_pwg_from_Fourier_features(B1, B2);
  double d1 = 0; double d2 = 0;
  for(int i = 0; i < M; i++){d1 += pow(B1[i].first,2) + pow(B1[i].second,2); d2 += pow(B2[i].first,2) + pow(B2[i].second,2);}
  return exp( -((1.0/M)*(d1+d2)-2*d3) / (2*pow(tau,2))  );
}

double compute_approximate_gaussian_RKHSdist_from_Fourier_features(const std::vector<std::pair<double,double> >& B1, \
                                                              const std::vector<std::pair<double,double> >& B2){
  int M  = B1.size();
  double d3 = compute_approximate_linear_pwg_from_Fourier_features(B1, B2);
  double d1 = 0; double d2 = 0;
  for(int i = 0; i < M; i++){d1 += pow(B1[i].first,2) + pow(B1[i].second,2); d2 += pow(B2[i].first,2) + pow(B2[i].second,2);}
  return std::sqrt((1.0/M)*(d1+d2)-2*d3);
}

std::vector<std::pair<double,double> > compute_Fourier_features(double C, int p, PD D, std::vector<std::pair<double,double> > Z){
  int m = D.size(); std::vector<std::pair<double,double> > B; int M = Z.size();
  for(int i = 0; i < M; i++){
    double d1 = 0; double d2 = 0; double zx = Z[i].first; double zy = Z[i].second;
    for(int j = 0; j < m; j++){
      double x = D[j].first; double y = D[j].second;
      d1 += atan(C*pow(y-x,p))*cos(x*zx + y*zy);
      d2 += atan(C*pow(y-x,p))*sin(x*zx + y*zy);
    }
    B.push_back(std::pair<double,double>(d1,d2));
  }
  return B;
}

std::vector<std::pair<double,double> > random_Fourier(double sigma, int M = 1000){
  std::normal_distribution<double> distrib(0,1); std::vector<std::pair<double,double> > Z;
  std::random_device rd;
  for(int i = 0; i < M; i++){
    //unsigned seedx = 2*i; unsigned seedy = 2*i+1;
    //std::default_random_engine generatorx(seedx); std::default_random_engine generatory(seedy);
    std::mt19937 e1(rd()); std::mt19937 e2(rd());
    double zx = distrib(e1/*generatorx*/); double zy = distrib(e2/*generatory*/);
    Z.push_back(std::pair<double,double>((1/sigma)*zx,(1/sigma)*zy));
  }
  return Z;
}

double compute_approximate_linear_pwg(PD PD1, PD PD2, double sigma, double C, int p, int M = 1000){
  std::vector<std::pair<double,double> > Z = random_Fourier(sigma, M);
  std::vector<std::pair<double,double> > B1 = compute_Fourier_features(C,p,PD1,Z);
  std::vector<std::pair<double,double> > B2 = compute_Fourier_features(C,p,PD2,Z);
  return compute_approximate_linear_pwg_from_Fourier_features(B1,B2);
}

double compute_approximate_gaussian_pwg(PD PD1, PD PD2, double sigma, double C, int p, double tau, int M = 1000){
  std::vector<std::pair<double,double> > Z = random_Fourier(sigma, M);
  std::vector<std::pair<double,double> > B1 = compute_Fourier_features(C,p,PD1,Z);
  std::vector<std::pair<double,double> > B2 = compute_Fourier_features(C,p,PD2,Z);
  return compute_approximate_gaussian_pwg_from_Fourier_features(B1,B2,tau);
}


} // namespace persistence_weighted_gaussian

} // namespace Gudhi

#endif //PWG_H_
