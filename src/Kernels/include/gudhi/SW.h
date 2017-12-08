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

#ifndef SW_H_
#define SW_H_

#define NUMPI 3.14159265359

#include <stdlib.h>
#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <set>
#include <map>
#include <limits>
#include <cmath>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <omp.h>
#include <assert.h>
#include <iomanip>
#include <gmp.h>
#include <gmpxx.h>
#include <random>
#include <chrono>
#include <ctime>
#include <math.h>

using PD = std::vector<std::pair<double,double> >;

std::vector<std::pair<double,double> > PDi, PDj;

bool compOri(const int& p, const int& q){
  if(PDi[p].second != PDi[q].second)
    return (PDi[p].second < PDi[q].second);
  else
    return (PDi[p].first > PDi[q].first);
}

bool compOrj(const int& p, const int& q){
  if(PDj[p].second != PDj[q].second)
    return (PDj[p].second < PDj[q].second);
  else
    return (PDj[p].first > PDj[q].first);
}

bool sortAngle(const std::pair<double, std::pair<int,int> >& p1, const std::pair<double, std::pair<int,int> >& p2){
  return p1.first < p2.first;
}

bool myComp(const std::pair<int,double> & P1, const std::pair<int,double> & P2){return P1.second < P2.second;}

namespace Gudhi {
namespace sliced_wasserstein {


double compute_approximate_SW(PD PD1, PD PD2, int N = 100){

  double step = NUMPI/N; double sw = 0;

  // Add projections onto diagonal.
  // ******************************
  int n1, n2; n1 = PD1.size(); n2 = PD2.size();
  for (int i = 0; i < n2; i++)
    PD1.push_back(std::pair<double,double>(   (PD2[i].first+PD2[i].second)/2,  (PD2[i].first+PD2[i].second)/2)   );
  for (int i = 0; i < n1; i++)
    PD2.push_back(std::pair<double,double>(   (PD1[i].first+PD1[i].second)/2,  (PD1[i].first+PD1[i].second)/2)   );
  int n = PD1.size();

  // Sort and compare all projections.
  // *********************************
  //#pragma omp parallel for
  for (int i = 0; i < N; i++){
    std::vector<std::pair<int,double> > L1, L2;
    for (int j = 0; j < n; j++){
      L1.push_back(   std::pair<int,double>(j, PD1[j].first*cos(-NUMPI/2+i*step) + PD1[j].second*sin(-NUMPI/2+i*step))   );
      L2.push_back(   std::pair<int,double>(j, PD2[j].first*cos(-NUMPI/2+i*step) + PD2[j].second*sin(-NUMPI/2+i*step))   );
    }
    std::sort(L1.begin(),L1.end(), myComp); std::sort(L2.begin(),L2.end(), myComp);
    double f = 0; for (int j = 0; j < n; j++)  f += std::abs(L1[j].second - L2[j].second);
    sw += f*step;
  }
  return sw/NUMPI;
}

double compute_int_cos(const double& alpha, const double& beta){ // Valid only if alpha is in [-pi,pi] and beta-alpha is in [0,pi]
  double res;
  assert((alpha >= 0 && alpha <= NUMPI) || (alpha >= -NUMPI && alpha <= 0));
  if (alpha >= 0 && alpha <= NUMPI){
    if (cos(alpha) >= 0){
      if(NUMPI/2 <= beta){res = 2-sin(alpha)-sin(beta);}
      else{res = sin(beta)-sin(alpha);}
    }
    else{
      if(1.5*NUMPI <= beta){res = 2+sin(alpha)+sin(beta);}
      else{res = sin(alpha)-sin(beta);}
    }
  }
  if (alpha >= -NUMPI && alpha <= 0){
    if (cos(alpha) <= 0){
      if(-NUMPI/2 <= beta){res = 2+sin(alpha)+sin(beta);}
      else{res = sin(alpha)-sin(beta);}
    }
    else{
      if(NUMPI/2 <= beta){res = 2-sin(alpha)-sin(beta);}
      else{res = sin(beta)-sin(alpha);}
    }
  }
  return res;
}

double compute_int(const double& theta1, const double& theta2, const int& p, const int& q){
  double norm = std::sqrt(pow(PDi[p].first-PDj[q].first,2) + pow(PDi[p].second-PDj[q].second,2));
  double angle1;
  if (PDi[p].first > PDj[q].first)
    angle1 = theta1 - asin( (PDi[p].second-PDj[q].second)/norm  );
  else
    angle1 = theta1 - asin( (PDj[q].second-PDi[p].second)/norm  );
  double angle2 = angle1+theta2-theta1;
  double integral = compute_int_cos(angle1,angle2);
  return norm*integral;
}

double compute_sw(const std::vector<std::vector<std::pair<int,double> > >& V1, \
                  const std::vector<std::vector<std::pair<int,double> > >& V2){
  int N = V1.size(); double sw = 0;
  for (int i = 0; i < N; i++){
    std::vector<std::pair<int,double> > U,V; U = V1[i]; V = V2[i];
    double theta1, theta2; theta1 = -NUMPI/2;
    int ku, kv; ku = 0; kv = 0; theta2 = std::min(U[ku].second,V[kv].second);
    while(theta1 != NUMPI/2){
      if(PDi[U[ku].first].first != PDj[V[kv].first].first || PDi[U[ku].first].second != PDj[V[kv].first].second)
        if(theta1 != theta2)
          sw += compute_int(theta1,theta2,U[ku].first,V[kv].first);
      theta1 = theta2;
      if (  (theta2 == U[ku].second)  &&  ku < U.size()-1  ){ku++;}
      if (  (theta2 == V[kv].second)  &&  kv < V.size()-1  ){kv++;}
      theta2 = std::min(U[ku].second, V[kv].second);
    }
  }
  return sw/NUMPI;
}

double compute_angle(const PD& PersDiag, const int& i, const int& j){
  std::pair<double,double> vect; double x1,y1, x2,y2;
  x1 = PersDiag[i].first; y1 = PersDiag[i].second;
  x2 = PersDiag[j].first; y2 = PersDiag[j].second;
  if (y1 - y2 > 0){
    vect.first = y1 - y2;
    vect.second = x2 - x1;}
  else{
    if(y1 - y2 < 0){
      vect.first = y2 - y1;
      vect.second = x1 - x2;
    }
    else{
      vect.first = 0;
      vect.second = abs(x1 - x2);}
  }
  double norm = std::sqrt(pow(vect.first,2) + pow(vect.second,2));
  return asin(vect.second/norm);
}

double compute_exact_SW(PD PD1, PD PD2){

  // Add projections onto diagonal.
  // ******************************
  int n1, n2; n1 = PD1.size(); n2 = PD2.size(); double max_ordinate = std::numeric_limits<double>::min();
  for (int i = 0; i < n2; i++){
    max_ordinate = std::max(max_ordinate, PD2[i].second);
    PD1.push_back(std::pair<double,double>(  ((PD2[i].first+PD2[i].second)/2), ((PD2[i].first+PD2[i].second)/2))   );
  }
  for (int i = 0; i < n1; i++){
    max_ordinate = std::max(max_ordinate, PD1[i].second);
    PD2.push_back(std::pair<double,double>(  ((PD1[i].first+PD1[i].second)/2), ((PD1[i].first+PD1[i].second)/2))   );
  }
  int N = PD1.size(); assert(N==PD2.size());

  // Slightly perturb the points so that the PDs are in generic positions.
  // *********************************************************************
  int mag = 0; while(max_ordinate > 10){mag++; max_ordinate/=10;}
  double thresh = pow(10,-5+mag);
  srand(time(NULL));
  for (int i = 0; i < N; i++){
    PD1[i].first += thresh*(1.0-2.0*rand()/RAND_MAX); PD1[i].second += thresh*(1.0-2.0*rand()/RAND_MAX);
    PD2[i].first += thresh*(1.0-2.0*rand()/RAND_MAX); PD2[i].second += thresh*(1.0-2.0*rand()/RAND_MAX);
  }

  // Compute all angles in both PDs.
  // *******************************
  std::vector<std::pair<double, std::pair<int,int> > > angles1, angles2;
  for (int i = 0; i < N; i++){
    for (int j = i+1; j < N; j++){
      double theta1 = compute_angle(PD1,i,j); double theta2 = compute_angle(PD2,i,j);
      angles1.push_back(std::pair<double, std::pair<int,int> >(theta1, std::pair<int,int>(i,j)));
      angles2.push_back(std::pair<double, std::pair<int,int> >(theta2, std::pair<int,int>(i,j)));
    }
  }

  // Sort angles.
  // ************
  std::sort(angles1.begin(), angles1.end(), sortAngle); std::sort(angles2.begin(), angles2.end(), sortAngle);

  // Initialize orders of the points of both PD (given by ordinates when theta = -pi/2).
  // ***********************************************************************************
  PDi = PD1; PDj = PD2;
  std::vector<int> orderp1, orderp2;
  for (int i = 0; i < N; i++){orderp1.push_back(i); orderp2.push_back(i);}
  std::sort(orderp1.begin(),orderp1.end(),compOri); std::sort(orderp2.begin(),orderp2.end(),compOrj);

  // Find the inverses of the orders.
  // ********************************
  std::vector<int> order1(N); std::vector<int> order2(N);
  for(int i = 0; i < N; i++){
    for (int j = 0; j < N; j++)
      if(orderp1[j] == i)
        order1[i] = j;
  }
  for(int i = 0; i < N; i++){
    for (int j = 0; j < N; j++)
      if(orderp2[j] == i)
        order2[i] = j;
  }

  // Record all inversions of points in the orders as theta varies along the positive half-disk.
  // *******************************************************************************************
  std::vector<std::vector<std::pair<int,double> > > anglePerm1(N);
  std::vector<std::vector<std::pair<int,double> > > anglePerm2(N);

  int M1 = angles1.size();
  for (int i = 0; i < M1; i++){
    double theta = angles1[i].first; int p = angles1[i].second.first; int q = angles1[i].second.second;
    anglePerm1[order1[p]].push_back(std::pair<int, double>(p,theta));
    anglePerm1[order1[q]].push_back(std::pair<int, double>(q,theta));
    int a = order1[p]; int b = order1[q]; order1[p] = b; order1[q] = a;
  }

  int M2 = angles2.size();
  for (int i = 0; i < M2; i++){
    double theta = angles2[i].first; int p = angles2[i].second.first; int q = angles2[i].second.second;
    anglePerm2[order2[p]].push_back(std::pair<int, double>(p,theta));
    anglePerm2[order2[q]].push_back(std::pair<int, double>(q,theta));
    int a = order2[p]; int b = order2[q]; order2[p] = b; order2[q] = a;
  }

  for (int i = 0; i < N; i++){
    anglePerm1[order1[i]].push_back(std::pair<int, double>(i,NUMPI/2));
    anglePerm2[order2[i]].push_back(std::pair<int, double>(i,NUMPI/2));
  }

  // Compute the SW distance with the list of inversions.
  // ****************************************************
  return compute_sw(anglePerm1,anglePerm2);

}

}}

#endif


