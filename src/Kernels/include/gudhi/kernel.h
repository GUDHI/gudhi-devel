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

#ifndef KERNEL_H_
#define KERNEL_H_

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
//#include <gmp.h>
//#include <gmpxx.h>
#include <random>
#include <chrono>
#include <ctime>

using PD = std::vector<std::pair<double,double> >;
bool sortAngle(const std::pair<double, std::pair<int,int> >& p1, const std::pair<double, std::pair<int,int> >& p2){return (p1.first < p2.first);}
bool myComp(const std::pair<int,double> & P1, const std::pair<int,double> & P2){return P1.second < P2.second;}

namespace Gudhi {
namespace kernel {


double pss_weight(std::pair<double,double> P){
  if(P.second > P.first)  return 1;
  else return -1;
}




// ********************************************************************
// Exact computation.
// ********************************************************************

/** \brief Computes the Linear Persistence Weighted Gaussian Kernel between two persistence diagrams.
 *  \ingroup kernel
 *
 * @param[in] PD1      first persistence diagram.
 * @param[in] PD2      second persistence diagram.
 * @param[in] sigma    bandwidth parameter of the Gaussian Kernel used for the Kernel Mean Embedding of the diagrams.
 * @param[in] weight   weight function for the points in the diagrams.
 *
 */
double lpwgk(PD PD1, PD PD2, double sigma, double (*weight)(std::pair<double,double>) = [](std::pair<double,double> P){return atan(P.second - P.first);}){
  int num_pts1 = PD1.size(); int num_pts2 = PD2.size(); double k = 0;
  for(int i = 0; i < num_pts1; i++)
    for(int j = 0; j < num_pts2; j++)
      k += (*weight)(PD1[i])*(*weight)(PD2[j])*exp(-(pow(PD1[i].first-PD2[j].first,2) + pow(PD1[i].second-PD2[j].second,2))/(2*pow(sigma,2)));
  return k;
}

/** \brief Computes the Persistence Scale Space Kernel between two persistence diagrams.
 * \ingroup kernel
 *
 * @param[in] PD1      first persistence diagram.
 * @param[in] PD2      second persistence diagram.
 * @param[in] sigma    bandwidth parameter of the Gaussian Kernel used for the Kernel Mean Embedding of the diagrams.
 *
 */
double pssk(PD PD1, PD PD2, double sigma){
  PD pd1 = PD1; int numpts = PD1.size();    for(int i = 0; i < numpts; i++)  pd1.push_back(std::pair<double,double>(PD1[i].second,PD1[i].first));
  PD pd2 = PD2;     numpts = PD2.size();    for(int i = 0; i < numpts; i++)  pd2.push_back(std::pair<double,double>(PD2[i].second,PD2[i].first));
  return lpwgk(pd1, pd2, 2*sqrt(sigma), &pss_weight) / (2*8*3.14159265359*sigma);
}

/** \brief Computes the Gaussian Persistence Weighted Gaussian Kernel between two persistence diagrams.
 * \ingroup kernel
 *
 * @param[in] PD1      first persistence diagram.
 * @param[in] PD2      second persistence diagram.
 * @param[in] sigma    bandwidth parameter of the Gaussian Kernel used for the Kernel Mean Embedding of the diagrams.
 * @param[in] tau      bandwidth parameter of the Gaussian Kernel used between the embeddings.
 * @param[in] weight   weight function for the points in the diagrams.
 *
 */
double gpwgk(PD PD1, PD PD2, double sigma, double tau, double (*weight)(std::pair<double,double>) = [](std::pair<double,double> P){return atan(P.second - P.first);}){
  double k1 = lpwgk(PD1,PD1,sigma,weight);
  double k2 = lpwgk(PD2,PD2,sigma,weight);
  double k3 = lpwgk(PD1,PD2,sigma,weight);
  return exp( - (k1+k2-2*k3) / (2*pow(tau,2))  );
}

/** \brief Computes the RKHS distance induced by the Gaussian Kernel Embedding between two persistence diagrams.
 * \ingroup kernel
 *
 * @param[in] PD1      first persistence diagram.
 * @param[in] PD2      second persistence diagram.
 * @param[in] sigma    bandwidth parameter of the Gaussian Kernel used for the Kernel Mean Embedding of the diagrams.
 * @param[in] weight   weight function for the points in the diagrams.
 *
 */
double dpwg(PD PD1, PD PD2, double sigma, double (*weight)(std::pair<double,double>) = [](std::pair<double,double> P){return atan(P.second - P.first);}){
  double k1 = lpwgk(PD1,PD1,sigma,weight);
  double k2 = lpwgk(PD2,PD2,sigma,weight);
  double k3 = lpwgk(PD1,PD2,sigma,weight);
  return std::sqrt(k1+k2-2*k3);
}

// Compute the angle formed by two points of a PD
double compute_angle(const PD & PersDiag, const int & i, const int & j){
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

// Compute the integral of |cos()| between alpha and beta
// Valid only if alpha is in [-pi,pi] and beta-alpha is in [0,pi]
double compute_int_cos(const double & alpha, const double & beta){
  double res = 0;
  assert((alpha >= 0 && alpha <= 3.14159265359) || (alpha >= -3.14159265359 && alpha <= 0));
  if (alpha >= 0 && alpha <= 3.14159265359){
    if (cos(alpha) >= 0){
      if(3.14159265359/2 <= beta){res = 2-sin(alpha)-sin(beta);}
      else{res = sin(beta)-sin(alpha);}
    }
    else{
      if(1.5*3.14159265359 <= beta){res = 2+sin(alpha)+sin(beta);}
      else{res = sin(alpha)-sin(beta);}
    }
  }
  if (alpha >= -3.14159265359 && alpha <= 0){
    if (cos(alpha) <= 0){
      if(-3.14159265359/2 <= beta){res = 2+sin(alpha)+sin(beta);}
      else{res = sin(alpha)-sin(beta);}
    }
    else{
      if(3.14159265359/2 <= beta){res = 2-sin(alpha)-sin(beta);}
      else{res = sin(beta)-sin(alpha);}
    }
  }
  return res;
}

double compute_int(const double & theta1, const double & theta2, const int & p, const int & q, const PD & PD1, const PD & PD2){
  double norm = std::sqrt(pow(PD1[p].first-PD2[q].first,2) + pow(PD1[p].second-PD2[q].second,2));
  double angle1;
  if (PD1[p].first > PD2[q].first)
    angle1 = theta1 - asin( (PD1[p].second-PD2[q].second)/norm  );
  else
    angle1 = theta1 - asin( (PD2[q].second-PD1[p].second)/norm  );
  double angle2 = angle1 + theta2 - theta1;
  double integral = compute_int_cos(angle1,angle2);
  return norm*integral;
}



double compute_sw(const std::vector<std::vector<std::pair<int,double> > > & V1, const std::vector<std::vector<std::pair<int,double> > > & V2, const PD & PD1, const PD & PD2){
  int N = V1.size(); double sw = 0;
  for (int i = 0; i < N; i++){
    std::vector<std::pair<int,double> > U,V; U = V1[i]; V = V2[i];
    double theta1, theta2; theta1 = -3.14159265359/2;
    unsigned int ku, kv; ku = 0; kv = 0; theta2 = std::min(U[ku].second,V[kv].second);
    while(theta1 != 3.14159265359/2){
      if(PD1[U[ku].first].first != PD2[V[kv].first].first || PD1[U[ku].first].second != PD2[V[kv].first].second)
        if(theta1 != theta2)
          sw += compute_int(theta1, theta2, U[ku].first, V[kv].first, PD1, PD2);
      theta1 = theta2;
      if (  (theta2 == U[ku].second)  &&  ku < U.size()-1  )  ku++;
      if (  (theta2 == V[kv].second)  &&  kv < V.size()-1  )  kv++;
      theta2 = std::min(U[ku].second, V[kv].second);
    }
  }
  return sw/3.14159265359;
}

/** \brief Computes the Sliced Wasserstein distance between two persistence diagrams.
 * \ingroup kernel
 *
 * @param[in] PD1 first persistence diagram.
 * @param[in] PD2 second persistence diagram.
 *
 */
 double sw(PD PD1, PD PD2){

  // Add projections onto diagonal.
  int n1, n2; n1 = PD1.size(); n2 = PD2.size(); double max_ordinate = std::numeric_limits<double>::lowest();
  for (int i = 0; i < n2; i++){
    max_ordinate = std::max(max_ordinate, PD2[i].second);
    PD1.push_back(  std::pair<double,double>(  ((PD2[i].first+PD2[i].second)/2), ((PD2[i].first+PD2[i].second)/2)  )   );
  }
  for (int i = 0; i < n1; i++){
    max_ordinate = std::max(max_ordinate, PD1[i].second);
    PD2.push_back(  std::pair<double,double>(  ((PD1[i].first+PD1[i].second)/2), ((PD1[i].first+PD1[i].second)/2)  )   );
  }
  int N = PD1.size(); assert(N==PD2.size());

  // Slightly perturb the points so that the PDs are in generic positions.
  int mag = 0; while(max_ordinate > 10){mag++; max_ordinate/=10;}
  double thresh = pow(10,-5+mag);
  srand(time(NULL));
  for (int i = 0; i < N; i++){
    PD1[i].first += thresh*(1.0-2.0*rand()/RAND_MAX); PD1[i].second += thresh*(1.0-2.0*rand()/RAND_MAX);
    PD2[i].first += thresh*(1.0-2.0*rand()/RAND_MAX); PD2[i].second += thresh*(1.0-2.0*rand()/RAND_MAX);
  }

  // Compute all angles in both PDs.
  std::vector<std::pair<double, std::pair<int,int> > > angles1, angles2;
  for (int i = 0; i < N; i++){
    for (int j = i+1; j < N; j++){
      double theta1 = compute_angle(PD1,i,j); double theta2 = compute_angle(PD2,i,j);
      angles1.push_back(std::pair<double, std::pair<int,int> >(theta1, std::pair<int,int>(i,j)));
      angles2.push_back(std::pair<double, std::pair<int,int> >(theta2, std::pair<int,int>(i,j)));
    }
  }

  // Sort angles.
  std::sort(angles1.begin(), angles1.end(), sortAngle); std::sort(angles2.begin(), angles2.end(), sortAngle);

  // Initialize orders of the points of both PDs (given by ordinates when theta = -pi/2).
  std::vector<int> orderp1, orderp2;
  for (int i = 0; i < N; i++){ orderp1.push_back(i); orderp2.push_back(i); }
  std::sort( orderp1.begin(), orderp1.end(), [=](int i, int j){ if(PD1[i].second != PD1[j].second) return (PD1[i].second < PD1[j].second); else return (PD1[i].first > PD1[j].first); } );
  std::sort( orderp2.begin(), orderp2.end(), [=](int i, int j){ if(PD2[i].second != PD2[j].second) return (PD2[i].second < PD2[j].second); else return (PD2[i].first > PD2[j].first); } );

  // Find the inverses of the orders.
  std::vector<int> order1(N); std::vector<int> order2(N);
  for(int i = 0; i < N; i++)  for (int j = 0; j < N; j++)  if(orderp1[j] == i){  order1[i] = j; break;  }
  for(int i = 0; i < N; i++)  for (int j = 0; j < N; j++)  if(orderp2[j] == i){  order2[i] = j; break;  }

  // Record all inversions of points in the orders as theta varies along the positive half-disk.
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
    anglePerm1[order1[i]].push_back(std::pair<int, double>(i,3.14159265359/2));
    anglePerm2[order2[i]].push_back(std::pair<int, double>(i,3.14159265359/2));
  }

  // Compute the SW distance with the list of inversions.
  return compute_sw(anglePerm1, anglePerm2, PD1, PD2);

}

 /** \brief Computes the Sliced Wasserstein Kernel between two persistence diagrams.
  * \ingroup kernel
  *
  * @param[in] PD1     first persistence diagram.
  * @param[in] PD2     second persistence diagram.
  * @param[in] sigma   bandwidth parameter.
  *
  */
  double swk(PD PD1, PD PD2, double sigma){
    return exp( - sw(PD1,PD2) / (2*pow(sigma, 2)) );
  }








// ********************************************************************
// Approximate computation.
// ********************************************************************

double approx_lpwg_Fourier(const std::vector<std::pair<double,double> >& B1, const std::vector<std::pair<double,double> >& B2){
  double d = 0; int M  = B1.size();
  for(int i = 0; i < M; i++) d += B1[i].first*B2[i].first + B1[i].second*B2[i].second;
  return (1.0/M)*d;
}

double approx_gpwg_Fourier(const std::vector<std::pair<double,double> >& B1, const std::vector<std::pair<double,double> >& B2, double tau){
  int M  = B1.size();
  double d3 = approx_lpwg_Fourier(B1, B2);
  double d1 = 0; double d2 = 0;
  for(int i = 0; i < M; i++){d1 += pow(B1[i].first,2) + pow(B1[i].second,2); d2 += pow(B2[i].first,2) + pow(B2[i].second,2);}
  return exp( -((1.0/M)*(d1+d2)-2*d3) / (2*pow(tau,2))  );
}

double approx_dpwg_Fourier(const std::vector<std::pair<double,double> >& B1, const std::vector<std::pair<double,double> >& B2){
  int M  = B1.size();
  double d3 = approx_lpwg_Fourier(B1, B2);
  double d1 = 0; double d2 = 0;
  for(int i = 0; i < M; i++){d1 += pow(B1[i].first,2) + pow(B1[i].second,2); d2 += pow(B2[i].first,2) + pow(B2[i].second,2);}
  return std::sqrt((1.0/M)*(d1+d2)-2*d3);
}

std::vector<std::pair<double,double> > Fourier_feat(PD D, std::vector<std::pair<double,double> > Z, double (*weight)(std::pair<double,double>) = [](std::pair<double,double> P){return atan(P.second - P.first);}){
  int m = D.size(); std::vector<std::pair<double,double> > B; int M = Z.size();
  for(int i = 0; i < M; i++){
    double d1 = 0; double d2 = 0; double zx = Z[i].first; double zy = Z[i].second;
    for(int j = 0; j < m; j++){
      double x = D[j].first; double y = D[j].second;
      d1 += (*weight)(D[j])*cos(x*zx + y*zy);
      d2 += (*weight)(D[j])*sin(x*zx + y*zy);
    }
    B.push_back(std::pair<double,double>(d1,d2));
  }
  return B;
}

std::vector<std::pair<double,double> > random_Fourier(double sigma, int M = 1000){
  std::normal_distribution<double> distrib(0,1); std::vector<std::pair<double,double> > Z; std::random_device rd;
  for(int i = 0; i < M; i++){
    std::mt19937 e1(rd()); std::mt19937 e2(rd());
    double zx = distrib(e1); double zy = distrib(e2);
    Z.push_back(std::pair<double,double>((1.0/sigma)*zx,(1.0/sigma)*zy));
  }
  return Z;
}


/** \brief Computes an approximation of the Linear Persistence Weighted Gaussian Kernel between two persistence diagrams with random Fourier features.
 * \ingroup kernel
 *
 * @param[in] PD1       first persistence diagram.
 * @param[in] PD2       second persistence diagram.
 * @param[in] sigma     bandwidth parameter of the Gaussian Kernel used for the Kernel Mean Embedding of the diagrams.
 * @param[in] weight    weight function for the points in the diagrams.
 * @param[in] M         number of Fourier features.
 *
 */
double approx_lpwgk(PD PD1, PD PD2, double sigma, double (*weight)(std::pair<double,double>) = [](std::pair<double,double> P){return atan(P.second - P.first);}, int M = 1000){
  std::vector<std::pair<double,double> > Z =  random_Fourier(sigma, M);
  std::vector<std::pair<double,double> > B1 = Fourier_feat(PD1,Z,weight);
  std::vector<std::pair<double,double> > B2 = Fourier_feat(PD2,Z,weight);
  return approx_lpwg_Fourier(B1,B2);
}

/** \brief Computes an approximation of the Persistence Scale Space Kernel between two persistence diagrams with random Fourier features.
 * \ingroup kernel
 *
 * @param[in] PD1       first persistence diagram.
 * @param[in] PD2       second persistence diagram.
 * @param[in] sigma     bandwidth parameter of the Gaussian Kernel used for the Kernel Mean Embedding of the diagrams.
 * @param[in] M         number of Fourier features.
 *
 */
double approx_pssk(PD PD1, PD PD2, double sigma, int M = 1000){
  PD pd1 = PD1; int numpts = PD1.size();    for(int i = 0; i < numpts; i++)  pd1.push_back(std::pair<double,double>(PD1[i].second,PD1[i].first));
  PD pd2 = PD2;     numpts = PD2.size();    for(int i = 0; i < numpts; i++)  pd2.push_back(std::pair<double,double>(PD2[i].second,PD2[i].first));
  return approx_lpwgk(pd1, pd2, 2*sqrt(sigma), &pss_weight, M) / (2*8*3.14159265359*sigma);
}


/** \brief Computes an approximation of the Gaussian Persistence Weighted Gaussian Kernel between two persistence diagrams with random Fourier features.
 * \ingroup kernel
 *
 * @param[in] PD1       first persistence diagram.
 * @param[in] PD2       second persistence diagram.
 * @param[in] sigma     bandwidth parameter of the Gaussian Kernel used for the Kernel Mean Embedding of the diagrams.
 * @param[in] tau       bandwidth parameter of the Gaussian Kernel used between the embeddings.
 * @param[in] weight    weight function for the points in the diagrams.
 * @param[in] M         number of Fourier features.
 *
 */
double approx_gpwgk(PD PD1, PD PD2, double sigma, double tau, double (*weight)(std::pair<double,double>) = [](std::pair<double,double> P){return atan(P.second - P.first);}, int M = 1000){
  std::vector<std::pair<double,double> > Z =  random_Fourier(sigma, M);
  std::vector<std::pair<double,double> > B1 = Fourier_feat(PD1,Z,weight);
  std::vector<std::pair<double,double> > B2 = Fourier_feat(PD2,Z,weight);
  return approx_gpwg_Fourier(B1,B2,tau);
}


/** \brief Computes an approximation of the Sliced Wasserstein distance between two persistence diagrams.
 * \ingroup kernel
 *
 * @param[in] PD1       first persistence diagram.
 * @param[in] PD2       second persistence diagram.
 * @param[in] N         number of points sampled on the circle.
 *
 */
double approx_sw(PD PD1, PD PD2, int N = 100){

  double step = 3.14159265359/N; double sw = 0;

  // Add projections onto diagonal.
  int n1, n2; n1 = PD1.size(); n2 = PD2.size();
  for (int i = 0; i < n2; i++)
    PD1.push_back(std::pair<double,double>(   (PD2[i].first + PD2[i].second)/2,  (PD2[i].first + PD2[i].second)/2)   );
  for (int i = 0; i < n1; i++)
    PD2.push_back(std::pair<double,double>(   (PD1[i].first + PD1[i].second)/2,  (PD1[i].first + PD1[i].second)/2)   );
  int n = PD1.size();

  // Sort and compare all projections.
  //#pragma omp parallel for
  for (int i = 0; i < N; i++){
    std::vector<std::pair<int,double> > L1, L2;
    for (int j = 0; j < n; j++){
      L1.push_back(   std::pair<int,double>(j, PD1[j].first*cos(-3.14159265359/2+i*step) + PD1[j].second*sin(-3.14159265359/2+i*step))   );
      L2.push_back(   std::pair<int,double>(j, PD2[j].first*cos(-3.14159265359/2+i*step) + PD2[j].second*sin(-3.14159265359/2+i*step))   );
    }
    std::sort(L1.begin(),L1.end(), myComp); std::sort(L2.begin(),L2.end(), myComp);
    double f = 0; for (int j = 0; j < n; j++)  f += std::abs(L1[j].second - L2[j].second);
    sw += f*step;
  }
  return sw/3.14159265359;
}

/** \brief Computes an approximation of the Sliced Wasserstein Kernel between two persistence diagrams.
 * \ingroup kernel
 *
 * @param[in] PD1       first persistence diagram.
 * @param[in] PD2       second persistence diagram.
 * @param[in] sigma     bandwidth parameter.
 * @param[in] N         number of points sampled on the circle.
 *
 */
double approx_swk(PD PD1, PD PD2, double sigma, int N = 100){
  return exp( - approx_sw(PD1,PD2,N) / (2*pow(sigma,2)));
}



} // namespace kernel

} // namespace Gudhi

#endif //KERNEL_H_
