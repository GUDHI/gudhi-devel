/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Mathieu Carrière
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

#ifndef KERNEL_H_
#define KERNEL_H_

#include <cstdlib>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include <limits>          //for numeric_limits<>
#include <utility>         //for pair<>

#include <boost/math/constants/constants.hpp>


namespace Gudhi {
namespace kernel {

using PD = std::vector<std::pair<double,double> >;
double pi = boost::math::constants::pi<double>();




// ********************************************************************
// Utils.
// ********************************************************************

bool sortAngle(const std::pair<double, std::pair<int,int> >& p1, const std::pair<double, std::pair<int,int> >& p2){return (p1.first < p2.first);}
bool myComp(const std::pair<int,double> & P1, const std::pair<int,double> & P2){return P1.second < P2.second;}

double pss_weight(std::pair<double,double> P){
  if(P.second > P.first)  return 1;
  else return -1;
}

double arctan_weight(std::pair<double,double> P){
  return atan(P.second - P.first);
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
  double norm = std::sqrt(vect.first*vect.first + vect.second*vect.second);
  return asin(vect.second/norm);
}

// Compute the integral of |cos()| between alpha and beta, valid only if alpha is in [-pi,pi] and beta-alpha is in [0,pi]
double compute_int_cos(const double & alpha, const double & beta){
  double res = 0;
  if (alpha >= 0 && alpha <= pi){
    if (cos(alpha) >= 0){
      if(pi/2 <= beta){res = 2-sin(alpha)-sin(beta);}
      else{res = sin(beta)-sin(alpha);}
    }
    else{
      if(1.5*pi <= beta){res = 2+sin(alpha)+sin(beta);}
      else{res = sin(alpha)-sin(beta);}
    }
  }
  if (alpha >= -pi && alpha <= 0){
    if (cos(alpha) <= 0){
      if(-pi/2 <= beta){res = 2+sin(alpha)+sin(beta);}
      else{res = sin(alpha)-sin(beta);}
    }
    else{
      if(pi/2 <= beta){res = 2-sin(alpha)-sin(beta);}
      else{res = sin(beta)-sin(alpha);}
    }
  }
  return res;
}

double compute_int(const double & theta1, const double & theta2, const int & p, const int & q, const PD & PD1, const PD & PD2){
  double norm = std::sqrt(  (PD1[p].first-PD2[q].first)*(PD1[p].first-PD2[q].first) + (PD1[p].second-PD2[q].second)*(PD1[p].second-PD2[q].second)  );
  double angle1;
  if (PD1[p].first > PD2[q].first)
    angle1 = theta1 - asin( (PD1[p].second-PD2[q].second)/norm  );
  else
    angle1 = theta1 - asin( (PD2[q].second-PD1[p].second)/norm  );
  double angle2 = angle1 + theta2 - theta1;
  double integral = compute_int_cos(angle1,angle2);
  return norm*integral;
}

template<class Weight = std::function<double (std::pair<double,double>) > >
std::vector<std::pair<double,double> > Fourier_feat(PD D, std::vector<std::pair<double,double> > Z, Weight weight = arctan_weight){
  int m = D.size(); std::vector<std::pair<double,double> > B; int M = Z.size();
  for(int i = 0; i < M; i++){
    double d1 = 0; double d2 = 0; double zx = Z[i].first; double zy = Z[i].second;
    for(int j = 0; j < m; j++){
      double x = D[j].first; double y = D[j].second;
      d1 += weight(D[j])*cos(x*zx + y*zy);
      d2 += weight(D[j])*sin(x*zx + y*zy);
    }
    B.emplace_back(d1,d2);
  }
  return B;
}

std::vector<std::pair<double,double> > random_Fourier(double sigma, int M = 1000){
  std::normal_distribution<double> distrib(0,1); std::vector<std::pair<double,double> > Z; std::random_device rd;
  for(int i = 0; i < M; i++){
    std::mt19937 e1(rd()); std::mt19937 e2(rd());
    double zx = distrib(e1); double zy = distrib(e2);
    Z.emplace_back(zx/sigma,zy/sigma);
  }
  return Z;
}










// ********************************************************************
// Kernel computation.
// ********************************************************************





/** \brief Computes the Linear Persistence Weighted Gaussian Kernel between two persistence diagrams with random Fourier features.
 * \ingroup kernel
 *
 * @param[in] PD1       first persistence diagram.
 * @param[in] PD2       second persistence diagram.
 * @param[in] sigma     bandwidth parameter of the Gaussian Kernel used for the Kernel Mean Embedding of the diagrams.
 * @param[in] weight    weight function for the points in the diagrams.
 * @param[in] M         number of Fourier features (set -1 for exact computation).
 *
 */
template<class Weight = std::function<double (std::pair<double,double>) > >
double linear_persistence_weighted_gaussian_kernel(const PD & PD1, const PD & PD2, double sigma, Weight weight = arctan_weight, int M = 1000){

  if(M == -1){
    int num_pts1 = PD1.size(); int num_pts2 = PD2.size(); double k = 0;
    for(int i = 0; i < num_pts1; i++)
      for(int j = 0; j < num_pts2; j++)
        k += weight(PD1[i])*weight(PD2[j])*exp(-((PD1[i].first-PD2[j].first)*(PD1[i].first-PD2[j].first) + (PD1[i].second-PD2[j].second)*(PD1[i].second-PD2[j].second))/(2*sigma*sigma));
    return k;
  }
  else{
    std::vector<std::pair<double,double> > Z =  random_Fourier(sigma, M);
    std::vector<std::pair<double,double> > B1 = Fourier_feat(PD1,Z,weight);
    std::vector<std::pair<double,double> > B2 = Fourier_feat(PD2,Z,weight);
    double d = 0; for(int i = 0; i < M; i++) d += B1[i].first*B2[i].first + B1[i].second*B2[i].second;
    return d/M;
  }
}

/** \brief Computes the Persistence Scale Space Kernel between two persistence diagrams with random Fourier features.
 * \ingroup kernel
 *
 * @param[in] PD1       first persistence diagram.
 * @param[in] PD2       second persistence diagram.
 * @param[in] sigma     bandwidth parameter of the Gaussian Kernel used for the Kernel Mean Embedding of the diagrams.
 * @param[in] M         number of Fourier features (set -1 for exact computation).
 *
 */
double persistence_scale_space_kernel(const PD & PD1, const PD & PD2, double sigma, int M = 1000){
  PD pd1 = PD1; int numpts = PD1.size();    for(int i = 0; i < numpts; i++)  pd1.emplace_back(PD1[i].second,PD1[i].first);
  PD pd2 = PD2;     numpts = PD2.size();    for(int i = 0; i < numpts; i++)  pd2.emplace_back(PD2[i].second,PD2[i].first);
  return linear_persistence_weighted_gaussian_kernel(pd1, pd2, 2*sqrt(sigma), pss_weight, M) / (2*8*pi*sigma);
}


/** \brief Computes the Gaussian Persistence Weighted Gaussian Kernel between two persistence diagrams with random Fourier features.
 * \ingroup kernel
 *
 * @param[in] PD1       first persistence diagram.
 * @param[in] PD2       second persistence diagram.
 * @param[in] sigma     bandwidth parameter of the Gaussian Kernel used for the Kernel Mean Embedding of the diagrams.
 * @param[in] tau       bandwidth parameter of the Gaussian Kernel used between the embeddings.
 * @param[in] weight    weight function for the points in the diagrams.
 * @param[in] M         number of Fourier features (set -1 for exact computation).
 *
 */
template<class Weight = std::function<double (std::pair<double,double>) > >
double gaussian_persistence_weighted_gaussian_kernel(const PD & PD1, const PD & PD2, double sigma, double tau, Weight weight = arctan_weight, int M = 1000){
  double k1 = linear_persistence_weighted_gaussian_kernel(PD1,PD1,sigma,weight,M);
  double k2 = linear_persistence_weighted_gaussian_kernel(PD2,PD2,sigma,weight,M);
  double k3 = linear_persistence_weighted_gaussian_kernel(PD1,PD2,sigma,weight,M);
  return exp( - (k1+k2-2*k3) / (2*tau*tau)  );
}


/** \brief Computes the Sliced Wasserstein Kernel between two persistence diagrams with sampled directions.
 * \ingroup kernel
 *
 * @param[in] PD1       first persistence diagram.
 * @param[in] PD2       second persistence diagram.
 * @param[in] sigma     bandwidth parameter.
 * @param[in] N         number of points sampled on the circle (set -1 for exact computation).
 *
 */
double sliced_wasserstein_kernel(PD PD1, PD PD2, double sigma, int N = 100){

  if(N == -1){

    // Add projections onto diagonal.
    int n1, n2; n1 = PD1.size(); n2 = PD2.size(); double max_ordinate = std::numeric_limits<double>::lowest();
    for (int i = 0; i < n2; i++){
      max_ordinate = std::max(max_ordinate, PD2[i].second);
      PD1.emplace_back(  (PD2[i].first+PD2[i].second)/2, (PD2[i].first+PD2[i].second)/2  );
    }
    for (int i = 0; i < n1; i++){
      max_ordinate = std::max(max_ordinate, PD1[i].second);
      PD2.emplace_back(  (PD1[i].first+PD1[i].second)/2, (PD1[i].first+PD1[i].second)/2  );
    }
    int num_pts_dgm = PD1.size();

    // Slightly perturb the points so that the PDs are in generic positions.
    int mag = 0; while(max_ordinate > 10){mag++; max_ordinate/=10;}
    double thresh = pow(10,-5+mag);
    srand(time(NULL));
    for (int i = 0; i < num_pts_dgm; i++){
      PD1[i].first += thresh*(1.0-2.0*rand()/RAND_MAX); PD1[i].second += thresh*(1.0-2.0*rand()/RAND_MAX);
      PD2[i].first += thresh*(1.0-2.0*rand()/RAND_MAX); PD2[i].second += thresh*(1.0-2.0*rand()/RAND_MAX);
    }

    // Compute all angles in both PDs.
    std::vector<std::pair<double, std::pair<int,int> > > angles1, angles2;
    for (int i = 0; i < num_pts_dgm; i++){
      for (int j = i+1; j < num_pts_dgm; j++){
        double theta1 = compute_angle(PD1,i,j); double theta2 = compute_angle(PD2,i,j);
        angles1.emplace_back(theta1, std::pair<int,int>(i,j));
        angles2.emplace_back(theta2, std::pair<int,int>(i,j));
      }
    }

    // Sort angles.
    std::sort(angles1.begin(), angles1.end(), sortAngle); std::sort(angles2.begin(), angles2.end(), sortAngle);

    // Initialize orders of the points of both PDs (given by ordinates when theta = -pi/2).
    std::vector<int> orderp1, orderp2;
    for (int i = 0; i < num_pts_dgm; i++){ orderp1.push_back(i); orderp2.push_back(i); }
    std::sort( orderp1.begin(), orderp1.end(), [=](int i, int j){ if(PD1[i].second != PD1[j].second) return (PD1[i].second < PD1[j].second); else return (PD1[i].first > PD1[j].first); } );
    std::sort( orderp2.begin(), orderp2.end(), [=](int i, int j){ if(PD2[i].second != PD2[j].second) return (PD2[i].second < PD2[j].second); else return (PD2[i].first > PD2[j].first); } );

    // Find the inverses of the orders.
    std::vector<int> order1(num_pts_dgm); std::vector<int> order2(num_pts_dgm);
    for(int i = 0; i < num_pts_dgm; i++)  for (int j = 0; j < num_pts_dgm; j++)  if(orderp1[j] == i){  order1[i] = j; break;  }
    for(int i = 0; i < num_pts_dgm; i++)  for (int j = 0; j < num_pts_dgm; j++)  if(orderp2[j] == i){  order2[i] = j; break;  }

    // Record all inversions of points in the orders as theta varies along the positive half-disk.
    std::vector<std::vector<std::pair<int,double> > > anglePerm1(num_pts_dgm);
    std::vector<std::vector<std::pair<int,double> > > anglePerm2(num_pts_dgm);

    int M1 = angles1.size();
    for (int i = 0; i < M1; i++){
      double theta = angles1[i].first; int p = angles1[i].second.first; int q = angles1[i].second.second;
      anglePerm1[order1[p]].emplace_back(p,theta);
      anglePerm1[order1[q]].emplace_back(q,theta);
      int a = order1[p]; int b = order1[q]; order1[p] = b; order1[q] = a;
    }

    int M2 = angles2.size();
    for (int i = 0; i < M2; i++){
      double theta = angles2[i].first; int p = angles2[i].second.first; int q = angles2[i].second.second;
      anglePerm2[order2[p]].emplace_back(p,theta);
      anglePerm2[order2[q]].emplace_back(q,theta);
      int a = order2[p]; int b = order2[q]; order2[p] = b; order2[q] = a;
    }

    for (int i = 0; i < num_pts_dgm; i++){
      anglePerm1[order1[i]].emplace_back(i,pi/2);
      anglePerm2[order2[i]].emplace_back(i,pi/2);
    }

    // Compute the SW distance with the list of inversions.
    double sw = 0;
    for (int i = 0; i < num_pts_dgm; i++){
      std::vector<std::pair<int,double> > U,V; U = anglePerm1[i]; V = anglePerm2[i];
      double theta1, theta2; theta1 = -pi/2;
      unsigned int ku, kv; ku = 0; kv = 0; theta2 = std::min(U[ku].second,V[kv].second);
      while(theta1 != pi/2){
        if(PD1[U[ku].first].first != PD2[V[kv].first].first || PD1[U[ku].first].second != PD2[V[kv].first].second)
          if(theta1 != theta2)
            sw += compute_int(theta1, theta2, U[ku].first, V[kv].first, PD1, PD2);
        theta1 = theta2;
        if (  (theta2 == U[ku].second)  &&  ku < U.size()-1  )  ku++;
        if (  (theta2 == V[kv].second)  &&  kv < V.size()-1  )  kv++;
        theta2 = std::min(U[ku].second, V[kv].second);
      }
    }

    return exp( -(sw/pi)/(2*sigma*sigma) );

  }


  else{
    double step = pi/N; double sw = 0;

    // Add projections onto diagonal.
    int n1, n2; n1 = PD1.size(); n2 = PD2.size();
    for (int i = 0; i < n2; i++)
      PD1.emplace_back(  (PD2[i].first + PD2[i].second)/2,  (PD2[i].first + PD2[i].second)/2  );
    for (int i = 0; i < n1; i++)
      PD2.emplace_back(  (PD1[i].first + PD1[i].second)/2,  (PD1[i].first + PD1[i].second)/2  );
    int n = PD1.size();

    // Sort and compare all projections.
    //#pragma omp parallel for
    for (int i = 0; i < N; i++){
      std::vector<std::pair<int,double> > L1, L2;
      for (int j = 0; j < n; j++){
        L1.emplace_back(   j, PD1[j].first*cos(-pi/2+i*step) + PD1[j].second*sin(-pi/2+i*step)   );
        L2.emplace_back(   j, PD2[j].first*cos(-pi/2+i*step) + PD2[j].second*sin(-pi/2+i*step)   );
      }
      std::sort(L1.begin(),L1.end(), myComp); std::sort(L2.begin(),L2.end(), myComp);
      double f = 0; for (int j = 0; j < n; j++)  f += std::abs(L1[j].second - L2[j].second);
      sw += f*step;
    }
    return exp( -(sw/pi)/(2*sigma*sigma) );
  }
}


} // namespace kernel

} // namespace Gudhi

#endif //KERNEL_H_
