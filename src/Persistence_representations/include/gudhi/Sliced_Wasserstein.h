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

#ifndef SLICED_WASSERSTEIN_H_
#define SLICED_WASSERSTEIN_H_

// gudhi include
#include <gudhi/read_persistence_from_file.h>
#include <gudhi/common_persistence_representations.h>
#include <gudhi/Debug_utils.h>

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

namespace Gudhi {
namespace Persistence_representations {

/**
 * \class Sliced_Wasserstein gudhi/Sliced_Wasserstein.h
 * \brief A class implementing the Sliced Wasserstein kernel.
 *
 * \ingroup Persistence_representations
 *
 * \details
 * The Sliced Wasserstein kernel is defined as a Gaussian-like kernel between persistence diagrams, where the distance used for
 * comparison is the Sliced Wasserstein distance \f$SW\f$ between persistence diagrams, defined as the integral of the 1-norm
 * between the sorted projections of the diagrams onto all lines passing through the origin:
 *
 * \f$ SW(D_1,D_2)=\int_{\theta\in\mathbb{S}}\,\|\pi_\theta(D_1\cup\pi_\Delta(D_2))-\pi_\theta(D_2\cup\pi_\Delta(D_1))\|_1{\rm d}\theta\f$,
 *
 * where \f$\pi_\theta\f$ is the projection onto the line defined with angle \f$\theta\f$ in the unit circle \f$\mathbb{S}\f$,
 * and \f$\pi_\Delta\f$ is the projection onto the diagonal.
 * The integral can be either computed exactly in \f$O(n^2{\rm log}(n))\f$ time, where \f$n\f$ is the number of points
 * in the diagrams, or approximated by sampling \f$N\f$ lines in the circle in \f$O(Nn{\rm log}(n))\f$ time. The Sliced Wasserstein Kernel is then computed as:
 *
 * \f$ k(D_1,D_2) = {\rm exp}\left(-\frac{SW(D_1,D_2)}{2\sigma^2}\right).\f$
 * 
 * For more details, please see \cite pmlr-v70-carriere17a .
 *
**/

class Sliced_Wasserstein {

 protected:
    Persistence_diagram diagram;
    int approx;
    double sigma;
    std::vector<std::vector<double> > projections, projections_diagonal;

 public:

  void build_rep(){

    if(approx > 0){
      
      double step = pi/this->approx;
      int n = diagram.size();

      for (int i = 0; i < this->approx; i++){
        std::vector<double> l,l_diag;
        for (int j = 0; j < n; j++){

          double px = diagram[j].first; double py = diagram[j].second;
          double proj_diag = (px+py)/2;

          l.push_back         (   px         * cos(-pi/2+i*step) + py         * sin(-pi/2+i*step)   );
          l_diag.push_back    (   proj_diag  * cos(-pi/2+i*step) + proj_diag  * sin(-pi/2+i*step)   );
        }

        std::sort(l.begin(), l.end()); std::sort(l_diag.begin(), l_diag.end());
        projections.push_back(l); projections_diagonal.push_back(l_diag);

      }

    }

  }

  /** \brief Sliced Wasserstein kernel constructor.
   * \ingroup Sliced_Wasserstein
   *
   * @param[in] _diagram  persistence diagram.
   * @param[in] _sigma  bandwidth parameter.
   * @param[in] _approx number of directions used to approximate the integral in the Sliced Wasserstein distance, set to -1 for exact computation.
   *
   */
  Sliced_Wasserstein(const Persistence_diagram & _diagram, double _sigma = 1.0, int _approx = 100){diagram = _diagram; approx = _approx; sigma = _sigma; build_rep();}

  // **********************************
  // Utils.
  // **********************************

  // Compute the angle formed by two points of a PD
  double compute_angle(const Persistence_diagram & diag, int i, int j) const {
    std::pair<double,double> vect; double x1,y1, x2,y2;
    x1 = diag[i].first; y1 = diag[i].second;
    x2 = diag[j].first; y2 = diag[j].second;
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
  double compute_int_cos(double alpha, double beta) const {
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

  double compute_int(double theta1, double theta2, int p, int q, const Persistence_diagram & diag1, const Persistence_diagram & diag2) const {
    double norm = std::sqrt(  (diag1[p].first-diag2[q].first)*(diag1[p].first-diag2[q].first) + (diag1[p].second-diag2[q].second)*(diag1[p].second-diag2[q].second)  );
    double angle1;
    if (diag1[p].first > diag2[q].first)
      angle1 = theta1 - asin( (diag1[p].second-diag2[q].second)/norm  );
    else
      angle1 = theta1 - asin( (diag2[q].second-diag1[p].second)/norm  );
    double angle2 = angle1 + theta2 - theta1;
    double integral = compute_int_cos(angle1,angle2);
    return norm*integral;
  }




  // **********************************
  // Scalar product + distance.
  // **********************************

  /** \brief Evaluation of the Sliced Wasserstein Distance between a pair of diagrams.
   * \ingroup Sliced_Wasserstein
   *
   * @pre       approx attribute needs to be the same for both instances.
   * @param[in] second other instance of class Sliced_Wasserstein. 
   *
   *
   */
  double compute_sliced_wasserstein_distance(const Sliced_Wasserstein & second) const {

    GUDHI_CHECK(this->approx != second.approx, std::invalid_argument("Error: different approx values for representations"));

    Persistence_diagram diagram1 = this->diagram; Persistence_diagram diagram2 = second.diagram; double sw = 0;

    if(this->approx == -1){

      // Add projections onto diagonal.
      int n1, n2; n1 = diagram1.size(); n2 = diagram2.size(); double max_ordinate = std::numeric_limits<double>::lowest();
      for (int i = 0; i < n2; i++){
        max_ordinate = std::max(max_ordinate, diagram2[i].second);
        diagram1.emplace_back(  (diagram2[i].first+diagram2[i].second)/2, (diagram2[i].first+diagram2[i].second)/2  );
      }
      for (int i = 0; i < n1; i++){
        max_ordinate = std::max(max_ordinate, diagram1[i].second);
        diagram2.emplace_back(  (diagram1[i].first+diagram1[i].second)/2, (diagram1[i].first+diagram1[i].second)/2  );
      }
      int num_pts_dgm = diagram1.size();

      // Slightly perturb the points so that the PDs are in generic positions.
      int mag = 0; while(max_ordinate > 10){mag++; max_ordinate/=10;}
      double thresh = pow(10,-5+mag);
      srand(time(NULL));
      for (int i = 0; i < num_pts_dgm; i++){
        diagram1[i].first += thresh*(1.0-2.0*rand()/RAND_MAX); diagram1[i].second += thresh*(1.0-2.0*rand()/RAND_MAX);
        diagram2[i].first += thresh*(1.0-2.0*rand()/RAND_MAX); diagram2[i].second += thresh*(1.0-2.0*rand()/RAND_MAX);
      }

      // Compute all angles in both PDs.
      std::vector<std::pair<double, std::pair<int,int> > > angles1, angles2;
      for (int i = 0; i < num_pts_dgm; i++){
        for (int j = i+1; j < num_pts_dgm; j++){
          double theta1 = compute_angle(diagram1,i,j); double theta2 = compute_angle(diagram2,i,j);
          angles1.emplace_back(theta1, std::pair<int,int>(i,j));
          angles2.emplace_back(theta2, std::pair<int,int>(i,j));
        }
      }

      // Sort angles.
      std::sort(angles1.begin(), angles1.end(), [=](const std::pair<double, std::pair<int,int> >& p1, const std::pair<double, std::pair<int,int> >& p2){return (p1.first < p2.first);});
      std::sort(angles2.begin(), angles2.end(), [=](const std::pair<double, std::pair<int,int> >& p1, const std::pair<double, std::pair<int,int> >& p2){return (p1.first < p2.first);});

      // Initialize orders of the points of both PDs (given by ordinates when theta = -pi/2).
      std::vector<int> orderp1, orderp2;
      for (int i = 0; i < num_pts_dgm; i++){ orderp1.push_back(i); orderp2.push_back(i); }
      std::sort( orderp1.begin(), orderp1.end(), [=](int i, int j){ if(diagram1[i].second != diagram1[j].second) return (diagram1[i].second < diagram1[j].second); else return (diagram1[i].first > diagram1[j].first); } );
      std::sort( orderp2.begin(), orderp2.end(), [=](int i, int j){ if(diagram2[i].second != diagram2[j].second) return (diagram2[i].second < diagram2[j].second); else return (diagram2[i].first > diagram2[j].first); } );

      // Find the inverses of the orders.
      std::vector<int> order1(num_pts_dgm); std::vector<int> order2(num_pts_dgm);
      for(int i = 0; i < num_pts_dgm; i++)  for (int j = 0; j < num_pts_dgm; j++)  if(orderp1[j] == i){  order1[i] = j; break;  }
      for(int i = 0; i < num_pts_dgm; i++)  for (int j = 0; j < num_pts_dgm; j++)  if(orderp2[j] == i){  order2[i] = j; break;  }

      // Record all inversions of points in the orders as theta varies along the positive half-disk.
      std::vector<std::vector<std::pair<int,double> > > anglePerm1(num_pts_dgm);
      std::vector<std::vector<std::pair<int,double> > > anglePerm2(num_pts_dgm);

      int m1 = angles1.size();
      for (int i = 0; i < m1; i++){
        double theta = angles1[i].first; int p = angles1[i].second.first; int q = angles1[i].second.second;
        anglePerm1[order1[p]].emplace_back(p,theta);
        anglePerm1[order1[q]].emplace_back(q,theta);
        int a = order1[p]; int b = order1[q]; order1[p] = b; order1[q] = a;
      }

      int m2 = angles2.size();
      for (int i = 0; i < m2; i++){
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
      for (int i = 0; i < num_pts_dgm; i++){
        std::vector<std::pair<int,double> > u,v; u = anglePerm1[i]; v = anglePerm2[i];
        double theta1, theta2; theta1 = -pi/2;
        unsigned int ku, kv; ku = 0; kv = 0; theta2 = std::min(u[ku].second,v[kv].second);
        while(theta1 != pi/2){
          if(diagram1[u[ku].first].first != diagram2[v[kv].first].first || diagram1[u[ku].first].second != diagram2[v[kv].first].second)
            if(theta1 != theta2)
              sw += compute_int(theta1, theta2, u[ku].first, v[kv].first, diagram1, diagram2);
          theta1 = theta2;
          if (  (theta2 == u[ku].second)  &&  ku < u.size()-1  )  ku++;
          if (  (theta2 == v[kv].second)  &&  kv < v.size()-1  )  kv++;
          theta2 = std::min(u[ku].second, v[kv].second);
        }
      }
    }


    else{

      double step = pi/this->approx;
      for (int i = 0; i < this->approx; i++){

        std::vector<double> v1; std::vector<double> l1 = this->projections[i]; std::vector<double> l1bis = second.projections_diagonal[i]; std::merge(l1.begin(), l1.end(), l1bis.begin(), l1bis.end(), std::back_inserter(v1));
        std::vector<double> v2; std::vector<double> l2 = second.projections[i]; std::vector<double> l2bis = this->projections_diagonal[i]; std::merge(l2.begin(), l2.end(), l2bis.begin(), l2bis.end(), std::back_inserter(v2));
        int n = v1.size(); double f = 0;
        for (int j = 0; j < n; j++)  f += std::abs(v1[j] - v2[j]);
        sw += f*step;

      }
    }

    return sw/pi;
  }

  /** \brief Evaluation of the kernel on a pair of diagrams.
   * \ingroup Sliced_Wasserstein
   *
   * @pre       approx and sigma attributes need to be the same for both instances.
   * @param[in] second other instance of class Sliced_Wasserstein.
   *
   */
  double compute_scalar_product(const Sliced_Wasserstein & second) const {
    GUDHI_CHECK(this->sigma != second.sigma, std::invalid_argument("Error: different sigma values for representations"));
    return std::exp(-compute_sliced_wasserstein_distance(second)/(2*this->sigma*this->sigma));
  }

  /** \brief Evaluation of the distance between images of diagrams in the Hilbert space of the kernel.
   * \ingroup Sliced_Wasserstein
   *
   * @pre       approx and sigma attributes need to be the same for both instances.
   * @param[in] second  other instance of class Sliced_Wasserstein.
   *
   */
  double distance(const Sliced_Wasserstein & second) const {
    GUDHI_CHECK(this->sigma != second.sigma, std::invalid_argument("Error: different sigma values for representations"));
    return std::pow(this->compute_scalar_product(*this) + second.compute_scalar_product(second)-2*this->compute_scalar_product(second),  0.5);
  }




}; // class Sliced_Wasserstein
}  // namespace Persistence_representations
}  // namespace Gudhi

#endif  // SLICED_WASSERSTEIN_H_
