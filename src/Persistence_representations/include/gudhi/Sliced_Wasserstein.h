/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Mathieu Carriere
 *
 *    Copyright (C) 2018  Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef SLICED_WASSERSTEIN_H_
#define SLICED_WASSERSTEIN_H_

// gudhi include
#include <gudhi/read_persistence_from_file.h>
#include <gudhi/common_persistence_representations.h>
#include <gudhi/Debug_utils.h>

#include <vector>     // for std::vector<>
#include <utility>    // for std::pair<>, std::move
#include <algorithm>  // for std::sort, std::max, std::merge
#include <cmath>      // for std::abs, std::sqrt
#include <stdexcept>  // for std::invalid_argument
#include <random>     // for std::random_device
#include <math.h>     // for nextafter
#include <climits>

namespace Gudhi {
namespace Persistence_representations {

/**
 * \class Sliced_Wasserstein gudhi/Sliced_Wasserstein.h
 * \brief A class implementing the Sliced Wasserstein kernel.
 *
 * \ingroup Persistence_representations
 *
 * \details
 * In this class, we compute infinite-dimensional representations of persistence diagrams by using the
 * Sliced Wasserstein kernel (see \ref sec_persistence_kernels for more details on kernels). We recall that
 * infinite-dimensional representations are defined implicitly, so only scalar products and distances are available for
 * the representations defined in this class.
 * The Sliced Wasserstein kernel is defined as a Gaussian-like kernel between persistence diagrams, where the distance
 * used for comparison is the Sliced Wasserstein distance \f$SW\f$ between persistence diagrams, defined as the
 * integral of the 1-norm between the sorted projections of the diagrams onto all lines passing through the origin:
 *
 * \f$ SW(D_1,D_2)=\int_{\theta\in\mathbb{S}}\,\|\pi_\theta(D_1\cup\pi_\Delta(D_2))-\pi_\theta(D_2\cup\pi_\Delta(D_1))\
 * |_1{\rm d}\theta\f$,
 *
 * where \f$\pi_\theta\f$ is the projection onto the line defined with angle \f$\theta\f$ in the unit circle
 * \f$\mathbb{S}\f$, and \f$\pi_\Delta\f$ is the projection onto the diagonal.
 * Assuming that the diagrams are in general position (i.e. there is no collinear triple), the integral can be computed
 * exactly in \f$O(n^2{\rm log}(n))\f$ time, where \f$n\f$ is the number of points in the diagrams. We provide two
 * approximations of the integral: one in which we slightly perturb the diagram points so that they are in general
 * position, and another in which we approximate the integral by sampling \f$N\f$ lines in the circle in
 * \f$O(Nn{\rm log}(n))\f$ time. The Sliced Wasserstein Kernel is then computed as:
 *
 * \f$ k(D_1,D_2) = {\rm exp}\left(-\frac{SW(D_1,D_2)}{2\sigma^2}\right).\f$
 *
 * The first method is usually much more accurate but also
 * much slower. For more details, please see \cite pmlr-v70-carriere17a .
 *
 **/

class Sliced_Wasserstein {
 protected:
  Persistence_diagram diagram;
  int approx;
  double sigma;
  double genericity_factor;
  std::vector<std::vector<double> > projections, projections_diagonal;

  // **********************************
  // Utils.
  // **********************************

  void build_rep() {
    if (approx > 0) {
      int n = diagram.size();
      double step = pi / this->approx;
      for (int i = 0; i < this->approx; i++) {
        std::vector<double> l, l_diag;
        for (int j = 0; j < n; j++) {
          double px = diagram[j].first;
          double py = diagram[j].second;
          double proj_diag = (px + py) / 2;
          l.push_back(px * cos(-pi / 2 + i * step) + py * sin(-pi / 2 + i * step));
          l_diag.push_back(proj_diag * cos(-pi / 2 + i * step) + proj_diag * sin(-pi / 2 + i * step));
        }
        std::sort(l.begin(), l.end());
        std::sort(l_diag.begin(), l_diag.end());
        projections.push_back(std::move(l));
        projections_diagonal.push_back(std::move(l_diag));
      }
      diagram.clear();
    }
  }

  // Compute the angle formed by two points of a PD
  double compute_angle(const Persistence_diagram& diag, int i, int j) const {
    if (diag[i].second == diag[j].second)
      return pi / 2;
    else
      return atan((diag[j].first - diag[i].first) / (diag[i].second - diag[j].second));
  }

  // Compute the angle formed by two points of different PDs, as well as the sign of the difference of their ordinates
  std::pair<double, int> compute_inter_angle(const Persistence_diagram& diag1, const Persistence_diagram& diag2, int i, int j) const {
    if (diag1[i].second == diag2[j].second)
      return std::pair<double,int>(pi / 2, 1);
    else
      if (diag1[i].second > diag2[j].second)
        return std::pair<double,int>(atan((diag2[j].first - diag1[i].first) / (diag1[i].second - diag2[j].second)), 1);
      else
        return std::pair<double,int>(atan((diag2[j].first - diag1[i].first) / (diag1[i].second - diag2[j].second)), -1);
  }

  // compute the integral of <diff|e_theta> between theta1 and theta2
  double close_form_integral(double theta1, double theta2, double diff_x, double diff_y) const {
    double int_cos = sin(theta2) - sin(theta1);
    double int_sin = cos(theta1) - cos(theta2);
    return diff_x * int_cos + diff_y * int_sin;
  }

  double compute_integral(double theta1, double theta2, std::vector<int> orderp1, std::vector<int> orderp2, std::vector<std::vector<std::pair<double,int> > > angles12, const Persistence_diagram& diag1,
                          const Persistence_diagram& diag2) const {
    double integral = 0;
    for (int i = 0; i < orderp1.size(); i++){

      int dgm1_idx = orderp1[i];
      int dgm2_idx = orderp2[i];
      double diff_x = diag1[dgm1_idx].first - diag2[dgm2_idx].first;
      double diff_y = diag1[dgm1_idx].second - diag2[dgm2_idx].second;
      double angle12 = angles12[dgm1_idx][dgm2_idx].first;
      int sign12 = angles12[dgm1_idx][dgm2_idx].second;
      double small_integral;

      //std::cout << theta1 << " " << theta2 << " " << sign12 << " " << angle12 << std::endl;

      if (sign12 == 1){
        if (angle12 < theta1){ // |<p-q|e_theta>| = <p-q|e_theta> within [theta1, theta2]        
          small_integral = this->close_form_integral(theta1, theta2, diff_x, diff_y);
        }
        else{
          if (theta2 < angle12){ // |<p-q|e_theta>| = -<p-q|e_theta> within [theta1, theta2]
            small_integral = -this->close_form_integral(theta1, theta2, diff_x, diff_y);
          }
          else{ // |<p-q|e_theta>| = -<p-q|e_theta> within [theta1, angle12] and |<p-q|e_theta>| = <p-q|e_theta> within [angle12, theta2]
            //std::cout << "breakdown: " << -this->close_form_integral(theta1, angle12, diff_x, diff_y) << " + " << this->close_form_integral(angle12, theta2, diff_x, diff_y) << std::endl;            
            small_integral = -this->close_form_integral(theta1, angle12, diff_x, diff_y) + this->close_form_integral(angle12, theta2, diff_x, diff_y);
          }
        }
      }
      else{
        if (angle12 < theta1){ // |<p-q|e_theta>| = -<p-q|e_theta> within [theta1, theta2]
          small_integral = -this->close_form_integral(theta1, theta2, diff_x, diff_y);
        }
        else{
          if (theta2 < angle12){ // |<p-q|e_theta>| = <p-q|e_theta> within [theta1, theta2]
            small_integral = this->close_form_integral(theta1, theta2, diff_x, diff_y);
          }
          else{ // |<p-q|e_theta>| = <p-q|e_theta> within [theta1, angle12] and |<p-q|e_theta>| = -<p-q|e_theta> within [angle12, theta2]
            //std::cout << "breakdown: " << this->close_form_integral(theta1, angle12, diff_x, diff_y) << " + " << - this->close_form_integral(angle12, theta2, diff_x, diff_y) << std::endl;            
            small_integral = this->close_form_integral(theta1, angle12, diff_x, diff_y) - this->close_form_integral(angle12, theta2, diff_x, diff_y);
          }
        }
      }

      //std::cout << "small int = " << small_integral << std::endl; 
      integral += small_integral;

    }

    return integral;

  }

  // Evaluation of the Sliced Wasserstein Distance between a pair of diagrams.
  double compute_sliced_wasserstein_distance(const Sliced_Wasserstein& second) const {
    GUDHI_CHECK(this->approx == second.approx,
                std::invalid_argument("Error: different approx values for representations"));

    Persistence_diagram diagram1 = this->diagram;
    Persistence_diagram diagram2 = second.diagram;

    double sw = 0;

    if (this->approx == -1) {

      int n1, n2;
      n1 = diagram1.size();
      n2 = diagram2.size();

      // Put diagrams in generic positions.
      for (int i = 0; i < n1; i++) {diagram1[i].first += (2*i)*genericity_factor*DBL_EPSILON; diagram1[i].second += (2*i+1)*genericity_factor*DBL_EPSILON;}
      for (int i = n1; i < n1+n2; i++) {diagram2[i-n1].first += (2*i)*genericity_factor*DBL_EPSILON; diagram2[i-n1].second += (2*i+1)*genericity_factor*DBL_EPSILON;}
//      std::cout << (2*(n1+n2-1)+1)*10*DBL_EPSILON << std::endl;

      // Add projections onto diagonal.
      for (int i = 0; i < n2; i++) {
        double proj_i = (diagram2[i].first + diagram2[i].second) / 2;
        diagram1.emplace_back(proj_i,proj_i);
      }
      for (int i = 0; i < n1; i++) {
        double proj_i = (diagram1[i].first + diagram1[i].second) / 2;
        diagram2.emplace_back(proj_i, proj_i);
      }
      int num_pts_dgm = diagram1.size();

      // Compute all intra-PD angles with the norms of the corresponding points (for sorting these angles later).
      std::vector<std::pair<double, std::pair<std::pair<int,double>, std::pair<int,double> > > > angles1, angles2;
      for (int i = 0; i < num_pts_dgm; i++) {
        for (int j = i + 1; j < num_pts_dgm; j++) {
          double theta1 = compute_angle(diagram1, i, j);
          double theta2 = compute_angle(diagram2, i, j);
          double p1ix = diagram1[i].first; double p1iy = diagram1[i].second; double p1jx = diagram1[j].first; double p1jy = diagram1[j].second; 
          double p2ix = diagram2[i].first; double p2iy = diagram2[i].second; double p2jx = diagram2[j].first; double p2jy = diagram2[j].second; 
          double n1i = p1ix*p1ix + p1iy*p1iy; double n1j = p1jx*p1jx + p1jy*p1jy; double n2i = p2ix*p2ix + p2iy*p2iy; double n2j = p2jx*p2jx + p2jy*p2jy;
          double m1 = std::min(n1i,n1j); double M1 = std::max(n1i,n1j); double m2 = std::min(n2i,n2j); double M2 = std::max(n2i,n2j);
          std::pair<int,double> iM1(i,M1); std::pair<int,double> jm1(j,m1); std::pair<int,double> iM2(i,M2); std::pair<int,double> jm2(j,m2);
          angles1.emplace_back(theta1, std::pair<std::pair<int,double>, std::pair<int,double>>(iM1, jm1));
          angles2.emplace_back(theta2, std::pair<std::pair<int,double>, std::pair<int,double>>(iM2, jm2));
        }
      }
      int num_angles = angles1.size();

      // Compute all inter-PD angles.
      std::vector<std::vector<std::pair<double,int> > > angles12(num_pts_dgm, std::vector<std::pair<double,int> >(num_pts_dgm));
      for (int i = 0; i < num_pts_dgm; i++) {
        for (int j = 0; j < num_pts_dgm; j++) {
          std::pair<double,int> theta12 = compute_inter_angle(diagram1, diagram2, i, j);
          angles12[i][j] = theta12;
        }
      }

      // Sort angles.
      std::sort(angles1.begin(), angles1.end(),
                [](const std::pair<double, std::pair<std::pair<int,double>, std::pair<int,double>> >& p1,
                   const std::pair<double, std::pair<std::pair<int,double>, std::pair<int,double>> >& p2) {
        if (p1.first == p2.first){
          if (p1.second.first.second == p2.second.first.second)
            return (p1.second.second.second < p2.second.second.second);
          else
            return (p1.second.first.second < p2.second.first.second);
        } else 
          return (p1.first < p2.first); 
      });
      std::sort(angles2.begin(), angles2.end(),
                [](const std::pair<double, std::pair<std::pair<int,double>, std::pair<int,double>> >& p1,
                   const std::pair<double, std::pair<std::pair<int,double>, std::pair<int,double>> >& p2) { 
        if (p1.first == p2.first){
          if (p1.second.first.second == p2.second.first.second)
            return (p1.second.second.second < p2.second.second.second);
          else
            return (p1.second.first.second < p2.second.first.second);
        } else
          return (p1.first < p2.first); 
      });

      // Initialize orders of the points of both PDs (given by ordinates when theta = -pi/2).
      std::vector<int> orderp1, orderp2;
      for (int i = 0; i < num_pts_dgm; i++) {
        orderp1.push_back(i);
        orderp2.push_back(i);
      }
      std::sort(orderp1.begin(), orderp1.end(), [&](int i, int j) {
        if (diagram1[i].second != diagram1[j].second)
          return (diagram1[i].second > diagram1[j].second);
        else
          return (diagram1[i].first < diagram1[j].first);
      });
      std::sort(orderp2.begin(), orderp2.end(), [&](int i, int j) {
        if (diagram2[i].second != diagram2[j].second)
          return (diagram2[i].second > diagram2[j].second);
        else
          return (diagram2[i].first < diagram2[j].first);
      });

      // Find the inverses of the orders.
      std::vector<int> order1(num_pts_dgm);
      std::vector<int> order2(num_pts_dgm);
      for (int i = 0; i < num_pts_dgm; i++) {
        order1[orderp1[i]] = i;
        order2[orderp2[i]] = i;
      }

      double theta1 = -pi/2;
      double theta2;
      int try_idx1 = 0;
      int try_idx2 = 0;
      int order_to_swap, i_to_swap, j_to_swap;

      while( (try_idx1 < num_angles) || (try_idx2 < num_angles)){

        // Find theta2 + next indices to swap
        if ((try_idx1 < num_angles) && (try_idx2 < num_angles)){
          if ((angles1[try_idx1].first < angles2[try_idx2].first) && (try_idx1 < num_angles)){
            theta2 = angles1[try_idx1].first;
            order_to_swap = 1;
            i_to_swap = angles1[try_idx1].second.first.first;
            j_to_swap = angles1[try_idx1].second.second.first;
            try_idx1 += 1;
          }
          else{
            theta2 = angles2[try_idx2].first;
            order_to_swap = 2;
            i_to_swap = angles2[try_idx2].second.first.first;
            j_to_swap = angles2[try_idx2].second.second.first;
            try_idx2 += 1;
          }
        }
        else{
          if (try_idx1 == num_angles){
            theta2 = angles2[try_idx2].first;
            order_to_swap = 2;
            i_to_swap = angles2[try_idx2].second.first.first;
            j_to_swap = angles2[try_idx2].second.second.first;
            try_idx2 += 1;
          }
          else{
            theta2 = angles1[try_idx1].first;
            order_to_swap = 1;
            i_to_swap = angles1[try_idx1].second.first.first;
            j_to_swap = angles1[try_idx1].second.second.first;
            try_idx1 += 1;
          }
        }

        // compute integral
        double integral_value = compute_integral(theta1, theta2, orderp1, orderp2, angles12, diagram1, diagram2);
        sw += integral_value;

//        std::cout << theta1 << " " << theta2 << std::endl;
//        std::cout << "orderp1 = [ ";
//        for (int i = 0; i < orderp1.size(); i++)  std::cout << orderp1[i] << " ";
//        std::cout << "]" << std::endl;
//        std::cout << "orderp2 = [ ";
//        for (int i = 0; i < orderp2.size(); i++)  std::cout << orderp2[i] << " ";
//        std::cout << "]" << std::endl;
//        std::cout << "int = " << integral_value << std::endl;

        // swap order
        if (order_to_swap == 1){
          orderp1[order1[i_to_swap]] = j_to_swap;
          orderp1[order1[j_to_swap]] = i_to_swap;
          int tmp_pos = order1[i_to_swap];
          order1[i_to_swap] = order1[j_to_swap];
          order1[j_to_swap] = tmp_pos;
 
        }
        else{
          orderp2[order2[i_to_swap]] = j_to_swap;
          orderp2[order2[j_to_swap]] = i_to_swap;
          int tmp_pos = order2[i_to_swap];
          order2[i_to_swap] = order2[j_to_swap];
          order2[j_to_swap] = tmp_pos;
        }

        // update theta1
        theta1 = theta2;

      }

      // Last integral between last angle and pi / 2
      theta2 = pi / 2;
      double integral_value = compute_integral(theta1, theta2, orderp1, orderp2, angles12, diagram1, diagram2); 
      sw += integral_value;

//      std::cout << theta1 << " " << theta2 << std::endl;
//      std::cout << "orderp1 = [ ";
//      for (int i = 0; i < orderp1.size(); i++)  std::cout << orderp1[i] << " ";
//      std::cout << "]" << std::endl;
//      std::cout << "orderp2 = [ ";
//      for (int i = 0; i < orderp2.size(); i++)  std::cout << orderp2[i] << " ";
//      std::cout << "]" << std::endl;
//      std::cout << "int = " << integral_value << std::endl;

    } else {
      double step = pi / this->approx;
      std::vector<double> v1, v2;
      for (int i = 0; i < this->approx; i++) {
        v1.clear();
        v2.clear();
        std::merge(this->projections[i].begin(), this->projections[i].end(), second.projections_diagonal[i].begin(),
                   second.projections_diagonal[i].end(), std::back_inserter(v1));
        std::merge(second.projections[i].begin(), second.projections[i].end(), this->projections_diagonal[i].begin(),
                   this->projections_diagonal[i].end(), std::back_inserter(v2));

        int n = v1.size();
        double f = 0;
        for (int j = 0; j < n; j++) f += std::abs(v1[j] - v2[j]);
        sw += f * step;
      }
    }

    return sw / pi;

  }

 public:
  /** \brief Sliced Wasserstein kernel constructor.
   * \implements Topological_data_with_distances, Real_valued_topological_data, Topological_data_with_scalar_product
   * \ingroup Sliced_Wasserstein
   *
   * @param[in] _diagram  persistence diagram.
   * @param[in] _sigma    bandwidth parameter.
   * @param[in] _approx   number of directions used to approximate the integral in the Sliced Wasserstein distance, set
   *                      to -1 for random perturbation. If positive, then projections of the diagram points on all
   *                      directions are stored in memory to reduce computation time.
   *
   */
  Sliced_Wasserstein(const Persistence_diagram& _diagram, double _sigma = 1.0, int _approx = 10, double _genericity_factor=10)
      : diagram(_diagram), approx(_approx), sigma(_sigma), genericity_factor(_genericity_factor) {
    build_rep();
  }

  /** \brief Evaluation of the kernel on a pair of diagrams.
   * \ingroup Sliced_Wasserstein
   *
   * @pre       approx and sigma attributes need to be the same for both instances.
   * @param[in] second other instance of class Sliced_Wasserstein.
   *
   */
  double compute_scalar_product(const Sliced_Wasserstein& second) const {
    GUDHI_CHECK(this->sigma == second.sigma,
                std::invalid_argument("Error: different sigma values for representations"));
    return std::exp(-compute_sliced_wasserstein_distance(second) / (2 * this->sigma * this->sigma));
  }

  /** \brief Evaluation of the distance between images of diagrams in the Hilbert space of the kernel.
   * \ingroup Sliced_Wasserstein
   *
   * @pre       approx and sigma attributes need to be the same for both instances.
   * @param[in] second  other instance of class Sliced_Wasserstein.
   *
   */
  double distance(const Sliced_Wasserstein& second) const {
    GUDHI_CHECK(this->sigma == second.sigma,
                std::invalid_argument("Error: different sigma values for representations"));
    return std::sqrt(this->compute_scalar_product(*this) + second.compute_scalar_product(second) -
                     2 * this->compute_scalar_product(second));
  }

  double sw_distance(const Sliced_Wasserstein& second) const {
    return this->compute_sliced_wasserstein_distance(second);
  }

};  // class Sliced_Wasserstein
}  // namespace Persistence_representations
}  // namespace Gudhi

#endif  // SLICED_WASSERSTEIN_H_
