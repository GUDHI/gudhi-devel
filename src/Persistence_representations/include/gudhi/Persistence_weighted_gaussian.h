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

// gudhi include
#include <gudhi/read_persistence_from_file.h>
#include <gudhi/common_persistence_representations.h>

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

using PD = std::vector<std::pair<double,double> >;
using Weight = std::function<double (std::pair<double,double>) >;

namespace Gudhi {
namespace Persistence_representations {
/**
 * \class Persistence_weighted_gaussian gudhi/Persistence_weighted_gaussian.h
 * \brief A class implementing the Persistence Weighted Gaussian Kernel and a specific case of it called the Persistence Scale Space Kernel.
 *
 * \ingroup Persistence_representations
 *
 * \details
 * The Persistence Weighted Gaussian Kernel is built with Gaussian Kernel Mean Embedding, meaning that each persistence diagram is first
 * sent to the Hilbert space of a Gaussian kernel with bandwidth parameter \f$\sigma >0\f$ using a weighted mean embedding \f$\Phi\f$:
 *
 * \f$ \Phi\,:\,D\,\rightarrow\,\sum_{p\in D}\,w(p)\,{\rm exp}\left(-\frac{\|p-\cdot\|_2^2}{2\sigma^2}\right) \f$,
 *
 * Usually, the weight function is chosen to be an arctan function of the distance of the point to the diagonal:
 * \f$w(p) = {\rm arctan}(C\,|y-x|^\alpha)\f$, for some parameters \f$C,\alpha >0\f$.
 * Then, their scalar product in this space is computed:
 *
 * \f$ k(D_1,D_2)=\langle\Phi(D_1),\Phi(D_2)\rangle
 * \,=\,\sum_{p\in D_1}\,\sum_{q\in D_2}\,w(p)\,w(q)\,{\rm exp}\left(-\frac{\|p-q\|_2^2}{2\sigma^2}\right).\f$
 *
 * Note that one may apply a second Gaussian kernel to their distance in this space and still get a kernel.
 *
 * It follows that the computation time is \f$O(n^2)\f$ where \f$n\f$ is the number of points
 * in the diagrams. This time can be improved by computing approximations of the kernel
 * with \f$m\f$ Fourier features \cite Rahimi07randomfeatures. In that case, the computation time becomes \f$O(mn)\f$.
 *
 * The Persistence Scale Space Kernel is a Persistence Weighted Gaussian Kernel between modified diagrams:
 * the symmetric of each point with respect to the diagonal is first added in each diagram, and then the weight function
 * is set to be +1 if the point is above the diagonal and -1 otherwise.
 * 
 * For more details, please consult <i>Persistence Weighted Kernel for Topological Data Analysis</i>\cite Kusano_Fukumizu_Hiraoka_PWGK 
 * and <i>A Stable Multi-Scale Kernel for Topological Machine Learning</i>\cite Reininghaus_Huber_ALL_PSSK . 
 * It implements the following concepts: Topological_data_with_distances, Topological_data_with_scalar_product.
 *
**/
class Persistence_weighted_gaussian{

 protected:
    PD diagram;
    Weight weight;
    double sigma;
    int approx;

 public:

  /** \brief Persistence Weighted Gaussian Kernel constructor.
   * \ingroup Persistence_weighted_gaussian
   *
   * @param[in] _diagram       persistence diagram.
   * @param[in] _sigma         bandwidth parameter of the Gaussian Kernel used for the Kernel Mean Embedding of the diagrams.
   * @param[in] _approx        number of random Fourier features in case of approximate computation, set to -1 for exact computation.
   * @param[in] _weight        weight function for the points in the diagrams.
   *
   */
  Persistence_weighted_gaussian(PD _diagram, double _sigma = 1.0, int _approx = 1000, Weight _weight = arctan_weight(1,1)){diagram = _diagram; sigma = _sigma; approx = _approx; weight = _weight;}
  
  PD get_diagram() const {return this->diagram;}
  double get_sigma() const {return this->sigma;}
  int get_approx() const {return this->approx;}
  Weight get_weight() const {return this->weight;}


  // **********************************
  // Utils.
  // **********************************

  /** \brief Specific weight of Persistence Scale Space Kernel.
   * \ingroup Persistence_weighted_gaussian
   *
   * @param[in] p point in 2D.
   *
   */
  static double pss_weight(std::pair<double,double> p)     {if(p.second > p.first)  return 1; else return -1;}
  static double linear_weight(std::pair<double,double> p)  {return std::abs(p.second - p.first);}
  static double const_weight(std::pair<double,double> p)   {return 1;}
  static std::function<double (std::pair<double,double>) > arctan_weight(double C, double power)  {return [=](std::pair<double,double> p){return C * atan(std::pow(std::abs(p.second - p.first), power));};}


  std::vector<std::pair<double,double> > Fourier_feat(PD diag, std::vector<std::pair<double,double> > z, Weight weight = arctan_weight(1,1)){
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

  /** \brief Evaluation of the kernel on a pair of diagrams.
   * \ingroup Persistence_weighted_gaussian
   *
   * @param[in] second other instance of class Persistence_weighted_gaussian. Warning: sigma, approx and weight parameters need to be the same for both instances!!! 
   *
   */
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

  /** \brief Evaluation of the distance between images of diagrams in the Hilbert space of the kernel.
   * \ingroup Persistence_weighted_gaussian
   *
   * @param[in] second other instance of class Persistence_weighted_gaussian. Warning: sigma, approx and weight parameters need to be the same for both instances!!! 
   *
   */
  double distance(Persistence_weighted_gaussian second) {
    if(this->sigma != second.get_sigma() || this->approx != second.get_approx()){
      std::cout << "Error: different representations!" << std::endl; return 0;
    }
    else return std::pow(this->compute_scalar_product(*this) + second.compute_scalar_product(second)-2*this->compute_scalar_product(second), 0.5);
  }


};

}  // namespace Persistence_weighted_gaussian
}  // namespace Gudhi

#endif  // PERSISTENCE_WEIGHTED_GAUSSIAN_H_
