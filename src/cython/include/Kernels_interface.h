/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Mathieu Carriere
 *
 *    Copyright (C) 2018 INRIA
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

#ifndef INCLUDE_KERNELS_INTERFACE_H_
#define INCLUDE_KERNELS_INTERFACE_H_

#include <gudhi/Sliced_Wasserstein.h>

#include <iostream>
#include <vector>
#include <utility>  // for std::pair

namespace Gudhi {

namespace persistence_diagram {

  double sw(const std::vector<std::pair<double, double>>& diag1,
            const std::vector<std::pair<double, double>>& diag2,
            double sigma, int N) {
    Gudhi::Persistence_representations::Sliced_Wasserstein sw1(diag1, sigma, N);
    Gudhi::Persistence_representations::Sliced_Wasserstein sw2(diag2, sigma, N);
    return sw1.compute_scalar_product(sw2);
  }

  double pwg(const std::vector<std::pair<double, double>>& diag1,
             const std::vector<std::pair<double, double>>& diag2,
             double sigma, int N) {
    Gudhi::Persistence_representations::Persistence_weighted_gaussian pwg1(diag1, sigma, N, Gudhi::Persistence_representations::arctan_weight);
    Gudhi::Persistence_representations::Persistence_weighted_gaussian pwg2(diag2, sigma, N, Gudhi::Persistence_representations::arctan_weight);
    return pwg1.compute_scalar_product(pwg2);
  }

  double pss(const std::vector<std::pair<double, double>>& diag1,
             const std::vector<std::pair<double, double>>& diag2,
             double sigma, int N) {

    std::vector<std::pair<double, double>> pd1 = diag1; int numpts = diag1.size();    for(int i = 0; i < numpts; i++)  pd1.emplace_back(diag1[i].second,diag1[i].first);
    std::vector<std::pair<double, double>> pd2 = diag2;     numpts = diag2.size();    for(int i = 0; i < numpts; i++)  pd2.emplace_back(diag2[i].second,diag2[i].first);

    Gudhi::Persistence_representations::Persistence_weighted_gaussian pwg1(pd1, 2*std::sqrt(sigma), N, Gudhi::Persistence_representations::pss_weight);
    Gudhi::Persistence_representations::Persistence_weighted_gaussian pwg2(pd2, 2*std::sqrt(sigma), N, Gudhi::Persistence_representations::pss_weight);

    return pwg1.compute_scalar_product  (pwg2) / (16*pi*sigma);
  }

  double pss_sym(const std::vector<std::pair<double, double>>& diag1,
             const std::vector<std::pair<double, double>>& diag2,
             double sigma, int N) {

    Gudhi::Persistence_representations::Persistence_weighted_gaussian pwg1(pd1, 2*std::sqrt(sigma), N, Gudhi::Persistence_representations::pss_weight);
    Gudhi::Persistence_representations::Persistence_weighted_gaussian pwg2(pd2, 2*std::sqrt(sigma), N, Gudhi::Persistence_representations::pss_weight);

    return pwg1.compute_scalar_product  (pwg2) / (16*pi*sigma);
  }


  std::vector<std::vector<double> > sw_matrix(const std::vector<std::vector<std::pair<double, double> > >& s1,
                                              const std::vector<std::vector<std::pair<double, double> > >& s2,
                                              double sigma, int N){
    std::vector<std::vector<double> > matrix;
    std::vector<Gudhi::Persistence_representations::Sliced_Wasserstein> ss1;
    int num_diag_1 = s1.size(); for(int i = 0; i < num_diag_1; i++){Gudhi::Persistence_representations::Sliced_Wasserstein sw1(s1[i], sigma, N); ss1.push_back(sw1);}
    std::vector<Gudhi::Persistence_representations::Sliced_Wasserstein> ss2;
    int num_diag_2 = s2.size(); for(int i = 0; i < num_diag_2; i++){Gudhi::Persistence_representations::Sliced_Wasserstein sw2(s2[i], sigma, N); ss2.push_back(sw2);}
    for(int i = 0; i < num_diag_1; i++){
      std::cout << 100.0*i/num_diag_1 << " %" << std::endl;
      std::vector<double> ps; for(int j = 0; j < num_diag_2; j++) ps.push_back(ss1[i].compute_scalar_product(ss2[j])); matrix.push_back(ps);
    }
    return matrix;
  }

  std::vector<std::vector<double> > pwg_matrix(const std::vector<std::vector<std::pair<double, double> > >& s1,
                                               const std::vector<std::vector<std::pair<double, double> > >& s2,
                                               double sigma, int N){
    std::vector<std::vector<double> > matrix; int num_diag_1 = s1.size(); int num_diag_2 = s2.size();
    for(int i = 0; i < num_diag_1; i++){
      std::cout << 100.0*i/num_diag_1 << " %" << std::endl;
      std::vector<double> ps; for(int j = 0; j < num_diag_2; j++) ps.push_back(pwg(s1[i], s2[j], sigma, N)); matrix.push_back(ps);
    }
    return matrix;
  }

  std::vector<std::vector<double> > pss_matrix(const std::vector<std::vector<std::pair<double, double> > >& s1,
                                               const std::vector<std::vector<std::pair<double, double> > >& s2,
                                               double sigma, int N){
    std::vector<std::vector<std::pair<double, double> > > ss1, ss2;
    std::vector<std::vector<double> > matrix; int num_diag_1 = s1.size(); int num_diag_2 = s2.size();
    for(int i = 0; i < num_diag_1; i++){
      std::vector<std::pair<double, double>> pd1 = s1[i]; int numpts = s1[i].size();    
      for(int j = 0; j < numpts; j++)  pd1.emplace_back(s1[i][j].second,s1[i][j].first);
      ss1.push_back(pd1);
    
    for(int i = 0; i < num_diag_2; i++){
      std::vector<std::pair<double, double>> pd2 = s2[i]; int numpts = s2[i].size();    
      for(int j = 0; j < numpts; j++)  pd2.emplace_back(s2[i][j].second,s2[i][j].first);
      ss2.push_back(pd2);

    for(int i = 0; i < num_diag_1; i++){
      std::cout << 100.0*i/num_diag_1 << " %" << std::endl;
      std::vector<double> ps; for(int j = 0; j < num_diag_2; j++) ps.push_back(pss_sym(ss1[i], ss2[j], sigma, N)); matrix.push_back(ps);
    }
    return matrix;
  }

}  // namespace persistence_diagram

}  // namespace Gudhi


#endif  // INCLUDE_KERNELS_INTERFACE_H_
