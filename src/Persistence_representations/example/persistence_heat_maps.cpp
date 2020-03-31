/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Pawel Dlotko and Mathieu Carriere
 *
 *    Copyright (C) 2019 Inria
 *
 *    Modifications:
 *      - 2018/04 MC: Add persistence heat maps computation
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Persistence_heat_maps.h>
#include <gudhi/common_persistence_representations.h>

#include <iostream>
#include <vector>
#include <utility>
#include <functional>
#include <cmath>

std::function<double(std::pair<double, double>, std::pair<double, double>)> Gaussian_function(double sigma) {
  return [=](std::pair<double, double> p, std::pair<double, double> q) {
    return std::exp(-((p.first - q.first) * (p.first - q.first) + (p.second - q.second) * (p.second - q.second)) /
                    (sigma));
  };
}

using constant_scaling_function = Gudhi::Persistence_representations::constant_scaling_function;
using Persistence_heat_maps = Gudhi::Persistence_representations::Persistence_heat_maps<constant_scaling_function>;

int main(int argc, char** argv) {
  // create two simple vectors with birth--death pairs:

  std::vector<std::pair<double, double> > persistence1;
  std::vector<std::pair<double, double> > persistence2;

  persistence1.push_back(std::make_pair(1, 2));
  persistence1.push_back(std::make_pair(6, 8));
  persistence1.push_back(std::make_pair(0, 4));
  persistence1.push_back(std::make_pair(3, 8));

  persistence2.push_back(std::make_pair(2, 9));
  persistence2.push_back(std::make_pair(1, 6));
  persistence2.push_back(std::make_pair(3, 5));
  persistence2.push_back(std::make_pair(6, 10));

  // over here we define a function we sill put on a top on every birth--death pair in the persistence interval. It can
  // be anything. Over here we will use standard Gaussian
  std::vector<std::vector<double> > filter = Gudhi::Persistence_representations::create_Gaussian_filter(5, 1);

  // creating two heat maps.
  Persistence_heat_maps hm1(persistence1, filter, false, 20, 0, 11);
  Persistence_heat_maps hm2(persistence2, filter, false, 20, 0, 11);

  std::vector<Persistence_heat_maps*> vector_of_maps;
  vector_of_maps.push_back(&hm1);
  vector_of_maps.push_back(&hm2);

  // compute median/mean of a vector of heat maps:
  Persistence_heat_maps mean;
  mean.compute_mean(vector_of_maps);
  Persistence_heat_maps median;
  median.compute_median(vector_of_maps);

  // to compute L^1 distance between hm1 and hm2:
  std::clog << "The L^1 distance is : " << hm1.distance(hm2, 1) << std::endl;

  // to average of hm1 and hm2:
  std::vector<Persistence_heat_maps*> to_average;
  to_average.push_back(&hm1);
  to_average.push_back(&hm2);
  Persistence_heat_maps av;
  av.compute_average(to_average);

  // to compute scalar product of hm1 and hm2:
  std::clog << "Scalar product is : " << hm1.compute_scalar_product(hm2) << std::endl;

  Persistence_heat_maps hm1k(persistence1, Gaussian_function(1.0));
  Persistence_heat_maps hm2k(persistence2, Gaussian_function(1.0));
  Persistence_heat_maps hm1i(persistence1, Gaussian_function(1.0), 20, 20, 0, 11, 0, 11);
  Persistence_heat_maps hm2i(persistence2, Gaussian_function(1.0), 20, 20, 0, 11, 0, 11);
  std::clog << "Scalar product computed with exact 2D kernel on grid is : " << hm1i.compute_scalar_product(hm2i)
            << std::endl;
  std::clog << "Scalar product computed with exact 2D kernel is : " << hm1k.compute_scalar_product(hm2k) << std::endl;

  return 0;
}
