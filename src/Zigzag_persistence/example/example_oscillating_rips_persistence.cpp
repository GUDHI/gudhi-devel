/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <iostream>
#include <vector>

#include <gudhi/Simplex_tree.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Oscillating_rips_persistence.h>

using Filtration_value = double;
using Point = std::vector<double>;
using Barcode = std::vector<Gudhi::persistence_matrix::Persistence_interval<int, Filtration_value> >;

void print_barcode(const Barcode& bars) {
  std::clog << "Resulting barcode:" << std::endl;
  for (auto& bar : bars) {
    std::clog << bar << std::endl;
  }
}

void print_points(const std::vector<Point>& points) {
  std::clog << "Number of points: " << points.size() << std::endl;
  for (const Point& p : points) {
    std::clog << "(" << p[0] << ", " << p[1] << ")" << std::endl;
  }
  std::clog << std::endl;
}

std::vector<Point> build_point_cloud(unsigned int numberOfPoints, int seed) {
  std::vector<Point> finalPoints;
  std::set<Point> points;
  std::random_device dev;
  std::mt19937 rng(dev());
  if (seed > -1) rng.seed(seed);
  std::uniform_real_distribution<double> dist(0, 10);

  for (unsigned int i = 0; i < numberOfPoints; ++i) {
    auto res = points.insert({dist(rng), dist(rng)});
    while (!res.second) {
      res = points.insert({dist(rng), dist(rng)});
    }
    finalPoints.push_back(*res.first);
  }

  print_points(finalPoints);

  return finalPoints;
}

int main(int argc, char* const argv[]) {
  if (argc != 5 && argc != 6) {
    std::clog << "Usage: ./comp nu mu max_dim numberOfPoints [seed]" << std::endl;
    return 0;
  }

  double nu = std::stod(argv[1]);
  double mu = std::stod(argv[2]);
  int maxDim = std::stoi(argv[3]);
  unsigned int numberOfPoints = std::stoi(argv[4]);
  int seed = -1;

  if (argc == 6) seed = std::stoi(argv[5]);

  std::clog << "nu, mu: " << nu << ", " << mu << std::endl;
  std::clog << "max dimension: " << maxDim << std::endl;
  std::clog << "number of points: " << numberOfPoints << std::endl;
  std::clog << "seed: " << seed << std::endl;

  std::clog << "********** Computing " << numberOfPoints << " random points in a 10 x 10 square" << std::endl;
  std::vector<Point> points = build_point_cloud(numberOfPoints, seed);

  std::clog << "********** Computing oscillating rips filtration and persistence" << std::endl;
  //with default templates and parameters. See documentation for more information.
  Barcode res = Gudhi::zigzag_persistence::compute_oscillating_rips_persistence(points, nu, mu, maxDim);
  print_barcode(res);

  return 0;
}
