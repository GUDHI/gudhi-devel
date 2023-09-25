/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s): Marc Glisse
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/choose_n_farthest_points.h>

#include <CGAL/Epick_d.h>
#include <CGAL/Random.h>

#include <boost/version.hpp>

#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <chrono>

// ./choose_n_farthest_points [n_pts [seed]]
int main(int argc, char**argv) {
#ifdef HIGH_DIM
  // Less favorable for choose_n_farthest_points_metric
  typedef CGAL::Epick_d<CGAL::Dimension_tag<4>> K;
#else
  typedef CGAL::Epick_d<CGAL::Dimension_tag<2>> K;
#endif
  typedef typename K::Point_d Point_d;

  // Random seed
  CGAL::Random rd((argc >= 3) ? atoi(argv[2]) : 0);

  // Total number of points
  const int n2 = (argc >= 2) ? atoi(argv[1]) : 10000;
  // Size of subsample
  const int n1 = n2;
  const int j = 0;
  std::vector<Point_d> points; points.reserve(n2);
  for (int i = 0; i < n2; ++i) {
#ifdef HIGH_DIM
    points.emplace_back(rd.get_double(-1., 1), rd.get_double(-1., 1), rd.get_double(-1., 1), rd.get_double(-1., 1));
#else
    points.emplace_back(rd.get_double(-1., 1), rd.get_double(-1., 1));
#endif
  }

  K k;
  auto dis = [&](auto&p, auto&q){return sqrt(k.squared_distance_d_object()(p,q));};

  std::vector<Point_d> results; results.reserve(points.size());
  std::vector<K::FT> dists; dists.reserve(points.size());
#ifndef PROFIL
  std::vector<Point_d> results2; results2.reserve(points.size());
  std::vector<K::FT> dists2; dists2.reserve(points.size());
  auto time_start1 = std::chrono::system_clock::now();
  Gudhi::subsampling::choose_n_farthest_points(dis, points, n1,
                                               j,
                                               std::back_inserter(results2),
                                               std::back_inserter(dists2)
                                               );
  auto time_stop1 = std::chrono::system_clock::now();
  auto time_start2 = std::chrono::system_clock::now();
#endif
  Gudhi::subsampling::choose_n_farthest_points_metric(dis, points, n1,
                                               j,
                                               std::back_inserter(results),
                                               std::back_inserter(dists)
                                               );
#ifndef PROFIL
  auto time_stop2 = std::chrono::system_clock::now();
  std::cerr << "Time (in msec.) generic " << std::chrono::duration_cast<std::chrono::milliseconds>((time_stop1 - time_start1)).count()
            << " vs metric "              << std::chrono::duration_cast<std::chrono::milliseconds>((time_stop2 - time_start2)).count()
            << "  (Boost version " << BOOST_VERSION << ")\n";
  if (dists != dists2 || results != results2) {
    // Note that with many points, it often happens that 2 points with the same distance are swapped in the output.
    std::cerr << "Results differ\n";
#ifdef LOG_DIFF
    std::ofstream log_gen("log_gen");
    std::ofstream log_met("log_met");
    for(std::size_t i = 0; i < results.size(); ++i){
      log_gen << dists2[i] << '\t' << results2[i] << '\n';
      log_met << dists [i] << '\t' << results [i] << '\n';
    }
#endif
    return -1;
  }
#endif
  return 0;
}
