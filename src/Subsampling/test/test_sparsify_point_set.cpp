/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clement Jamin
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Subsampling - test sparsify_point_set
#include <boost/test/unit_test.hpp>

#include <gudhi/sparsify_point_set.h>

#include <CGAL/Epick_d.h>
#include <CGAL/Random.h>

#include <vector>
#include <iterator>

BOOST_AUTO_TEST_CASE(test_sparsify_point_set) 
{
  typedef CGAL::Epick_d<CGAL::Dimension_tag<4> >   K;
  typedef typename K::Point_d                      Point_d;
  
  CGAL::Random rd;

  std::vector<Point_d> points;
  for (int i = 0 ; i < 500 ; ++i)
    points.push_back(Point_d(rd.get_double(-1.,1),rd.get_double(-1.,1),rd.get_double(-1.,1),rd.get_double(-1.,1)));

  K k;
  std::vector<Point_d> results;
  Gudhi::subsampling::sparsify_point_set(k, points, 0.5, std::back_inserter(results));
  std::cout << "Before sparsification: " << points.size() << " points.\n";
  std::cout << "After  sparsification: " << results.size() << " points.\n";
  //for (auto p : results)
  //  std::cout << p << "\n";

  BOOST_CHECK(points.size() > results.size());
}
