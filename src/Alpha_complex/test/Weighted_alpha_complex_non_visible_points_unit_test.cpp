/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "weighted_alpha_complex_non_visible_points"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <CGAL/Epick_d.h>
#include <CGAL/Epeck_d.h>

#include <vector>

#include <gudhi/Alpha_complex.h>
#include <gudhi/Simplex_tree.h>


using list_of_1d_kernel_variants = boost::mpl::list<CGAL::Epeck_d< CGAL::Dynamic_dimension_tag >,
                                                    CGAL::Epeck_d< CGAL::Dimension_tag<1>>,
                                                    CGAL::Epick_d< CGAL::Dynamic_dimension_tag >,
                                                    CGAL::Epick_d< CGAL::Dimension_tag<1>>
                                                    >;

BOOST_AUTO_TEST_CASE_TEMPLATE(Weighted_alpha_complex_non_visible_points, Kernel, list_of_1d_kernel_variants) {
  // check that for 2 closed weighted 1-d points, one with a high weight to hide the second one with a small weight,
  // that the point with a small weight has the same high filtration value than the edge formed by the 2 points
  using Point_d = typename Kernel::Point_d;
  std::vector<Point_d> points;
  std::vector<double> p1 {0.};
  points.emplace_back(p1.begin(), p1.end());
  // closed enough points
  std::vector<double> p2 {0.1};
  points.emplace_back(p2.begin(), p2.end());
  std::vector<typename Kernel::FT> weights {100., 0.01};

  Gudhi::alpha_complex::Alpha_complex<Kernel, true> alpha_complex(points, weights);
  Gudhi::Simplex_tree<> stree;
  BOOST_CHECK(alpha_complex.create_complex(stree));

  std::clog << "Iterator on weighted alpha complex simplices in the filtration order, with [filtration value]:"
            << std::endl;
  for (auto f_simplex : stree.filtration_simplex_range()) {
    std::clog << "   ( ";
    for (auto vertex : stree.simplex_vertex_range(f_simplex)) {
      std::clog << vertex << " ";
    }
    std::clog << ") -> " << "[" << stree.filtration(f_simplex) << "] " << std::endl;
  }

  BOOST_CHECK(stree.filtration(stree.find({0})) == -100.);
  BOOST_CHECK(stree.filtration(stree.find({1})) == stree.filtration(stree.find({0, 1})));
  BOOST_CHECK(stree.filtration(stree.find({1})) > 100000);
}