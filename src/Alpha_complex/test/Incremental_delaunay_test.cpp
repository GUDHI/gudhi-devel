/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Marc Glisse
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "incremental_delaunay"
#include <boost/test/unit_test.hpp>
#include <boost/range/adaptor/sliced.hpp>

#include <CGAL/Epick_d.h>
#include <CGAL/Epeck_d.h>

#include <vector>

#include <gudhi/Incremental_delaunay.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/random_point_generators.h>
#include <gudhi/Alpha_complex.h>

using Static_kernel_3 = CGAL::Epick_d< CGAL::Dimension_tag<3> >;
using Dynamic_kernel = CGAL::Epick_d< CGAL::Dynamic_dimension_tag >;
// Epick_d should be all that users need, but still test Epeck_d
using Exact_kernel = CGAL::Epeck_d< CGAL::Dynamic_dimension_tag >;
using list_of_kernel_3_variants = boost::mpl::list<Static_kernel_3, Dynamic_kernel, Exact_kernel>;
using Complex = Gudhi::Simplex_tree<>;

template<class K, class R> Complex check_contains_sub_delaunay(K const&, Complex&cplx, R const& points, int i) {
  Gudhi::alpha_complex::Alpha_complex<K> ac(boost::adaptors::slice(points,0,i));
  Complex delaunay;
  ac.create_complex(delaunay, std::numeric_limits<double>::infinity(), false, true /* no filtration */);
  delaunay.for_each_simplex([&](auto sh, int){
      BOOST_CHECK(cplx.find(delaunay.simplex_vertex_range(sh)) != cplx.null_simplex()); });
  return delaunay;
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Incremental_delaunay_random_points, K, list_of_kernel_3_variants) {
  K k(3);
  auto points = Gudhi::generate_points_in_ball_d<K>(100, 3, 1);
  Complex cplx;
  Gudhi::construct_incremental_delaunay(k, cplx, points);
  BOOST_CHECK(cplx.dimension() == 4); // 3 could theoretically happen, but unlikely

  Complex delaunay;
  for(int i=0;i<=100;++i) delaunay = check_contains_sub_delaunay(k, cplx, points, i);

  cplx.for_each_simplex([&](auto sh, int dim){
      if (dim != 4) return;
      // Cannot use slice because our iterators don't satisfy the concepts
      auto verts = cplx.simplex_vertex_range(sh);
      std::vector<Complex::Vertex_handle> face(std::next(verts.begin()), verts.end());
      BOOST_CHECK(delaunay.find(face) == delaunay.null_simplex());
      // 4-simplices are created when adding the last point breaks the Delaunay simplex of the others
      });
}

BOOST_AUTO_TEST_CASE(Incremental_delaunay_degenerate) {
  typedef Static_kernel_3 K;
  K k(3);
  std::vector<K::Point_d> points;
  Complex cplx;
  Gudhi::construct_incremental_delaunay(k, cplx, points);
  BOOST_CHECK(cplx.is_empty());

  points.emplace_back(0,0,0);
  points.emplace_back(0,0,1);
  points.emplace_back(0,0,3);
  Gudhi::construct_incremental_delaunay(k, cplx, points);
  BOOST_CHECK(cplx.num_simplices() == 5);
  BOOST_CHECK(cplx.dimension() == 1);

  cplx.clear();
  points.emplace_back(0,0,2);
  Gudhi::construct_incremental_delaunay(k, cplx, points);
  Complex cplx2;
  cplx2.insert_simplex_and_subfaces({0,1});
  cplx2.insert_simplex_and_subfaces({1,2,3});
  BOOST_CHECK(cplx == cplx2);

  cplx.clear();
  points.emplace_back(0,1,0);
  Gudhi::construct_incremental_delaunay(k, cplx, points);
  BOOST_CHECK(cplx.num_simplices() == 17);
}
