/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <sstream>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "skeleton_blocker_geometric_complex"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/Skeleton_blocker.h>

struct Geometry_trait {
  typedef std::vector<double> Point;
};

typedef Geometry_trait::Point Point;
typedef Gudhi::skeleton_blocker::Skeleton_blocker_simple_geometric_traits<Geometry_trait> Complex_geometric_traits;
typedef Gudhi::skeleton_blocker::Skeleton_blocker_geometric_complex< Complex_geometric_traits > Complex;
typedef Complex::Vertex_handle Vertex_handle;
typedef Complex_geometric_traits::Root_vertex_handle Root_vertex_handle;

BOOST_AUTO_TEST_CASE(test_skeleton_blocker_off_reader_writer) {
  Complex complex;
  Gudhi::skeleton_blocker::Skeleton_blocker_off_reader<Complex> off_reader("test2.off", complex);
  BOOST_CHECK(off_reader.is_valid());
  
  std::clog << "complex has " <<
      complex.num_vertices() << " vertices, " <<
      complex.num_blockers() << " blockers, " <<
      complex.num_edges() << " edges and " <<
      complex.num_triangles() << " triangles.";

  BOOST_CHECK(complex.num_vertices() == 7);
  BOOST_CHECK(complex.num_edges() == 12);
  BOOST_CHECK(complex.num_triangles() == 6);

  Gudhi::skeleton_blocker::Skeleton_blocker_off_writer<Complex> off_writer("tmp.off", complex);
  Complex same;
  Gudhi::skeleton_blocker::Skeleton_blocker_off_reader<Complex> off_reader2("tmp.off", same);

  std::clog << "\ncomplex:" << complex.to_string() << std::endl;
  std::clog << "\nsame:" << same.to_string() << std::endl;

  BOOST_CHECK(complex == same);
}

BOOST_AUTO_TEST_CASE(test_skeleton_blocker_abstract_link) {
  Complex complex;
  Gudhi::skeleton_blocker::Skeleton_blocker_off_reader<Complex> off_reader("test2.off", complex);
  BOOST_CHECK(off_reader.is_valid());

  std::clog << "complex has " <<
      complex.num_vertices() << " vertices, " <<
      complex.num_blockers() << " blockers, " <<
      complex.num_edges() << " edges and " <<
      complex.num_triangles() << " triangles.";

  BOOST_CHECK(complex.num_vertices() == 7);
  BOOST_CHECK(complex.num_edges() == 12);
  BOOST_CHECK(complex.num_triangles() == 6);

  auto link_0 = complex.abstract_link(Vertex_handle(0));

  std::clog << "\n link(0):" << link_0.to_string() << std::endl;

  BOOST_CHECK(link_0.num_vertices() == 2);
  BOOST_CHECK(link_0.num_edges() == 1);
  BOOST_CHECK(link_0.num_blockers() == 0);
  
  // Check the 2 link vertices
  auto vertex_handle = link_0.vertex_range().begin();
  BOOST_CHECK(link_0[*vertex_handle].get_id() == Root_vertex_handle(1));
  vertex_handle++;
  BOOST_CHECK(link_0[*(vertex_handle)].get_id() == Root_vertex_handle(4));
  
  // Check the lonely link edge
  auto edge_handle = link_0.edge_range().begin();
  BOOST_CHECK(link_0[*edge_handle].first() == Root_vertex_handle(1));
  BOOST_CHECK(link_0[*(edge_handle)].second() == Root_vertex_handle(4));

  auto link_geometric_0 = complex.link(Vertex_handle(0));
  std::clog << "\n link_geometric(0):" << link_geometric_0.to_string() << std::endl;

  BOOST_CHECK(link_0 == link_geometric_0);

  auto print_point = [&](Vertex_handle v) {
    for (auto x : link_geometric_0.point(v)) std::clog << x << " ";
    std::clog << std::endl;
  };

  std::for_each(link_geometric_0.vertex_range().begin(), link_geometric_0.vertex_range().end(), print_point);

  // Check the 2 link vertices
  vertex_handle = link_geometric_0.vertex_range().begin();
  std::vector<double> point_1 = {0,2,0};
  std::vector<double> point_4 = {-1,1,0};
  BOOST_CHECK(link_geometric_0.point(*vertex_handle) == point_1);
  vertex_handle++;
  BOOST_CHECK(link_geometric_0.point(*vertex_handle) == point_4);

}
