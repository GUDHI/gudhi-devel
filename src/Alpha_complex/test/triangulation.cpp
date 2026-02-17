#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "triangulation"
#include <boost/test/unit_test.hpp>

#include <CGAL/Epeck_d.h>

#include <vector>
#include <iostream>

// Use static dimension_tag for the user not to be able to set dimension
typedef CGAL::Epeck_d< CGAL::Dimension_tag<4> > Kernel_4;
typedef Kernel_4::Point_d Point_4;
typedef std::vector<Point_4> Vector_4_Points;

BOOST_AUTO_TEST_CASE(triangulation) {
  // ----------------------------------------------------------------------------
  // Init of a list of points
  // ----------------------------------------------------------------------------
  Vector_4_Points points;
  std::vector<double> coords = { 0.0, 0.0, 0.0, 1.0 };
  points.push_back(Point_4(coords.begin(), coords.end()));
  coords = { 0.0, 0.0, 1.0, 0.0 };
  points.push_back(Point_4(coords.begin(), coords.end()));
  coords = { 0.0, 1.0, 0.0, 0.0 };
  points.push_back(Point_4(coords.begin(), coords.end()));
  coords = { 1.0, 0.0, 0.0, 0.0 };
  points.push_back(Point_4(coords.begin(), coords.end()));

  Kernel_4 k;
  auto f = k.compute_squared_radius_d_object()(points.begin(), points.end());
  std::clog << "compute_squared_radius_d_object = " << f << "\n";

}
