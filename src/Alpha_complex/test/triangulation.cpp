#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "triangulation"
#include <boost/test/unit_test.hpp>

#include <boost/range/size.hpp>
#include <boost/range/combine.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include <CGAL/Epeck_d.h>
#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/Spatial_sort_traits_adapter_d.h>
#include <CGAL/property_map.h>  // for CGAL::Identity_property_map
#include <CGAL/version.h>  // for CGAL_VERSION_NR
#include <CGAL/NT_converter.h>

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

  using Internal_vertex_handle = std::ptrdiff_t;
  using Triangulation_full_cell = CGAL::Triangulation_ds_full_cell<void, CGAL::TDS_full_cell_mirror_storage_policy>;
  using TDS = CGAL::Triangulation_data_structure<typename Kernel_4::Dimension,
                                                 CGAL::Triangulation_vertex<Kernel_4, Internal_vertex_handle>,
                                                 Triangulation_full_cell >;
  using Triangulation = CGAL::Delaunay_triangulation<Kernel_4, TDS>;
  
  std::vector<Internal_vertex_handle> indices(boost::counting_iterator<Internal_vertex_handle>(0),
                                              boost::counting_iterator<Internal_vertex_handle>(points.size()));
  using Point_property_map = boost::iterator_property_map<typename std::vector<Point_4>::iterator,
                                                          CGAL::Identity_property_map<Internal_vertex_handle>>;
  using Search_traits_d = CGAL::Spatial_sort_traits_adapter_d<Kernel_4, Point_property_map>;

  CGAL::spatial_sort(indices.begin(), indices.end(), Search_traits_d(std::begin(points)));

  Kernel_4 k;
  Triangulation triangulation(k.point_dimension_d_object()(points[0]));
  
  typename Triangulation::Full_cell_handle hint;
  for (auto index : indices) {
    typename Triangulation::Vertex_handle pos = triangulation.insert(points[index], hint);
    if (pos != nullptr) {
      // Save index value as data to retrieve it after insertion
      pos->data() = index;
      hint = pos->full_cell();
    }
  }
  auto f = k.compute_squared_radius_d_object()(points.begin(), points.end());
  
  std::clog << "compute_squared_radius_d_object = " << f << "\n";

}
