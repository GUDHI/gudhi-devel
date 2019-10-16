#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "permutahedral_representation"
#include <boost/test/unit_test.hpp>

#include <gudhi/Unitary_tests_utils.h>
#include <gudhi/Freudenthal_triangulation.h>
#include <gudhi/Coxeter_triangulation.h>

BOOST_AUTO_TEST_CASE(permutahedral_representation) {
  
  // Point location check
  typedef std::vector<double> Point;
  typedef Gudhi::coxeter_triangulation::Freudenthal_triangulation<> FK_triangulation;
  typedef typename FK_triangulation::Simplex_handle Simplex_handle;
  typedef typename FK_triangulation::Vertex_handle Vertex_handle;
  typedef typename Simplex_handle::OrderedSetPartition Ordered_set_partition; 
  typedef typename Ordered_set_partition::value_type Part;

  FK_triangulation tr(3);

  // Point location check
  {
    Point point({3, -1, 0});
    Simplex_handle s = tr.locate_point(point);
    BOOST_CHECK( s.vertex() == Vertex_handle({3, -1, 0}) );
    BOOST_CHECK( s.partition() == Ordered_set_partition({Part({0, 1, 2, 3})}) );
  }

  {
    Point point({3.5, -1.5, 0.5});
    Simplex_handle s = tr.locate_point(point);
    BOOST_CHECK( s.vertex() == Vertex_handle({3, -2, 0}) );
    BOOST_CHECK( s.partition() == Ordered_set_partition({Part({0, 1, 2}), Part({3})}) );
  }

  {
    Point point({3.5, -1.8, 0.5});
    Simplex_handle s = tr.locate_point(point);
    BOOST_CHECK( s.vertex() == Vertex_handle({3, -2, 0}) );
    BOOST_CHECK( s.partition() == Ordered_set_partition({Part({0, 2}), Part({1}), Part({3})}) );
  }

  {
    Point point({3.5, -1.8, 0.3});
    Simplex_handle s = tr.locate_point(point);
    BOOST_CHECK( s.vertex() == Vertex_handle({3, -2, 0}) );
    BOOST_CHECK( s.partition() == Ordered_set_partition({Part({0}), Part({2}), Part({1}), Part({3})}) );
  }

  // Dimension check
  BOOST_CHECK( tr.dimension() == 3 );
  // Matrix check
  Eigen::MatrixXd default_matrix = Eigen::MatrixXd::Identity(3, 3);  
  BOOST_CHECK( tr.matrix() == default_matrix );
  // Vector check
  Eigen::MatrixXd default_offset = Eigen::VectorXd::Zero(3);  
  BOOST_CHECK( tr.offset() == default_offset );

  // Barycenter check
  Point point({3.5, -1.8, 0.3});
  Simplex_handle s = tr.locate_point(point);
  Eigen::Vector3d barycenter_cart = Eigen::Vector3d::Zero();
  for (auto v: s.vertex_range())
    for (std::size_t i = 0; i < v.size(); i++)
      barycenter_cart(i) += v[i];
  barycenter_cart /= 4.; // simplex is three-dimensional
  Eigen::Vector3d barycenter = tr.barycenter(s);
  for (std::size_t i = 0; i < barycenter.size(); i++)
    GUDHI_TEST_FLOAT_EQUALITY_CHECK(barycenter(i), barycenter_cart(i), 1e-7);

  // Barycenter check for twice the scale
  s = tr.locate_point(point, 2);
  barycenter_cart = Eigen::Vector3d::Zero();
  for (auto v: s.vertex_range())
    for (std::size_t i = 0; i < v.size(); i++)
      barycenter_cart(i) += v[i];
  barycenter_cart /= 3.; // simplex is now a two-dimensional face
  barycenter_cart /= 2.; // scale
  barycenter = tr.barycenter(s, 2);
  for (std::size_t i = 0; i < barycenter.size(); i++)
    GUDHI_TEST_FLOAT_EQUALITY_CHECK(barycenter(i), barycenter_cart(i), 1e-7);
  
  // Matrix and offset change check
  Eigen::MatrixXd new_matrix(3,3);
  new_matrix << 1, 0, 0, -1, 1, 0, -1, 0, 1;
  Eigen::Vector3d new_offset(1.5, 1, 0.5);
  tr.change_matrix(new_matrix);
  tr.change_offset(new_offset);
  
  BOOST_CHECK( tr.matrix() == new_matrix );
  BOOST_CHECK( tr.offset() == new_offset );
}
