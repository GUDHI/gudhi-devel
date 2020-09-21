/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2019 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "oracle"
#include <boost/test/unit_test.hpp>
#include <gudhi/Unitary_tests_utils.h>

#include <string>

#include <gudhi/Implicit_manifold_intersection_oracle.h>

#include <gudhi/Functions/Function_Sm_in_Rd.h>
#include <gudhi/Functions/Cartesian_product.h>

#include <gudhi/Coxeter_triangulation.h>

#include <random>
#include <cstdlib>

using namespace Gudhi::coxeter_triangulation;

BOOST_AUTO_TEST_CASE(oracle) {

  Function_Sm_in_Rd fun_sph(5.1111, 2);
  auto oracle = make_oracle(fun_sph);
  Coxeter_triangulation<> cox_tr(oracle.amb_d());
  // cox_tr.change_offset(Eigen::VectorXd::Random(oracle.amb_d()));
  
  Eigen::VectorXd seed = fun_sph.seed();
  auto s = cox_tr.locate_point(seed);

  std::size_t num_intersected_edges = 0;
  for (auto f: s.face_range(oracle.cod_d())) {
    auto qr = oracle.intersects(f, cox_tr);
    if (qr.success)
      num_intersected_edges++;
    auto vertex_it = f.vertex_range().begin();
    Eigen::Vector3d p1 = cox_tr.cartesian_coordinates(*vertex_it++);
    Eigen::Vector3d p2 = cox_tr.cartesian_coordinates(*vertex_it++);
    BOOST_CHECK( vertex_it == f.vertex_range().end() );
    Eigen::MatrixXd m(3,3);
    if (qr.success) {      
      m.col(0) = qr.intersection;
      m.col(1) = p1;
      m.col(2) = p2;
      GUDHI_TEST_FLOAT_EQUALITY_CHECK(m.determinant(), 0.0, 1e-10);
    }
  }
  BOOST_CHECK( num_intersected_edges == 3 || num_intersected_edges == 4 );
}
