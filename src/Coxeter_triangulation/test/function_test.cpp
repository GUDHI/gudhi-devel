/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2019 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

// workaround for the annoying boost message in boost 1.69
#define BOOST_PENDING_INTEGER_LOG2_HPP
#include <boost/integer/integer_log2.hpp>
// end workaround 

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "function"
#include <boost/test/unit_test.hpp>
#include <gudhi/Unitary_tests_utils.h>

#include <gudhi/Functions/Function_Sm_in_Rd.h>
#include <gudhi/Functions/Function_affine_plane_in_Rd.h>
#include <gudhi/Functions/Constant_function.h>
#include <gudhi/Functions/Function_chair_in_R3.h>
#include <gudhi/Functions/Function_torus_in_R3.h>
#include <gudhi/Functions/Function_whitney_umbrella_in_R3.h>
#include <gudhi/Functions/Function_lemniscate_revolution_in_R3.h>
#include <gudhi/Functions/Function_iron_in_R3.h>
#include <gudhi/Functions/Function_moment_curve_in_Rd.h>
#include <gudhi/Functions/random_orthogonal_matrix.h>
#include <gudhi/Functions/Embed_in_Rd.h>
#include <gudhi/Functions/Translate.h>
#include <gudhi/Functions/Linear_transformation.h>
#include <gudhi/Functions/Negation.h>
#include <gudhi/Functions/Cartesian_product.h>
#include <gudhi/Functions/PL_approximation.h>

#include <gudhi/Coxeter_triangulation.h>

#include <string>

#include <random>
#include <cstdlib>

using namespace Gudhi::coxeter_triangulation;

template <class Function>
void test_function(const Function& fun) {
  Eigen::VectorXd seed = fun.seed();
  Eigen::VectorXd res_seed = fun(fun.seed());
  BOOST_CHECK( seed.size() == fun.amb_d() );
  BOOST_CHECK( res_seed.size() == fun.cod_d() );
  for (std::size_t i = 0; i < fun.cod_d(); i++)
    GUDHI_TEST_FLOAT_EQUALITY_CHECK(res_seed(i), 0., 1e-10);
}

BOOST_AUTO_TEST_CASE(function) {

  {
    // the sphere testing part
    std::size_t m = 3, d = 5;
    Eigen::VectorXd center(d); center << 2, 1.5, -0.5, 4.5, -1;
    double radius = 5;
    typedef Function_Sm_in_Rd Function_sphere;
    Function_sphere fun_sphere(radius, m, d, center);
    test_function(fun_sphere);
  }
  {
    // the affine plane testing part
    std::size_t m = 0, d = 5;
    Eigen::MatrixXd normal_matrix = Eigen::MatrixXd::Zero(d, d-m);
    for (std::size_t i = 0; i < d-m; ++i)
      normal_matrix(i,i) = 1;
    typedef Function_affine_plane_in_Rd Function_plane;
    Function_plane fun_plane(normal_matrix);
    test_function(fun_plane);
  }
  {
    // the constant function testing part
    std::size_t k = 2, d = 5;
    auto x = Eigen::VectorXd::Constant(k, 1);
    Constant_function fun_const(d, k, x);
    Eigen::VectorXd res_zero = fun_const(Eigen::VectorXd::Zero(d));
    for (std::size_t i = 0; i < k; ++i)
      GUDHI_TEST_FLOAT_EQUALITY_CHECK(res_zero(i), x(i), 1e-10);
  }
  {
    // the chair function
    Function_chair_in_R3 fun_chair;
    test_function(fun_chair);
  }
  {
    // the torus function
    Function_torus_in_R3 fun_torus;
    test_function(fun_torus);
  }
  {
    // the whitney umbrella function
    Function_whitney_umbrella_in_R3 fun_umbrella;
    test_function(fun_umbrella);
  }
  {
    // the lemniscate revolution function
    Function_lemniscate_revolution_in_R3 fun_lemniscate;
    test_function(fun_lemniscate);
  }
  {
    // the iron function
    Function_iron_in_R3 fun_iron;
    test_function(fun_iron);
  }
  {
    Function_moment_curve_in_Rd fun_moment_curve(3, 5);
    test_function(fun_moment_curve);
  }
  {
    // random orthogonal matrix
    Eigen::MatrixXd matrix = random_orthogonal_matrix(5);
    Eigen::MatrixXd id_matrix = matrix.transpose() * matrix;
    for (std::size_t i = 0; i < 5; ++i)
      for (std::size_t j = 0; j < 5; ++j)
	if (i == j)
	  GUDHI_TEST_FLOAT_EQUALITY_CHECK(id_matrix(i,j), 1.0, 1e-10);
	else
	  GUDHI_TEST_FLOAT_EQUALITY_CHECK(id_matrix(i,j), 0.0, 1e-10);
  }
  {
    // function embedding
    Function_iron_in_R3 fun_iron;
    auto fun_embed = make_embedding(fun_iron, 5);
    test_function(fun_iron);

    // function translation
    Eigen::VectorXd off = Eigen::VectorXd::Random(5);
    auto fun_trans = translate(fun_embed, off);
    test_function(fun_trans);

    // function linear transformation
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Random(5, 5);
    BOOST_CHECK( matrix.determinant() != 0. );
    auto fun_lin = make_linear_transformation(fun_trans, matrix);
    test_function(fun_lin);

    // function negative
    auto fun_neg = negation(fun_lin);
    test_function(fun_neg);

    // function product
    typedef Function_Sm_in_Rd Function_sphere;
    Function_sphere fun_sphere(1, 1);
    auto fun_prod = make_product_function(fun_sphere, fun_sphere, fun_sphere);
    test_function(fun_prod);
    
    // function PL approximation
    Coxeter_triangulation<> cox_tr(6);
    typedef Coxeter_triangulation<>::Vertex_handle Vertex_handle;
    auto fun_pl = make_pl_approximation(fun_prod, cox_tr);
    Vertex_handle v0 = Vertex_handle(cox_tr.dimension(), 0);
    Eigen::VectorXd x0 = cox_tr.cartesian_coordinates(v0);
    Eigen::VectorXd value0 = fun_prod(x0);
    Eigen::VectorXd pl_value0 = fun_pl(x0);
    for (int i = 0; i < fun_pl.cod_d(); i++)
      GUDHI_TEST_FLOAT_EQUALITY_CHECK(value0(i), pl_value0(i), 1e-10);
    Vertex_handle v1 = v0;
    v1[0] += 1;
    Eigen::VectorXd x1 = cox_tr.cartesian_coordinates(v1);
    Eigen::VectorXd value1 = fun_prod(x1);
    Eigen::VectorXd pl_value1 = fun_pl(x1);
    for (int i = 0; i < fun_pl.cod_d(); i++)
      GUDHI_TEST_FLOAT_EQUALITY_CHECK(value1(i), pl_value1(i), 1e-10);
    Eigen::VectorXd pl_value_mid = fun_pl(0.5*x0 + 0.5*x1);
    for (int i = 0; i < fun_pl.cod_d(); i++)
      GUDHI_TEST_FLOAT_EQUALITY_CHECK(0.5*value0(i) + 0.5*value1(i), pl_value_mid(i), 1e-10);
  }
}
