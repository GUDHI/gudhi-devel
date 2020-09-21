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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "manifold_tracing"
#include <boost/test/unit_test.hpp>
#include <gudhi/Unitary_tests_utils.h>

#include <iostream>

#include <gudhi/Coxeter_triangulation.h>
#include <gudhi/Functions/Function_Sm_in_Rd.h>
#include <gudhi/Implicit_manifold_intersection_oracle.h>
#include <gudhi/Manifold_tracing.h>

using namespace Gudhi::coxeter_triangulation;

BOOST_AUTO_TEST_CASE(manifold_tracing) {
  // manifold without boundary
  Function_Sm_in_Rd fun_sph(5.1111, 2);
  auto oracle = make_oracle(fun_sph);
  Coxeter_triangulation<> cox_tr(oracle.amb_d());
  // cox_tr.change_offset(Eigen::VectorXd::Random(oracle.amb_d()));
  
  using MT = Manifold_tracing<Coxeter_triangulation<> >;
  Eigen::VectorXd seed = fun_sph.seed();
  std::vector<Eigen::VectorXd> seed_points(1, seed);
  typename MT::Out_simplex_map out_simplex_map;
  manifold_tracing_algorithm(seed_points, cox_tr, oracle, out_simplex_map);

  for (auto si_pair: out_simplex_map) {
    BOOST_CHECK ( si_pair.first.dimension() == oracle.function().cod_d() );
    BOOST_CHECK ( si_pair.second.size() == (long int)oracle.function().amb_d() );
  }
  std::cout << "out_simplex_map.size() = " << out_simplex_map.size() << "\n";
  BOOST_CHECK ( out_simplex_map.size() == 1140 );
  

  // manifold with boundary
  Function_Sm_in_Rd fun_boundary(3.0, 2, fun_sph.seed());
  auto oracle_with_boundary = make_oracle(fun_sph, fun_boundary);
  typename MT::Out_simplex_map interior_simplex_map, boundary_simplex_map;
  manifold_tracing_algorithm(seed_points, cox_tr, oracle_with_boundary, interior_simplex_map, boundary_simplex_map);
  for (auto si_pair: interior_simplex_map) {
    BOOST_CHECK ( si_pair.first.dimension() == oracle.function().cod_d() );
    BOOST_CHECK ( si_pair.second.size() == (long int)oracle.function().amb_d() );
  }
  std::cout << "interior_simplex_map.size() = " << interior_simplex_map.size() << "\n";
  BOOST_CHECK ( interior_simplex_map.size() == 96 );
  for (auto si_pair: boundary_simplex_map) {
    BOOST_CHECK ( si_pair.first.dimension() == oracle.function().cod_d()+1 );
    BOOST_CHECK ( si_pair.second.size() == (long int)oracle.function().amb_d() );
  }
  std::cout << "boundary_simplex_map.size() = " << boundary_simplex_map.size() << "\n";
  BOOST_CHECK ( boundary_simplex_map.size() == 55 );
  
}
