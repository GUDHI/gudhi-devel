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
#define BOOST_TEST_MODULE "cell_complex"
#include <boost/test/unit_test.hpp>
#include <gudhi/Unitary_tests_utils.h>

#include <gudhi/Debug_utils.h>
#include <gudhi/IO/output_debug_traces_to_html.h>
#include <iostream>

#include <gudhi/Coxeter_triangulation.h>
#include <gudhi/Functions/Function_Sm_in_Rd.h>
#include <gudhi/Functions/Function_torus_in_R3.h>
#include <gudhi/Implicit_manifold_intersection_oracle.h>
#include <gudhi/Manifold_tracing.h>
#include <gudhi/Cell_complex.h>

using namespace Gudhi::coxeter_triangulation;

BOOST_AUTO_TEST_CASE(cell_complex) {
  double radius = 1.1111;
  Function_torus_in_R3 fun_torus(radius, 3 * radius);
  Eigen::VectorXd seed = fun_torus.seed();
  Function_Sm_in_Rd fun_bound(2.5 * radius, 2, seed);

  auto oracle = make_oracle(fun_torus, fun_bound);
  double lambda = 0.2;
  Coxeter_triangulation<> cox_tr(oracle.amb_d());
  cox_tr.change_offset(Eigen::VectorXd::Random(oracle.amb_d()));
  cox_tr.change_matrix(lambda * cox_tr.matrix());

  using MT = Manifold_tracing<Coxeter_triangulation<> >;
  using Out_simplex_map = typename MT::Out_simplex_map;
  std::vector<Eigen::VectorXd> seed_points(1, seed);
  Out_simplex_map interior_simplex_map, boundary_simplex_map;
  manifold_tracing_algorithm(seed_points, cox_tr, oracle, interior_simplex_map, boundary_simplex_map);

  std::size_t intr_d = oracle.amb_d() - oracle.cod_d();
  Cell_complex<Out_simplex_map> cell_complex(intr_d);
  cell_complex.construct_complex(interior_simplex_map, boundary_simplex_map);

  std::size_t interior_sc_map_size0 = cell_complex.interior_simplex_cell_map(0).size();
  std::size_t interior_sc_map_size1 = cell_complex.interior_simplex_cell_map(1).size();
  std::size_t interior_sc_map_size2 = cell_complex.interior_simplex_cell_map(2).size();
  std::size_t boundary_sc_map_size0 = cell_complex.boundary_simplex_cell_map(0).size();
  std::size_t boundary_sc_map_size1 = cell_complex.boundary_simplex_cell_map(1).size();
  BOOST_CHECK(interior_simplex_map.size() == interior_sc_map_size0);
  BOOST_CHECK(boundary_sc_map_size0 - boundary_sc_map_size1 == 0);
  BOOST_CHECK(interior_sc_map_size0 - interior_sc_map_size1 + interior_sc_map_size2 == 0);
}
