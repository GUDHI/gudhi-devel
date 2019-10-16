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
  Function_Sm_in_Rd fun_sph(5.1111, 2);
  auto oracle = make_oracle(fun_sph, 0.3);
  Coxeter_triangulation<> cox_tr(oracle.amb_d());
  cox_tr.change_offset(Eigen::VectorXd::Random(oracle.amb_d()));
  
  using MT = Manifold_tracing<Coxeter_triangulation<> >;
  Eigen::VectorXd seed = fun_sph.seed();
  std::vector<Eigen::VectorXd> seed_points(1, seed);
  typename MT::Out_simplex_map out_simplex_map;
  manifold_tracing_algorithm(seed_points, cox_tr, oracle, out_simplex_map);

  for (auto si_pair: out_simplex_map) {
    BOOST_CHECK ( si_pair.first.dimension() == oracle.function().cod_d() );
    BOOST_CHECK ( si_pair.second.size() == oracle.function().amb_d() );
  }
  BOOST_CHECK ( out_simplex_map.size() == 1054 );
}
