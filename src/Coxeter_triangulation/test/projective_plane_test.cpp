// #include <gudhi/Debug_utils.h>
// #include <gudhi/IO/output_debug_traces_to_html.h>
#include <iostream>

#include <gudhi/Coxeter_triangulation.h>
#include <gudhi/Functions/Function_RP2_in_R4.h>
#include <gudhi/Implicit_manifold_intersection_oracle.h>
#include <gudhi/Manifold_tracing.h>
#include <gudhi/Cell_complex.h>

#include <gudhi/IO/build_mesh_from_cell_complex.h>
#include <gudhi/IO/output_meshes_to_medit.h>

using namespace Gudhi::coxeter_triangulation;

int main(int argc, char** argv) {

#ifdef GUDHI_DEBUG    
  std::cout << "Debug run\n";
#else
  std::cout << "Release run\n";
#endif

#ifdef DEBUG_TRACES    
  std::cout << "Debug traces are on\n";
#else
  std::cout << "Debug traces are off\n";
#endif

  double radius = 3.1111;
  Function_RP2_in_R4 fun_sph(radius);
  Eigen::VectorXd seed, result;
  fun_sph.seed(seed); 
  fun_sph.evaluate(seed, result);
  std::cout << "result =\n" << result << "\n";
    
  auto oracle = make_oracle(fun_sph);
  double lambda = 0.2;
  if (argc >= 2)
    lambda = atof(argv[1]);
  Coxeter_triangulation<> cox_tr(oracle.amb_d());
  cox_tr.change_offset(Eigen::VectorXd::Random(oracle.amb_d()));
  cox_tr.change_matrix(lambda * cox_tr.matrix());
  
  using MT = Manifold_tracing<Coxeter_triangulation<> >;
  using Out_simplex_map = typename MT::Out_simplex_map;
  std::vector<Eigen::VectorXd> seed_points(1, seed);
  Out_simplex_map interior_simplex_map;
  manifold_tracing_algorithm(seed_points, cox_tr, oracle, interior_simplex_map);
  std::cout << "Interior_simplex_map:\n";
  // for (auto si_pair: interior_simplex_map)
  //   std::cout << "Simplex = \033[1;33m" << si_pair.first << "\033[0m"
  // 	      << " point:\n"  << si_pair.second << "\n";
  std::cout << "\nSize of the initial output = " << interior_simplex_map.size() << "\n\n";
  
  std::size_t intr_d = oracle.amb_d() - oracle.cod_d();
  Cell_complex<Out_simplex_map> cc(intr_d);
  cc.construct_complex(interior_simplex_map);
  
  output_meshes_to_medit(3,
			 "RP2_test",
			 build_mesh_from_cell_complex(cc,
						      Configuration({true, true, true, 1, 5, 3}),
						      Configuration({true, true, true, 2, 13, 14})));
#ifdef GUDHI_COX_OUTPUT_TO_HTML
  write_to_html("RP2_test");
#endif
  return 0;
}
