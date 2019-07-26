#include <iostream>

#include <gudhi/Coxeter_triangulation.h>
#include <gudhi/Functions/Function_Sm_in_Rd.h>
#include <gudhi/Functions/Function_torus_in_R3.h>
#include <gudhi/Implicit_manifold_intersection_oracle.h>
#include <gudhi/Manifold_tracing.h>
#include <gudhi/Cell_complex.h>

#include <gudhi/IO/build_mesh_from_cell_complex.h>
#include <gudhi/IO/output_meshes_to_medit.h>

using namespace Gudhi::coxeter_triangulation;

template <class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vector) {
    os << "(";  
  if (vector.empty()) {
    std::cout << ")";
    return os;
  }
  auto v_it = vector.begin();
  os << *v_it++;
  for (; v_it != vector.end(); ++v_it)
    os << ", " << *v_it;
  os << ")";
  return os;
}

int main(int argc, char** argv) {
  double radius = 1.1111;
  Function_torus_in_R3 fun_sph(radius, 3*radius);
  Eigen::VectorXd seed;
  fun_sph.seed(seed);
  Function_Sm_in_Rd fun_bound(2.5*radius, 2, seed);
    
  double thr = 0.1;
  if (argc >= 2)
    thr = atof(argv[1]);
  auto oracle = make_oracle(fun_sph, fun_bound, thr);
  double lambda = 0.2;
  if (argc >= 3)
    lambda = atof(argv[2]);
  Coxeter_triangulation<> cox_tr(oracle.amb_d());
  cox_tr.change_offset(Eigen::VectorXd::Random(oracle.amb_d()));
  cox_tr.change_matrix(lambda * cox_tr.matrix());
  
  using MT = Manifold_tracing<Coxeter_triangulation<> >;
  using Out_simplex_map = typename MT::Out_simplex_map;
  std::vector<Eigen::VectorXd> seed_points(1, seed);
  Out_simplex_map interior_simplex_map, boundary_simplex_map;
  manifold_tracing_algorithm(seed_points, cox_tr, oracle, interior_simplex_map, boundary_simplex_map);
  // std::cout << "Interior_simplex_map:\n";
  // for (auto si_pair: interior_simplex_map)
  //   std::cout << "Simplex = \033[1;33m" << si_pair.first << "\033[0m"
  // 	      << " point:\n"  << si_pair.second << "\n";
  // std::cout << "\nSize of the initial output = " << interior_simplex_map.size() << "\n\n";
  // std::cout << "Boundary_simplex_map:\n";
  // for (auto si_pair: boundary_simplex_map)
  //   std::cout << "Simplex = \033[1;32m" << si_pair.first << "\033[0m"
  // 	      << " point:\n"  << si_pair.second << "\n";
  // std::cout << "\nSize of the initial output = " << boundary_simplex_map.size() << "\n\n";
  

  
  std::size_t intr_d = oracle.amb_d() - oracle.cod_d();
  Cell_complex<Out_simplex_map> cc(intr_d);
  cc.construct_complex(interior_simplex_map, boundary_simplex_map);
  // auto& sc_map = cc.simplex_cell_map(0);
  // for (auto si_pair: out_simplex_map)
  //   if (sc_map.find(si_pair.first) != sc_map.end())
  //     std::cout << "Simplex = \033[1;31m" << si_pair.first << "\033[0m"
  // 		<< " point:\n"  << si_pair.second << "\n";
  //   else
  //     std::cout << "Simplex = \033[1;34m" << si_pair.first << "\033[0m"
  // 		<< " point:\n"  << si_pair.second << "\n";     
  // std::cout << "\nSize of the initial output = " << out_simplex_map.size() << "\n";

  // std::cout << "The state of the map:\n";
  // for (auto& sc_pair: sc_map)
  //   std::cout << "Simplex = \033[1;31m" << sc_pair.first << "\033[0m"
  // 	      << " point:\n"  << cc.cell_point_map().at(sc_pair.second) << "\n";
  // std::cout << "Size of the final output = " << sc_map.size() << "\n";

  // using Simplex_handle = Coxeter_triangulation<>::Simplex_handle;
  // using Vertex = typename Simplex_handle::Vertex;
  // using Partition = typename Simplex_handle::OrderedSetPartition;
  // using Part = typename Partition::value_type;
  // // using SC = Simplex_comparator<Simplex_handle>;
  // Vertex v1 = {1,1,0};
  // Partition omega1 = {Part({0,1,2,3})};
  // Simplex_handle s1(v1, omega1);
  // Vertex v2 = {1,1,0};
  // Partition omega2 = {Part({2}), Part({0,1,3})};
  // Simplex_handle s2(v2, omega2);
  
  // std::cout << "Vertex comparison: " << v1 << " "
  // 	    << (v1 < v2? "<": "") << (v1 > v2? ">": "") << " " << v2 << "\n";
  // std::cout << "Partition comparison: " << omega1 << " "
  // 	    << (omega1 < omega2? "<": "") << (omega1 > omega2? ">": "") << " " << omega2 << "\n";
  // std::cout << "Simplex comparison: " << s1 << " "
  // 	    << (SC()(s1, s2)? "<": "") << (SC()(s2, s1)? ">": "") << " " << s2 << "\n";

  output_meshes_to_medit(3,
			 "test",
			 build_mesh_from_cell_complex(cc,
						      Configuration({true, true, true, 1, 5, 3}),
						      Configuration({true, true, true, 2, 13, 14})));
  return 0;
}
