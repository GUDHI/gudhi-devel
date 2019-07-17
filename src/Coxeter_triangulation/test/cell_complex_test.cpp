#include <iostream>

#include <gudhi/Coxeter_triangulation.h>
#include <gudhi/Functions/Function_Sm_in_Rd.h>
#include <gudhi/Implicit_manifold_intersection_oracle.h>
#include <gudhi/Manifold_tracing.h>
#include <gudhi/Cell_complex.h>

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

int main() {
  Function_Sm_in_Rd fun_sph(1.1111, 2);
  auto oracle = make_oracle(fun_sph, 0.0);
  Coxeter_triangulation<> cox_tr(oracle.amb_d());
  cox_tr.change_offset(Eigen::VectorXd::Random(oracle.amb_d()));
  
  using MT = Manifold_tracing<Coxeter_triangulation<> >;
  using Out_simplex_map = typename MT::Out_simplex_map;
  std::vector<Eigen::VectorXd> seed_points(1, fun_sph.seed());
  Out_simplex_map out_simplex_map;
  manifold_tracing_algorithm(seed_points, cox_tr, oracle, out_simplex_map);

  std::size_t intr_d = oracle.amb_d() - oracle.cod_d();
  Cell_complex<Out_simplex_map> cc(intr_d);
  cc.construct_complex(out_simplex_map);
  auto& sc_map = cc.simplex_cell_map(0);
  for (auto si_pair: out_simplex_map)
    if (sc_map.find(si_pair.first) != sc_map.end())
      std::cout << "Simplex = \033[1;31m" << si_pair.first << "\033[0m\n";
    else
      std::cout << "Simplex = \033[1;34m" << si_pair.first << "\033[0m\n";      
  std::cout << "\nSize of the initial output = " << out_simplex_map.size() << "\n";

  std::cout << "The state of the map:\n";
  for (auto& sc_pair: sc_map)
    std::cout << "Simplex = \033[1;31m" << sc_pair.first << "\033[0m\n";
  std::cout << "Size of the final output = " << sc_map.size() << "\n";

  using Simplex_handle = Coxeter_triangulation<>::Simplex_handle;
  using Vertex = typename Simplex_handle::Vertex;
  using Partition = typename Simplex_handle::OrderedSetPartition;
  using Part = typename Partition::value_type;
  using SC = Simplex_comparator<Simplex_handle>;
  Vertex v1 = {1,1,0};
  Partition omega1 = {Part({0,1,2,3})};
  Simplex_handle s1(v1, omega1);
  Vertex v2 = {1,1,0};
  Partition omega2 = {Part({2}), Part({0,1,3})};
  Simplex_handle s2(v2, omega2);
  
  std::cout << "Vertex comparison: " << v1 << " "
	    << (v1 < v2? "<": "") << (v1 > v2? ">": "") << " " << v2 << "\n";
  std::cout << "Partition comparison: " << omega1 << " "
	    << (omega1 < omega2? "<": "") << (omega1 > omega2? ">": "") << " " << omega2 << "\n";
  std::cout << "Simplex comparison: " << s1 << " "
	    << (SC()(s1, s2)? "<": "") << (SC()(s2, s1)? ">": "") << " " << s2 << "\n";
  return 0;
}
