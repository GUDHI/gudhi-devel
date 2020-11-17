#include <iostream>

#include <gudhi/Coxeter_triangulation.h>
#include <gudhi/Implicit_manifold_intersection_oracle.h>  // for Gudhi::coxeter_triangulation::make_oracle
#include <gudhi/Manifold_tracing.h>
#include <gudhi/Cell_complex.h>
#include <gudhi/Functions/Function_Sm_in_Rd.h>

using namespace Gudhi::coxeter_triangulation;

int main(int argc, char** argv) {
  // Oracle is a circle of radius 1
  double radius = 1.;
  auto oracle = make_oracle(Function_Sm_in_Rd(radius, 1)); 

  // Define a Coxeter triangulation.
  Coxeter_triangulation<> cox_tr(oracle.amb_d());
  // Theory forbids that a vertex of the triangulation lies exactly on the circle.
  // Add some offset to avoid algorithm degeneracies.
  cox_tr.change_offset(-Eigen::VectorXd::Random(oracle.amb_d()));
  // For a better manifold approximation, one can change the circle radius value or change the linear transformation
  // matrix.
  // The number of points and edges will increase with a better resolution.
  //cox_tr.change_matrix(0.5 * cox_tr.matrix());

  // Manifold tracing algorithm
  using Out_simplex_map = typename Manifold_tracing<Coxeter_triangulation<> >::Out_simplex_map;

  std::vector<Eigen::VectorXd> seed_points(1, oracle.seed());
  Out_simplex_map interior_simplex_map;
  manifold_tracing_algorithm(seed_points, cox_tr, oracle, interior_simplex_map);

  // Constructing the cell complex
  std::size_t intr_d = oracle.amb_d() - oracle.cod_d();
  Cell_complex<Out_simplex_map> cell_complex(intr_d);
  cell_complex.construct_complex(interior_simplex_map);

  // List of Hasse_cell pointers to retrieve vertices values from edges
  std::map<Cell_complex<Out_simplex_map>::Hasse_cell*, std::size_t> vi_map;
  std::size_t index = 0;

  std::clog << "Vertices:" << std::endl;
  for (const auto& cp_pair : cell_complex.cell_point_map()) {
    std::clog << index << " : (" << cp_pair.second(0) << ", " << cp_pair.second(1) << ")" << std::endl;
    vi_map.emplace(std::make_pair(cp_pair.first, index++));
  }

  std::clog << "Edges:" << std::endl;
  for (const auto& sc_pair : cell_complex.interior_simplex_cell_map(1)) {
    Cell_complex<Out_simplex_map>::Hasse_cell* edge_cell = sc_pair.second;
    for (const auto& vi_pair : edge_cell->get_boundary()) std::clog << vi_map[vi_pair.first] << " ";
    std::clog << std::endl;
  }
  return 0;
}
