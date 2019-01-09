#include <iostream>
#include <vector>
#include <list>
#include <queue>
#include <fstream>
#include <cstdlib>
#include <unordered_set>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>

#include <gudhi_patches/CGAL/Epick_d.h>
#include <gudhi/random_point_generators.h> // construct_point
#include <gudhi/Coxeter_triangulation_ds.h>
#include <gudhi/Hasse_diagram_persistence.h>
#include "output_max_cells_to_medit.h"
#include "output_allgowerschmidt.h"
#include "output_transition_graph_to_medit.h"

#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <Eigen/SVD>

#include "functions/sphere_S1_in_R2.h"
#include "functions/sphere_S2_in_R3.h"
#include "functions/chair_in_R3.h"

using Cell_id = Gudhi::Cell_id;
// using Point_d = Eigen::VectorXd;
using Kernel = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
using Point_d = Kernel::Point_d;
using CT = Gudhi::Coxeter_triangulation_ds;


template <class Function>
bool intersects(const Cell_id& c,
		const Function& fun,
		const CT& ct) {
  std::size_t cod_d = fun.cod_d();
  Eigen::MatrixXd matrix(cod_d + 1, cod_d + 1);
  for (std::size_t i = 0; i < cod_d + 1; ++i)
    matrix(0, i) = 1;
  std::size_t j = 0;
  for (auto v: ct.vertex_range(c)) {
    Eigen::VectorXd v_coords = fun(ct.cartesian_coordinates(v));
    for (std::size_t i = 1; i < cod_d + 1; ++i)
      matrix(i, j) = v_coords(i-1);
    j++;
  }
  Eigen::VectorXd z(cod_d + 1);
  z(0) = 1;
  for (std::size_t i = 1; i < cod_d + 1; ++i)
    z(i) = 0;
  Eigen::VectorXd lambda = matrix.colPivHouseholderQr().solve(z);
  for (std::size_t i = 0; i < cod_d + 1; ++i)
    if (lambda(i) < 0 || lambda(i) > 1)
      return false;
  return true;
}
 
template <class Point_range,
	  class Function>
void compute_complex(const Point_range& seed_points,
		     double level,
		     std::unordered_set<Cell_id>& max_cells,
		     const Function& fun,
		     bool output_to_medit = true,
		     std::string file_name_prefix = "reconstruction") {
  std::size_t amb_d = fun.amb_d();
  std::size_t cod_d = fun.cod_d();
  CT ct(amb_d);
  std::unordered_set<Cell_id> facet_cells;
  std::unordered_map<Cell_id, std::size_t> order_map;
    
  std::queue<Cell_id> queue;
  for (const Point_d& p: seed_points) {
    Cell_id c = ct.locate_point(p, level);
    for (auto f: ct.face_range(c, cod_d))
      if (intersects(f, fun, ct) && max_cells.emplace(f).second && order_map.emplace(std::make_pair(f, order_map.size())).second)
	queue.emplace(f);
  }

  // std::size_t snapshot_num = 0;
  while (!queue.empty()) {
    Cell_id s = queue.front();
    queue.pop();
    // Graph_node node = boost::add_vertex(graph);
    // nc_map.emplace(std::make_pair(node, s));
    Gudhi::Coxeter_triangulation_ds::Coface_iterator cof_it(s, ct, cod_d+1);
    Gudhi::Coxeter_triangulation_ds::Coface_iterator cof_end;
    for (; cof_it != cof_end; ++cof_it)
      if (facet_cells.emplace(*cof_it).second)
	for (auto f: ct.face_range(*cof_it, cod_d))
	  if (intersects(f, fun, ct) && max_cells.emplace(f).second && order_map.emplace(std::make_pair(f, order_map.size())).second)
	    queue.emplace(f);
    // std::cout << "queue.size() = " << queue.size() << "\n";
    // output_max_cells_to_medit(max_cells, ct, file_name_prefix + std::to_string(snapshot_num++));
  }
  std::cout << "#max_cells = " << max_cells.size() << "\n";
  std::cout << "#facet_cells = " << facet_cells.size() << "\n";

  std::unordered_set<Cell_id> half_set;
  // for (auto c: max_cells)
  //   if (c.value(0) == 0)
  //     half_set.emplace(c);
  // output_max_cells_to_medit(half_set, ct, order_map, file_name_prefix);
  // output_transition_graph_to_medit(half_set, ct, file_name_prefix + "_tg");
  if (output_to_medit)
    output_max_cells_to_medit(max_cells, ct, order_map, cod_d, file_name_prefix);
  // output_allgowerschmidt_to_medit(max_cells, ct, order_map, file_name_prefix + "_as");
  output_transition_graph_to_medit(max_cells, ct, file_name_prefix + "_tg");
}
  

int main(int argc, char * const argv[]) {
  Kernel k;
  std::unordered_set<Cell_id> max_cells; 
  std::size_t exp_number = 1;

  switch (exp_number) {
  // Circle
  case 0: {
    double r = 5;
    Function_S1_in_R2 fun(r);
    std::vector<Point_d> seed_points = {Gudhi::construct_point(k, r+fun.off_[0], fun.off_[1])};
    double level = 15;
    if (argc > 1)
      level = std::atof(argv[1]);
    compute_complex(seed_points, level, max_cells, fun, true, "circle_reconstruction");
    break;
  }  
  // Sphere
  case 1: {
    double r = 5;
    Function_S2_in_R3 fun(r);
    std::vector<Point_d> seed_points = {Gudhi::construct_point(k, r+fun.off_[0], fun.off_[1], fun.off_[2])};
    double level = 1.513;
    if (argc > 1)
      level = std::atof(argv[1]);
    compute_complex(seed_points, level, max_cells, fun, true, "sphere_reconstruction");
    break;
  }
  // Chair
  case 2: {
    Function_chair_in_R3 fun;
    Eigen::VectorXd seed = fun.seed();
    std::vector<Point_d> seed_points {Gudhi::construct_point(k, seed(0), seed(1), seed(2))};
    std::cout << "fun(fun.seed()) = " << fun(fun.seed()) << "\n";
    double level = 35.13;
    if (argc > 1)
      level = std::atof(argv[1]);
    std::cout << "level = " << level << "\n";
    compute_complex(seed_points, level, max_cells, fun, true, "chair_reconstruction");
    break;
  }
  }
}
