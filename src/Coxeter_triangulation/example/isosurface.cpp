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
// #include <gudhi/Hasse_diagram_persistence.h>
// #include "output_hasse_to_medit.h"

#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <Eigen/SVD>

#include "functions/sphere_S1_in_R2.h"

using Cell_id = Gudhi::Cell_id;
// using Point_d = Eigen::VectorXd;
using Kernel = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
using Point_d = Kernel::Point_d;
using CT = Gudhi::Coxeter_triangulation_ds;

namespace std {
  template<>
  struct hash<Cell_id> {
    typedef Cell_id argument_type;
    typedef std::size_t result_type;
    result_type operator()(const argument_type& c) const noexcept {
      return c.value(0);
    }
  };
}

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
		     const Function& fun) {
  std::size_t amb_d = fun.amb_d();
  std::size_t cod_d = fun.cod_d();
  CT ct(amb_d);
  std::unordered_set<Cell_id> facet_cells;
    
  std::queue<Cell_id> queue;
  for (const Point_d& p: seed_points) {
    Cell_id c = ct.locate_point(p, level);
    for (auto f: ct.face_range(c, cod_d))
      if (max_cells.emplace(f).second)
	queue.emplace(f);
  }
  
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
	  if (intersects(f, fun, ct) && max_cells.emplace(f).second)
	    queue.emplace(f);
    std::cout << "queue.size() = " << queue.size() << "\n";
  }
  std::cout << "#max_cells = " << max_cells.size() << "\n";
  std::cout << "#facet_cells = " << facet_cells.size() << "\n";
}

  

int main(int argc, char * const argv[]) {
  Kernel k;
  std::unordered_set<Cell_id> max_cells; 
  
  // Circle
  {
    double r = 5;
    Function_S1_in_R2 fun(5);
    std::vector<Point_d> seed_points = {Gudhi::construct_point(k, r, 0)};
    double level = 1.5;
    compute_complex(seed_points, level, max_cells, fun);
  }
}
