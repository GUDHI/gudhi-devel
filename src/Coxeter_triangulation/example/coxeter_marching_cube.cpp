#define CC_STAR_COMPLETION
// #define CC_A_V_VISITORS

#include <iostream>
#include <vector>
#include <fstream>

#include <gudhi/Simple_coxeter_system_remastered.h>

// #include <gudhi/Points_off_io.h>
// #include <gudhi/Coxeter_system.h>
// #include <gudhi/Coxeter_complex.h>
// #include <gudhi/Coxeter_complex/Off_point_range.h>
// #include <gudhi/Clock.h>
#include <gudhi/Hasse_diagram_persistence.h>

#include <CGAL/Epick_d.h>

// #include "memory_usage.h"
// #include "cxx-prettyprint/prettyprint.hpp"
// #include "output_points_to_medit.h"

using Cell_id = typename Simple_coxeter_system::Alcove_id;

// using K = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
// using FT = K::FT;
// using Point_d = K::Point_d;

using Point_d = std::vector<double>;

using Hasse_cell = Gudhi::Hasse_diagram::Hasse_diagram_cell<int, double, double>;
using Hasse_boundary = std::vector<std::pair<Hasse_cell*, int> >;
using Hasse_diagram = std::map<Cell_id, Hasse_cell*>;
using Vertex_map = std::map<Hasse_cell*, Point_d>;

const unsigned amb_d = 2; // Ambient (domain) dimension
const unsigned cod_d = 1; // Codomain dimension
const Simple_coxeter_system cs('A', amb_d);
Hasse_diagram hd;
Vertex_map vm;

// using Point_vector = std::vector< Point_d >;

// using Coxeter_complex = Gudhi::Coxeter_complex<Point_vector, Coxeter_system>;

void mark(const Cell_id& c_id) {
}

bool is_marked(const Cell_id& c_id) {
  return true;
}

void insert(const Cell_id& c_id) {
}

// The function
Eigen::VectorXd f(Eigen::VectorXd p) {
  double x = p(0), y = p(1);
  Eigen::VectorXd coords(cod_d);
  coords(0) = x*x + y*y - 25;
  return coords;
}

bool intersects(const Cell_id& f_id, Eigen::MatrixXd& li_matrix, Eigen::VectorXd& last_column) {
  std::cout << " Intersects called for face " << f_id << "\n";
  std::size_t curr_row = cod_d;
  for (std::size_t k: cs.normal_basis_range(f_id)) {
    std::size_t j = std::floor(0.5*(1 + std::sqrt(1+8*k)));
    std::size_t i = (j*j + j - 2)/2 - k;
    std::cout << "  (k, i, j) = (" << k << ", " << i << ", " << j << ")\n";
    for (std::size_t col = 0; col < amb_d; ++col)
      li_matrix(curr_row, col) = 0;
    for (std::size_t l = i; l < j; ++l)
      for (std::size_t col = 0; col < amb_d; ++col)
        li_matrix(curr_row, col) += cs.simple_root_matrix()(l, col);
    last_column(curr_row) = f_id[k] / f_id.level();
    curr_row++;
  }
  // for (; curr_row < amb_d + cod_d; ++curr_row) {
  //   for (std::size_t col = 0; col < amb_d; ++col)
  //     li_matrix(curr_row, col) = 0;
  // } 
  std::cout << "  lin_interpolation_matrix (after completion):\n" << li_matrix << "\n";
  std::cout << "  last_column:\n" << last_column << "\n";
  Eigen::VectorXd intersection = li_matrix.colPivHouseholderQr().solve(last_column);
  std::cout << "  intersection:\n" << intersection << "\n";
  std::cout << "  rel_error = " << (li_matrix * intersection - last_column).norm() / last_column.norm() << "\n";
  if ((li_matrix * intersection - last_column).norm() / last_column.norm() > 1e-5 / f_id.level())
    return false;
  Cell_id i_id = cs.query_point_location(intersection, f_id.level());
  std::cout << "  cell containing the intersection: " << i_id << "\n";
  for (std::size_t k = 0; k < i_id.size(); ++k)
    if (!f_id.is_fixed(k) && i_id[k] != f_id[k])
      return false;
  return true;
}

void seed_expansion(const Cell_id& c_id) {
  std::cout << "Simplex: "  << c_id << "\n";
  mark(c_id);
  std::vector<Cell_id> meet_faces;
  Eigen::MatrixXd
    point_matrix(amb_d + 1, amb_d + 1),
    value_matrix(cod_d, amb_d + 1);
  int i = 0;
  for (Cell_id v_id: cs.face_range(c_id, 0)) {
    Eigen::VectorXd cart_coords = cs.cartesian_coordinates_of_vertex(v_id);
    std::cout << " Vertex: "  << v_id << "\n" << cart_coords << "\n";
    for (std::size_t j = 0; j < amb_d; ++j)
      point_matrix(j, i) = cart_coords(j);
    point_matrix(amb_d, i) = 1;
    Eigen::VectorXd val_vector = f(cart_coords);
    for (std::size_t j = 0; j < cod_d; ++j)
      value_matrix(j, i) = val_vector(j);
    ++i;
  }
  std::cout << " point_matrix:\n" << point_matrix << "\n";
  Eigen::MatrixXd point_matrix_inv = point_matrix.colPivHouseholderQr().inverse();
  std::cout << " point_matrix_inv:\n" << point_matrix_inv << "\n";
  std::cout << " value_matrix:\n" << value_matrix << "\n";
  Eigen::MatrixXd lin_interpolation_matrix = value_matrix * point_matrix_inv;
  Eigen::VectorXd last_column = lin_interpolation_matrix.col(amb_d);
  std::cout << " lin_interpolation_matrix (before completion):\n" << lin_interpolation_matrix << "\n";
  lin_interpolation_matrix.conservativeResize(amb_d, amb_d);
  last_column.conservativeResize(amb_d);
  last_column = -last_column;
  for (Cell_id f_id: cs.face_range(c_id, cod_d))
    if (intersects(f_id, lin_interpolation_matrix, last_column)) {
      std::cout << " Result = true\n";
      meet_faces.push_back(f_id);
    }
    else
      std::cout << " Result = false\n";
  for (const Cell_id& f_id: meet_faces)
    if (!is_marked(f_id))
      for (Cell_id cf_id: cs.coface_range(c_id, amb_d))
        if (!is_marked(cf_id))
          seed_expansion(cf_id);
}

int main(int argc, char * const argv[]) {
  double level = 2;
  std::vector<Point_d> seed_points = {{5,0}};
  std::cout << "root_t_:\n" << cs.simple_root_matrix() << "\n";
  for (const Point_d& p: seed_points) {
    seed_expansion(cs.query_point_location(p, level));
  }
  Simple_coxeter_system scs_test('A', 4);
  Cell_id a_id(1, 1);
  // a_id.push_back(0, true);
  // a_id.push_back(0);
  // a_id.push_back(0);
  // a_id.push_back(0, true);
  // a_id.push_back(0);
  // a_id.push_back(0);
  // a_id.push_back(0);
  // a_id.push_back(0);
  // a_id.push_back(1, true);
  // a_id.push_back(1, true);

  // a_id.push_back(0, true);
  // a_id.push_back(0, true);
  // a_id.push_back(0, true);
  // a_id.push_back(0, true);
  // a_id.push_back(0, true);
  // a_id.push_back(0, true);
  // a_id.push_back(0, true);
  // a_id.push_back(0, true);
  // a_id.push_back(0, true);
  // a_id.push_back(0, true);
  
  // std::cout << a_id << " in dimension 4.\n";
  // for (std::size_t i: scs_test.normal_basis_range(a_id)) {
  //   std::cout << i << "\n";
  // }
}
