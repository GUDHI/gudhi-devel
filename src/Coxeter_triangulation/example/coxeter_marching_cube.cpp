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
Point_d f(Point_d p) {
  double x = p[0], y = p[1];
  std::vector<double> coords(cod_d);
  coords[0] = x*x + y*y - 25;
  return coords;
}

bool intersects(const Cell_id& f_id, Eigen::MatrixXd li_matrix) {
  std::size_t curr_row = 0;
  Eigen::VectorXd val_vector(amb_d);
  for (; curr_row < cod_d; ++curr_row)
    val_vector(curr_row) = 0;
  for (std::size_t k: cs.normal_basis_range(f_id)) {
    std::size_t j = std::floor(0.5*(1 + std::sqrt(1+8*k)));
    std::size_t i = (j*j + j - 2)/2 - k;      
    for (std::size_t col = 0; i < amb_d; ++i)
      li_matrix(curr_row, col) = 0;
    for (std::size_t l = i; l < j; ++l)
      for (std::size_t col = 0; col < amb_d; ++col)
        li_matrix(curr_row, col) += cs.simple_root_matrix()(l, col);
    val_vector(curr_row) = f_id[k];
    curr_row++;
  }
  Eigen::VectorXd intersection = li_matrix.colPivHouseholderQr().solve(val_vector);
  Cell_id i_id = cs.query_point_location(intersection, f_id.level());
  for (std::size_t k = 0; k < i_id.size(); ++k)
    if (!f_id.is_fixed(k) && i_id[k] != f_id[k])
      return false;
  return true;
}

void seed_expansion(const Cell_id& c_id) {
  mark(c_id);
  std::vector<Cell_id> meet_faces;
  Eigen::MatrixXd
    point_matrix(amb_d+1, amb_d),
    value_matrix(amb_d+1, cod_d);
  int i = 0;
  for (Cell_id v_id: cs.face_range(c_id, 0)) {
    Point_d cart_coords = cs.barycenter(v_id);
    for (unsigned j = 0; j < amb_d; ++j)
      point_matrix(i, j) = cart_coords[j];
    Point_d val_vector = f(cart_coords);
    for (unsigned j = 0; j < cod_d; ++j)
      value_matrix(i, j) = val_vector[j];
    ++i;
  }
  Eigen::MatrixXd lin_interpolation_matrix = point_matrix.colPivHouseholderQr().solve(value_matrix);
  lin_interpolation_matrix.transpose();
  lin_interpolation_matrix.resize(amb_d, Eigen::NoChange);
  for (Cell_id f_id: cs.face_range(c_id, cod_d))
    if (intersects(c_id, lin_interpolation_matrix))
      meet_faces.push_back(f_id);
  for (const Cell_id& f_id: meet_faces)
    if (!is_marked(f_id))
      for (Cell_id cf_id: cs.coface_range(c_id, amb_d))
        if (!is_marked(cf_id))
          seed_expansion(cf_id);
}

int main(int argc, char * const argv[]) {
  double level = 1;
  std::vector<Point_d> seed_points;
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
