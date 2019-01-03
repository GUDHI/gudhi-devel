#ifndef SPHERE_S2_IN_R3_H_
#define SPHERE_S2_IN_R3_H_

#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <Eigen/SVD>

struct Function_S1_in_R2 {
  Eigen::VectorXd operator()(const Eigen::VectorXd& p) const {
    double x = p(0), y = p(1);
    Eigen::VectorXd coords(cod_d());
    coords(0) = x*x + y*y - r_*r_;
    return coords;
  }

  std::size_t amb_d() const {return 2;};
  std::size_t cod_d() const {return 1;};
  
  Function_S1_in_R2(double r) : r_(r) {}
  double r_;
};


// void test_circle(double level) {
//   const unsigned amb_d = 2; // Ambient (domain) dimension
//   const unsigned cod_d = 1; // Codomain dimension
//   double r = 5;
//   Eigen::Vector2d point1(r, 0.0);
//   std::vector<Point_d> seed_points = {point1};
//   std::string name = "circle";
//   std::cout << "Test " << test_no++ << ": " << name << ", level = " << level <<  "...\n";

//   Hasse_diagram hd;
//   VP_map vp_map;
//   Gudhi::Clock t;
//   compute_hasse_diagram(seed_points, level, amb_d, cod_d, hd, vp_map, f);
//   t.end();
//   std::vector<unsigned> dimensions(amb_d-cod_d+1, 0);
//   int chi = 0;
//   for (auto cell: hd) {
//     dimensions[cell->get_dimension()]++;
//     chi += 1-2*(cell->get_dimension()%2);
//   }
//   std::cout << "Simplices by dimension: " << dimensions << "\n";
//   std::cout << "Euler characteristic = " << chi << "\n";
//   std::cout << "Reconstruction time: " <<  t.num_seconds() << "s\n";
//   output_hasse_to_medit(hd, vp_map, "marching_cube_output_"+name);
//   std::cout << "Wrote the reconstruction in marching_cube_output_" << name << ".mesh\n";
// }

#endif
