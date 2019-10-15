// workaround for the annoying boost message in boost 1.69
#define BOOST_PENDING_INTEGER_LOG2_HPP
#include <boost/integer/integer_log2.hpp>
// end workaround 

#include <iostream>

#include <gudhi/Coxeter_triangulation.h>
#include <gudhi/Functions/Function_affine_plane_in_Rd.h>
#include <gudhi/Functions/Function_Sm_in_Rd.h>
#include <gudhi/Functions/Cartesian_product.h>
#include <gudhi/Functions/Linear_transformation.h>
#include <gudhi/Implicit_manifold_intersection_oracle.h>
#include <gudhi/Manifold_tracing.h>
#include <gudhi/Cell_complex.h>
#include <gudhi/Functions/random_orthogonal_matrix.h>

#include <gudhi/IO/build_mesh_from_cell_complex.h>
#include <gudhi/IO/output_meshes_to_medit.h>

using namespace Gudhi::coxeter_triangulation;

int main(int argc, char** argv) {
  
  // Creating a circle S1 in R2 of specified radius
  double radius = 1.0;
  Function_Sm_in_Rd fun_circle(radius, 1);

  // Creating a flat torus S1xS1 in R4 from two circle functions
  auto fun_flat_torus = make_product_function(fun_circle, fun_circle);

  // Apply a random rotation in R4
  auto matrix = random_orthogonal_matrix(4);
  auto fun_flat_torus_rotated = make_linear_transformation(fun_flat_torus, matrix);

  // Computing the seed of the function fun_flat_torus
  Eigen::VectorXd seed;
  fun_flat_torus_rotated.seed(seed);    
  
  // Defining a domain function that defines the boundary, which is a hyperplane passing by the origin and orthogonal to x.
  Eigen::MatrixXd normal_matrix = Eigen::MatrixXd::Zero(4, 1);
  for (std::size_t i = 0; i < 4; i++)
    normal_matrix(i,0) = -seed(i);
  Function_affine_plane_in_Rd fun_bound(normal_matrix, -seed/2);
    
  // Defining the intersection oracle
  auto oracle = make_oracle(fun_flat_torus_rotated, fun_bound);

  // Define a Coxeter triangulation scaled by a factor lambda.
  // The triangulation is translated by a random vector to avoid violating the genericity hypothesis.
  double lambda = 0.05;
  Coxeter_triangulation<> cox_tr(oracle.amb_d());
  cox_tr.change_offset(Eigen::VectorXd::Random(oracle.amb_d()));
  cox_tr.change_matrix(lambda * cox_tr.matrix());

  // Manifold tracing algorithm
  using MT = Manifold_tracing<Coxeter_triangulation<> >;
  using Out_simplex_map = typename MT::Out_simplex_map;
  std::vector<Eigen::VectorXd> seed_points(1, seed);
  Out_simplex_map interior_simplex_map, boundary_simplex_map;
  manifold_tracing_algorithm(seed_points, cox_tr, oracle, interior_simplex_map, boundary_simplex_map);

  // Constructing the cell complex
  std::size_t intr_d = oracle.amb_d() - oracle.cod_d();
  Cell_complex<Out_simplex_map> cell_complex(intr_d);
  cell_complex.construct_complex(interior_simplex_map, boundary_simplex_map);

  // Output the cell complex to a file readable by medit
  output_meshes_to_medit(3,
			 "flat_torus_with_boundary",
			 build_mesh_from_cell_complex(cell_complex,
						      Configuration({true, true, true, 1, 5, 3}),
						      Configuration({true, true, true, 2, 13, 14})));

  return 0;
}
