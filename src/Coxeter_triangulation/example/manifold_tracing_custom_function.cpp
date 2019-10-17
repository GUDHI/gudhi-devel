#include <iostream>

#include <gudhi/Coxeter_triangulation.h>
#include <gudhi/Functions/Function_Sm_in_Rd.h>
#include <gudhi/Implicit_manifold_intersection_oracle.h>
#include <gudhi/Manifold_tracing.h>
#include <gudhi/Cell_complex.h>
#include <gudhi/Functions/random_orthogonal_matrix.h>
#include <gudhi/Functions/Linear_transformation.h>

#include <gudhi/IO/build_mesh_from_cell_complex.h>
#include <gudhi/IO/output_meshes_to_medit.h>

using namespace Gudhi::coxeter_triangulation;

/* A definition of a function that defines a 2d surface embedded in R^4, but that normally
 * lives on a complex projective plane.
 * In terms of harmonic coordinates [x:y:z] of points on the complex projective plane,
 * the equation of the manifold is x^3*y + y^3*z + z^3*x = 0.
 * The embedding consists of restricting the manifold to the affine subspace z = 1.
 */
struct Function_surface_on_CP2_in_R4 : public Function {

  Eigen::VectorXd operator()(const Eigen::VectorXd& p) const {
    // The real and imaginary parts of the variables x and y
    double xr = p(0), xi = p(1), yr = p(2), yi = p(3);
    Eigen::VectorXd result(cod_d());
    
    // Squares and cubes of real and imaginary parts used in the computations
    double
      xr2 = xr*xr, xi2 = xi*xi, yr2 = yr*yr, yi2 = yi*yi,
      xr3 = xr2*xr, xi3 = xi2*xi, yr3 = yr2*yr, yi3 = yi2*yi;

    // The first coordinate of the output is Re(x^3*y + y^3 + x)
    result(0) = 
      xr3*yr - 3*xr*xi2*yr - 3*xr2*xi*yi + xi3*yi
      + yr3 - 3*yr*yi2 + xr;
    // The second coordinate of the output is Im(x^3*y + y^3 + x)
    result(1) =
      3*xr2*xi*yr + xr3*yi - 3*xr*xi2*yi - xi3*yr
      + 3*yr2*yi - yi3 + xi;
    return result;
  }

  std::size_t amb_d() const {return 4;};
  std::size_t cod_d() const {return 2;};

  Eigen::VectorXd seed() const {
    Eigen::VectorXd result = Eigen::VectorXd::Zero(4);
    return result;
  }

  Function_surface_on_CP2_in_R4() {}  
};

int main(int argc, char** argv) {

  // The function for the (non-compact) manifold
  Function_surface_on_CP2_in_R4 fun;

  // Seed of the function
  Eigen::VectorXd seed = fun.seed();

  // Creating the function that defines the boundary of a compact region on the manifold
  double radius = 3.0;
  Function_Sm_in_Rd fun_sph(radius, 3, seed);

  // Defining the intersection oracle
  auto oracle = make_oracle(fun, fun_sph);

  // Define a Coxeter triangulation scaled by a factor lambda.
  // The triangulation is translated by a random vector to avoid violating the genericity hypothesis.
  double lambda = 0.2;
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
			 "manifold_on_CP2_with_boundary",
			 build_mesh_from_cell_complex(cell_complex,
						      Configuration(true, true, true, 1, 5, 3),
						      Configuration(true, true, true, 2, 13, 14)));
  return 0;
}
