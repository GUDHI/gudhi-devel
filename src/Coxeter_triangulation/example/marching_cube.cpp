#include <iostream>
#include <vector>

#include <gudhi/Points_off_io.h>
#include <gudhi/Ad_simplex.h>

#include <CGAL/Epick_d.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>

using K = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
using FT = K::FT;
using Point_d = K::Point_d;
using Point_vector = std::vector< Point_d >;

using Matrix = Eigen::SparseMatrix<FT>;
using Triplet = Eigen::Triplet<FT>;

/** Current state of the algorithm.
 *  Input: a point cloud 'point_vector'
 *  Output: a reconstruction (a simplicial complex?, a Czech-like complex?)
 */

int main(int argc, char * const argv[]) {
  std::cout << "Marching cube adaptation for Coxeter triangulations\n";
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0]
        << " path_to_off_point_file\n";
    return 0;
  }
  Point_vector point_vector;
  Gudhi::Points_off_reader<Point_d> off_reader(argv[1]);
  if (!off_reader.is_valid()) {
      std::cerr << "Coxeter triangulations - Unable to read file " << argv[1] << "\n";
      exit(-1);  // ----- >>
    }
  point_vector = Point_vector(off_reader.get_point_cloud());
  int N = point_vector.size();
  // short d = point_vector[0].size();
  short d = 2;
  std::cout << "Successfully read " << N << " points in dimension " << d << std::endl;

  // The A root vectors, computed as a matrix

  std::vector<Triplet> cartan_triplets;
  cartan_triplets.reserve(3*d-2);
  for (int i = 0; i < d; i++) {
    cartan_triplets.push_back(Triplet(i,i,2.0));
  }
  for (int i = 1; i < d; i++) {
    cartan_triplets.push_back(Triplet(i-1,i,-1.0));
    cartan_triplets.push_back(Triplet(i,i-1,-1.0));
  }
  Matrix cartan(d,d);
  cartan.setFromTriplets(cartan_triplets.begin(), cartan_triplets.end());
  std::cout << "cartan =" << std::endl << cartan << std::endl;

  Eigen::SimplicialLLT<Matrix, Eigen::Lower> chol(cartan);
  Matrix do_vect = chol.matrixL();
  std::cout << "do_vect^t =" << std::endl << do_vect << std::endl;
  
  std::vector<Triplet> r_t_triplets;
  r_t_triplets.reserve(2*d);
  for (int i = 0; i < d; i++) {
    r_t_triplets.push_back(Triplet(i,i,1.0));
  }
  for (int i = 0; i < d; i++) {
    r_t_triplets.push_back(Triplet(d,i,1.0));
  }
  Matrix r_t(d+1,d);
  r_t.setFromTriplets(r_t_triplets.begin(), r_t_triplets.end()); 
  std::cout << "r_t =" << std::endl << r_t << std::endl;

  Matrix norm_t = r_t * do_vect;
  std::cout << "norm_t =" << std::endl << norm_t << std::endl;
  

}
