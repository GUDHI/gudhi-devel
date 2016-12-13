#include <gudhi/Tangential_complex.h>
#include <gudhi/sparsify_point_set.h>

#include <CGAL/Epick_d.h>
#include <CGAL/Random.h>

#include <array>
#include <vector>

namespace tc = Gudhi::tangential_complex;


int main(void) {
  typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> Kernel;
  typedef Kernel::FT FT;
  typedef Kernel::Point_d Point;
  typedef Kernel::Vector_d Vector;
  typedef tc::Tangential_complex<Kernel, CGAL::Dynamic_dimension_tag, CGAL::Parallel_tag> TC;


  const int INTRINSIC_DIM = 1;

  // Generate points on a 2-sphere
  std::vector<Point> points;
  // [[0, 0], [1, 0], [0, 1], [1, 1]]
  std::vector<double> point = {0.0, 0.0};
  points.push_back(Point(point.size(), point.begin(), point.end()));
  point = {1.0, 0.0};
  points.push_back(Point(point.size(), point.begin(), point.end()));
  point = {0.0, 1.0};
  points.push_back(Point(point.size(), point.begin(), point.end()));
  point = {1.0, 1.0};
  points.push_back(Point(point.size(), point.begin(), point.end()));
  std::cout << "points = " << points.size() << std::endl;
  Kernel k;

  // Compute the TC
  TC tc(points, INTRINSIC_DIM, k);
  tc.compute_tangential_complex();
  TC::Num_inconsistencies num_inc = tc.number_of_inconsistent_simplices();
  std::cout << "TC vertices = " << tc.number_of_vertices() << " - simplices = " << num_inc.num_simplices <<
               " - inconsistencies = " << num_inc.num_inconsistent_simplices << std::endl;

  // Export the TC into a Simplex_tree
  Gudhi::Simplex_tree<> stree;
  tc.create_complex(stree);

  std::cout << "********************************************************************\n";
  std::cout << "* The complex contains " << stree.num_simplices() << " simplices";
  std::cout << "   - dimension " << stree.dimension() << "   - filtration " << stree.filtration() << "\n";
  std::cout << "* Iterator on Simplices in the filtration, with [filtration value]:\n";
  for (auto f_simplex : stree.filtration_simplex_range()) {
    std::cout << "   " << "[" << stree.filtration(f_simplex) << "] ";
    for (auto vertex : stree.simplex_vertex_range(f_simplex)) {
      std::cout << static_cast<int>(vertex) << " ";
    }
    std::cout << std::endl;
  }

  tc.fix_inconsistencies_using_perturbation(0.01, 30.0);

  // Export the TC into a Simplex_tree
  tc.create_complex(stree);

  std::cout << "********************************************************************\n";
  std::cout << "* The complex contains " << stree.num_simplices() << " simplices\n";
  std::cout << "   - dimension " << stree.dimension() << "   - filtration " << stree.filtration() << "\n";
  std::cout << "* Iterator on Simplices in the filtration, with [filtration value]:\n";
  for (auto f_simplex : stree.filtration_simplex_range()) {
    std::cout << "   " << "[" << stree.filtration(f_simplex) << "] ";
    for (auto vertex : stree.simplex_vertex_range(f_simplex)) {
      std::cout << static_cast<int>(vertex) << " ";
    }
    std::cout << std::endl;
  }

  return 0;
}
