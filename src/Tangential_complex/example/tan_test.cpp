#include <gudhi/Tangential_complex.h>
#include <gudhi/sparsify_point_set.h>

#include <CGAL/Epick_d.h>
#include <CGAL/Random.h>

#include <array>
#include <vector>

namespace tc = Gudhi::tangential_complex;

typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_d Point;
typedef Kernel::Vector_d Vector;
typedef tc::Tangential_complex<
Kernel, CGAL::Dynamic_dimension_tag,
CGAL::Parallel_tag> TC;

int main(void) {
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

  Kernel k;
  for (int i = 0; i < 100; i++) {

    // Compute the TC
    TC tc(points, INTRINSIC_DIM, k);
    tc.compute_tangential_complex();
    TC::Num_inconsistencies num_inc = tc.number_of_inconsistent_simplices();
    std::cout << "TC vertices = " << tc.number_of_vertices() << " - simplices = " << num_inc.num_simplices << 
        " - inconsistencies = " << num_inc.num_inconsistent_simplices << std::endl;

    tc.fix_inconsistencies_using_perturbation(10.0, 60.0);
    // Export the TC into a Simplex_tree
    Gudhi::Simplex_tree<> stree;
    tc.create_complex(stree);
    std::cout << "ST vertices = " << stree.num_vertices() << " - simplices = " << stree.num_simplices() << std::endl;

  }

  return 0;
}
