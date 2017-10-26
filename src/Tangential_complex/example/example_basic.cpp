#include <gudhi/Tangential_complex.h>
#include <gudhi/sparsify_point_set.h>
#include <gudhi/Fake_simplex_tree.h>


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
  const int INTRINSIC_DIM = 2;
  const int AMBIENT_DIM = 3;
  const int NUM_POINTS = 100;

  Kernel k;

  // Generate points on a 2-sphere
  CGAL::Random_points_on_sphere_d<Point> generator(AMBIENT_DIM, 3.);
  std::vector<Point> points;
  points.reserve(NUM_POINTS);
  for (int i = 0; i < NUM_POINTS; ++i)
    points.push_back(*generator++);

  // Compute the TC
  TC tc(points, INTRINSIC_DIM, k);
  tc.compute_tangential_complex();

  // Export the TC into a Simplex_tree
  //Gudhi::Simplex_tree<> stree;
  Gudhi::Fake_simplex_tree stree;
  tc.create_complex(stree);

  // Display stats about inconsistencies
  tc.number_of_inconsistent_simplices(true);  // verbose

  return 0;
}
