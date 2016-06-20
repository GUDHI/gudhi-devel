// #ifdef _DEBUG
// # define TBB_USE_THREADING_TOOL
// #endif

#include <gudhi/sparsify_point_set.h>

#include <CGAL/Epick_d.h>
#include <CGAL/Random.h>

#include <array>
#include <vector>
#include <iterator>

int main() {
  typedef CGAL::Epick_d<CGAL::Dimension_tag<4> >   K;
  typedef typename K::FT                           FT;
  typedef typename K::Point_d                      Point_d;
  
  CGAL::Random rd;

  std::vector<Point_d> points;
  for (int i = 0 ; i < 500 ; ++i)
    points.push_back(Point_d(std::array<FT,4>({rd.get_double(-1.,1),rd.get_double(-1.,1),rd.get_double(-1.,1),rd.get_double(-1.,1)})));

  K k;
  std::vector<Point_d> results;
  Gudhi::subsampling::sparsify_point_set(k, points, 0.5, std::back_inserter(results));
  std::cout << "Before sparsification: " << points.size() << " points.\n";
  std::cout << "After  sparsification: " << results.size() << " points.\n";
  //for (auto p : results)
  //  std::cout << p << "\n";

  return 0;
}
