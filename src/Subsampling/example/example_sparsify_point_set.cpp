#include <gudhi/sparsify_point_set.h>

#include <CGAL/Epick_d.h>
#include <CGAL/Random.h>

#include <vector>
#include <iterator>

int main (void)
{
  typedef CGAL::Epick_d<CGAL::Dimension_tag<4> >   K;
  typedef typename K::FT                           FT;
  typedef typename K::Point_d                      Point_d;
  
  CGAL::Random rd;

  std::vector<Point_d> points;
  for (int i = 0 ; i < 500 ; ++i)
    points.push_back(Point_d(rd.get_double(-1.,1),rd.get_double(-1.,1),rd.get_double(-1.,1),rd.get_double(-1.,1)));

  K k;
  std::vector<Point_d> results;
  Gudhi::subsampling::sparsify_point_set(k, points, 0.4, std::back_inserter(results));
  std::cout << "Before sparsification: " << points.size() << " points.\n";
  std::cout << "After  sparsification: " << results.size() << " points.\n";

  return 0;
}
