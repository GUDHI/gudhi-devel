#include <gudhi/pick_n_random_points.h>

#include <CGAL/Epick_d.h>
#include <CGAL/Random.h>

#include <iostream>
#include <vector>
#include <iterator>

int main(void) {
  typedef CGAL::Epick_d<CGAL::Dimension_tag<4> > K;
  typedef typename K::Point_d Point_d;

  CGAL::Random rd;

  std::vector<Point_d> points;
  for (int i = 0; i < 500; ++i)
    points.push_back(Point_d(rd.get_double(-1., 1), rd.get_double(-1., 1),
                             rd.get_double(-1., 1), rd.get_double(-1., 1)));

  K k;
  std::vector<Point_d> results;
  Gudhi::subsampling::pick_n_random_points(points, 100, std::back_inserter(results));
  std::clog << "Before sparsification: " << points.size() << " points.\n";
  std::clog << "After  sparsification: " << results.size() << " points.\n";

  return 0;
}
