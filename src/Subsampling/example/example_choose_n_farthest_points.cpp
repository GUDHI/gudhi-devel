#include <gudhi/choose_n_farthest_points.h>

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
  Gudhi::subsampling::choose_n_farthest_points(k, points, 100,
                                               Gudhi::subsampling::random_starting_point,
                                               std::back_inserter(results));
  std::cout << "Before sparsification: " << points.size() << " points.\n";
  std::cout << "After  sparsification: " << results.size() << " points.\n";

  return 0;
}
