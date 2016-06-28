#include <gudhi/Spatial_tree_data_structure.h>

#include <CGAL/Epick_d.h>
#include <CGAL/Random.h>

#include <array>
#include <vector>

namespace gss = Gudhi::spatial_searching;

int main (void)
{
  typedef CGAL::Epick_d<CGAL::Dimension_tag<4> >      K;
  typedef typename K::FT                              FT;
  typedef typename K::Point_d                         Point;
  typedef std::vector<Point>                          Points;

  typedef gss::Spatial_tree_data_structure<K, Points> Points_ds;
  typedef typename Points_ds::KNS_range               KNS_range;
  typedef typename Points_ds::KNS_iterator            KNS_iterator;
  typedef typename Points_ds::INS_range               INS_range;
  typedef typename Points_ds::INS_iterator            INS_iterator;

  CGAL::Random rd;

  Points points;
  for (int i = 0; i < 500; ++i)
    points.push_back(Point(std::array<FT, 4>({ rd.get_double(-1.,1),rd.get_double(-1.,1),rd.get_double(-1.,1),rd.get_double(-1.,1) })));

  Points_ds points_ds(points);

  auto kns_range = points_ds.query_ANN(points[20], 10, true);
  for (auto const& nghb : kns_range)
    std::cout << nghb.first << " (sq. dist. = " << nghb.second << ")\n";

  return 0;
}
