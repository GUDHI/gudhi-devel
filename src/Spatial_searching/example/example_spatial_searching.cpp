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

  CGAL::Random rd;

  Points points;
  for (int i = 0; i < 500; ++i)
    points.push_back(Point(std::array<FT, 4>({ rd.get_double(-1.,1),rd.get_double(-1.,1),rd.get_double(-1.,1),rd.get_double(-1.,1) })));

  Points_ds points_ds(points);

  // 20-nearest neighbor query
  std::cout << "20 nearest neighbors:\n";
  auto kns_range = points_ds.query_k_nearest_neighbors(points[20], 10, true);
  for (auto const& nghb : kns_range)
    std::cout << nghb.first << " (sq. dist. = " << nghb.second << ")\n";

  // Incremental nearest neighbor query
  std::cout << "Incremental nearest neighbors:\n";
  auto ins_range = points_ds.query_incremental_nearest_neighbors(points[45]);
  // Get all the neighbors that are closer than 0.5
  for (auto ins_iterator = ins_range.begin(); ins_iterator->second < 0.5*0.5 ; ++ins_iterator)
    std::cout << ins_iterator->first << " (sq. dist. = " << ins_iterator->second << ")\n";

  // 20-farthest neighbor query
  std::cout << "20 farthest neighbors:\n";
  auto kfs_range = points_ds.query_k_farthest_neighbors(points[20], 10, true);
  for (auto const& nghb : kfs_range)
    std::cout << nghb.first << " (sq. dist. = " << nghb.second << ")\n";

  // Incremental farthest neighbor query
  std::cout << "Incremental farthest neighbors:\n";
  auto ifs_range = points_ds.query_incremental_farthest_neighbors(points[45]);
  // Get all the neighbors that are farthest than 2.3
  for (auto ifs_iterator = ifs_range.begin(); ifs_iterator->second > 2.3*2.3 ; ++ifs_iterator)
    std::cout << ifs_iterator->first << " (sq. dist. = " << ifs_iterator->second << ")\n";
  
  return 0;
}
