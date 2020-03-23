#include <gudhi/Kd_tree_search.h>

#include <CGAL/Epick_d.h>
#include <CGAL/Random.h>

#include <vector>

namespace gss = Gudhi::spatial_searching;

int main(void) {
  typedef CGAL::Epick_d<CGAL::Dimension_tag<4> > K;
  typedef typename K::Point_d Point;
  typedef std::vector<Point> Points;

  typedef gss::Kd_tree_search<K, Points> Points_ds;

  CGAL::Random rd;

  Points points;
  for (int i = 0; i < 500; ++i)
    points.push_back(Point(rd.get_double(-1., 1), rd.get_double(-1., 1), rd.get_double(-1., 1), rd.get_double(-1., 1)));

  Points_ds points_ds(points);

  // 10-nearest neighbor query
  std::clog << "10 nearest neighbors from points[20]:\n";
  auto knn_range = points_ds.k_nearest_neighbors(points[20], 10, true);
  for (auto const& nghb : knn_range)
    std::clog << nghb.first << " (sq. dist. = " << nghb.second << ")\n";

  // Incremental nearest neighbor query
  std::clog << "Incremental nearest neighbors:\n";
  auto inn_range = points_ds.incremental_nearest_neighbors(points[45]);
  // Get the neighbors in distance order until we hit the first point
  for (auto ins_iterator = inn_range.begin(); ins_iterator->first != 0; ++ins_iterator)
    std::clog << ins_iterator->first << " (sq. dist. = " << ins_iterator->second << ")\n";

  // 10-furthest neighbor query
  std::clog << "10 furthest neighbors from points[20]:\n";
  auto kfn_range = points_ds.k_furthest_neighbors(points[20], 10, true);
  for (auto const& nghb : kfn_range)
    std::clog << nghb.first << " (sq. dist. = " << nghb.second << ")\n";

  // Incremental furthest neighbor query
  std::clog << "Incremental furthest neighbors:\n";
  auto ifn_range = points_ds.incremental_furthest_neighbors(points[45]);
  // Get the neighbors in distance reverse order until we hit the first point
  for (auto ifs_iterator = ifn_range.begin(); ifs_iterator->first != 0; ++ifs_iterator)
    std::clog << ifs_iterator->first << " (sq. dist. = " << ifs_iterator->second << ")\n";

  // All-near-neighbors search
  std::clog << "All-near-neighbors search:\n";
  std::vector<std::size_t> rs_result;
  points_ds.all_near_neighbors(points[45], 0.5, std::back_inserter(rs_result));
  K k;
  for (auto const& p_idx : rs_result)
    std::clog << p_idx << " (sq. dist. = " << k.squared_distance_d_object()(points[p_idx], points[45]) << ")\n";

  return 0;
}
