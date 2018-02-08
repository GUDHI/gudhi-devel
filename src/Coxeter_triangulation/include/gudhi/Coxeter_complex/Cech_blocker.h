#ifndef CECH_BLOCKER_H_
#define CECH_BLOCKER_H_

#include <CGAL/Epick_d.h>
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Min_sphere_of_points_d_traits_d.h>

namespace Gudhi {

template<class SimplexTree, unsigned dimension_>
class Cech_blocker {
  using Simplex_handle = typename SimplexTree::Simplex_handle;
  using Filtration_value = typename SimplexTree::Filtration_value;
  using Kernel = CGAL::Epick_d< CGAL::Dimension_tag<dimension_> >;// CGAL::Dynamic_dimension_tag >;
  using Traits = CGAL::Min_sphere_of_points_d_traits_d<Kernel, Filtration_value, dimension_>;
  using Min_sphere = CGAL::Min_sphere_of_spheres_d<Traits>;
  using Point = typename Kernel::Point_d;
 public:
  bool operator()(Simplex_handle sh) {
    std::vector<Point> points;
    for (auto vertex : simplex_tree_.simplex_vertex_range(sh)) {
      points.push_back(point_cloud_[vertex]);
#if DEBUG_TRACES
      std::cout << "#(" << vertex << ")#";
#endif  // DEBUG_TRACES
    }
    Min_sphere ms(points.begin(),points.end());
    Filtration_value radius = ms.radius();
#if DEBUG_TRACES
    std::cout << "radius = " << radius << " - " << (radius > threshold_) << std::endl;
#endif  // DEBUG_TRACES
    simplex_tree_.assign_filtration(sh, radius);
    return (radius > threshold_);
  }
  Cech_blocker(SimplexTree& simplex_tree, Filtration_value threshold, const std::vector<Point>& point_cloud)
    : simplex_tree_(simplex_tree),
      threshold_(threshold),
      point_cloud_(point_cloud) { }
 private:
  SimplexTree simplex_tree_;
  Filtration_value threshold_;
  std::vector<Point> point_cloud_;
};

}
#endif
