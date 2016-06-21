// #ifdef _DEBUG
// # define TBB_USE_THREADING_TOOL
// #endif

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "witness_complex_points"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/choose_by_farthest_point.h>
#include <vector>
#include <iterator>

#include <CGAL/Epick_d.h>

typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag>                K;
typedef typename K::FT                                            FT;
typedef typename K::Point_d                                       Point_d;


BOOST_AUTO_TEST_CASE(test_choose_farthest_point) {
  std::vector< Point_d > points, landmarks;
  // Add grid points (625 points)
  for (FT i = 0; i < 5; i += 1.0)
    for (FT j = 0; j < 5; j += 1.0)
      for (FT k = 0; k < 5; k += 1.0)
        for (FT l = 0; l < 5; l += 1.0)
          points.push_back(Point_d(std::vector<FT>({i, j, k, l})));

  landmarks.clear();
  K k;
  Gudhi::choose_by_farthest_point(k, points, 100, std::back_inserter(landmarks));
  
  assert(landmarks.size() == 100);
}
