/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/Debug_utils.h>  // for GUDHI_CHECK

#include <boost/range/metafunctions.hpp>
#include <boost/range/size.hpp>

#include <iterator>  // for std::begin, std::end
#include <utility>  // for std::pair
#include <cmath>  // for std::abs
#include <vector>
#include <iostream>


#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "compute_edge_graph"
#include <boost/test/unit_test.hpp>

class Manhattan_distance {
 public:
  // boost::range_value is not SFINAE-friendly so we cannot use it in the return type
  template< typename Point >
  typename std::iterator_traits<typename boost::range_iterator<Point>::type>::value_type
  operator()(const Point& p1, const Point& p2) const {
    auto it1 = std::begin(p1);
    auto it2 = std::begin(p2);
    typedef typename boost::range_value<Point>::type NT;
    NT dist = 0;
    for (; it1 != std::end(p1); ++it1, ++it2) {
      GUDHI_CHECK(it2 != std::end(p2), "inconsistent point dimensions");
      NT tmp = *it1 - *it2;
      dist += std::abs(tmp);
    }
    GUDHI_CHECK(it2 == std::end(p2), "inconsistent point dimensions");
    return dist;
  }
};

struct Simplicial_complex {
  using Filtration_value = double;
  using Vertex_handle = int;
};

BOOST_AUTO_TEST_CASE( compute_manhattan_edge_graph )
{
  using Point = std::vector<Simplicial_complex::Filtration_value>;
  using Point_cloud = std::vector<Point>;

  Point_cloud point_cloud = {{1.0, 1.0 },
                             {7.0, 0.0 },
                             {4.0, 6.0 },
                             {9.0, 6.0 },
                             {0.0, 14.0},
                             {2.0, 19.0},
                             {9.0, 17.0}};
  
  Gudhi::Filtered_edges_container<Simplicial_complex> graph =
    Gudhi::compute_edge_graph<Gudhi::Filtered_edges_container, Simplicial_complex>(
    point_cloud,
    std::numeric_limits<typename Simplicial_complex::Filtration_value>::infinity(),
    Manhattan_distance());

  // Check it is sorted by filtration values
  Simplicial_complex::Filtration_value prev_filt = 0.;
  auto edges = graph.edges();
  for (auto value: edges) {
    std::cout << "filtration = " << std::get<0>(value) <<
                 " - edge = (" << std::get<1>(value) << ", " << std::get<2>(value) << ")" << std::endl;
    BOOST_CHECK(std::get<0>(value) >= prev_filt);
    prev_filt = std::get<0>(value);
  }

  BOOST_CHECK(graph.size() == 21);
  BOOST_CHECK(graph.get_filtration_min() == 5.);
  BOOST_CHECK(graph.get_filtration_max() == 24.);
  BOOST_CHECK(graph.get_filtration_at(0) == 5.);
  BOOST_CHECK(graph.get_filtration_at(20) == 24.);

  auto sub_graph_by_filt = graph.sub_filter_edges_by_filtration(17.);
  BOOST_CHECK(sub_graph_by_filt.size() == 14);
  auto sub_graph_by_idx = graph.sub_filter_edges_by_index(13);
  BOOST_CHECK(sub_graph_by_filt == sub_graph_by_idx);
}
