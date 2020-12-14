#include <gudhi/choose_n_farthest_points.h>

#include <iostream>
#include <vector>
#include <iterator>


typedef unsigned Point;

/* The class Distance contains a distance function defined on the set of points {0, 1, 2, 3}
 * and computes a distance according to the matrix:
 * 0 1 2 4
 * 1 0 4 2
 * 2 4 0 1
 * 4 2 1 0
 */
class Distance {
   private:
    std::vector<std::vector<double>> matrix_;

   public:
    Distance() {
      matrix_.push_back({0, 1, 2, 4});
      matrix_.push_back({1, 0, 4, 2});
      matrix_.push_back({2, 4, 0, 1});
      matrix_.push_back({4, 2, 1, 0});
    }

    double operator()(Point p1, Point p2) const {
      return matrix_[p1][p2];
    }
};

int main(void) {
  std::vector<Point> points = {0, 1, 2, 3};
  std::vector<Point> results;

  Gudhi::subsampling::choose_n_farthest_points(Distance(), points, 2,
                                               Gudhi::subsampling::random_starting_point,
                                               std::back_inserter(results));
  std::clog << "Before sparsification: " << points.size() << " points.\n";
  std::clog << "After  sparsification: " << results.size() << " points.\n";
  std::clog << "Result table: {" << results[0] << "," << results[1] << "}\n";
}
