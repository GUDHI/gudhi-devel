#include <gudhi/choose_n_farthest_points.h>

#include <iostream>
#include <vector>
#include <iterator>


/* The class Kernel contains a distance function defined on the set of points {0, 1, 2, 3}
 * and computes a distance according to the matrix:
 * 0 1 2 4
 * 1 0 4 2
 * 2 4 0 1
 * 4 2 1 0
 */
class Kernel {
 public:
  typedef double FT;
  typedef unsigned Point_d;

  // Class Squared_distance_d
  class Squared_distance_d {
   private:
    std::vector<std::vector<FT>> matrix_;

   public:
    Squared_distance_d() {
      matrix_.push_back(std::vector<FT>({0, 1, 2, 4}));
      matrix_.push_back(std::vector<FT>({1, 0, 4, 2}));
      matrix_.push_back(std::vector<FT>({2, 4, 0, 1}));
      matrix_.push_back(std::vector<FT>({4, 2, 1, 0}));
    }

    FT operator()(Point_d p1, Point_d p2) {
      return matrix_[p1][p2];
    }
  };

  // Constructor
  Kernel() {}

  // Object of type Squared_distance_d
  Squared_distance_d squared_distance_d_object() const {
    return Squared_distance_d();
  }
};

int main(void) {
  typedef Kernel K;
  typedef typename K::Point_d Point_d;

  K k;
  std::vector<Point_d> points = {0, 1, 2, 3};
  std::vector<Point_d> results;

  Gudhi::subsampling::choose_n_farthest_points(k, points, 2,
                                               Gudhi::subsampling::random_starting_point,
                                               std::back_inserter(results));
  std::clog << "Before sparsification: " << points.size() << " points.\n";
  std::clog << "After  sparsification: " << results.size() << " points.\n";
  std::clog << "Result table: {" << results[0] << "," << results[1] << "}\n";

  return 0;
}
