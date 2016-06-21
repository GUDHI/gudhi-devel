// #ifdef _DEBUG
// # define TBB_USE_THREADING_TOOL
// #endif

#include <gudhi/pick_random_points.h>
#include <gudhi/choose_by_farthest_point.h>
#include <vector>
#include <iterator>

#include <CGAL/Epick_d.h>


int main() {
  typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag>                K;
  typedef typename K::FT                                            FT;
  typedef typename K::Point_d                                       Point_d;
  
  std::vector<Point_d> vect;
  vect.push_back(Point_d(std::vector<FT>({0,0,0,0})));
  vect.push_back(Point_d(std::vector<FT>({0,0,0,1})));
  vect.push_back(Point_d(std::vector<FT>({0,0,1,0})));
  vect.push_back(Point_d(std::vector<FT>({0,0,1,1})));
  vect.push_back(Point_d(std::vector<FT>({0,1,0,0})));
  vect.push_back(Point_d(std::vector<FT>({0,1,0,1})));
  vect.push_back(Point_d(std::vector<FT>({0,1,1,0})));
  vect.push_back(Point_d(std::vector<FT>({0,1,1,1})));
  vect.push_back(Point_d(std::vector<FT>({1,0,0,0})));
  vect.push_back(Point_d(std::vector<FT>({1,0,0,1})));
  vect.push_back(Point_d(std::vector<FT>({1,0,1,0})));
  vect.push_back(Point_d(std::vector<FT>({1,0,1,1})));
  vect.push_back(Point_d(std::vector<FT>({1,1,0,0})));
  vect.push_back(Point_d(std::vector<FT>({1,1,0,1})));
  vect.push_back(Point_d(std::vector<FT>({1,1,1,0})));
  vect.push_back(Point_d(std::vector<FT>({1,1,1,1})));
  
  
  std::vector<Point_d> landmarks;
  Gudhi::pick_random_points(vect, 5, std::back_inserter(landmarks));
  std::cout << "landmark vector contains: ";
  for (auto l: landmarks)
    std::cout << l << "\n";
  assert(landmarks_size() == 5);
  
}
