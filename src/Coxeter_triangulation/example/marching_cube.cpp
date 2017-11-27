#include <iostream>
#include <vector>

#include <gudhi/Points_off_io.h>

#include <CGAL/Epick_d.h>

using K = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
using Point_d = K::Point_d;
using Point_vector = std::vector< Point_d >;

int main(int argc, char * const argv[]) {
  std::cout << "Marching cube adaptation for Coxeter triangulations\n";
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0]
        << " path_to_off_point_file\n";
    return 0;
  }
  Point_vector point_vector;
  Gudhi::Points_off_reader<Point_d> off_reader(argv[1]);
  if (!off_reader.is_valid()) {
      std::cerr << "Coxeter triangulations - Unable to read file " << argv[1] << "\n";
      exit(-1);  // ----- >>
    }
  point_vector = Point_vector(off_reader.get_point_cloud());
  std::cout << "Successfully read " << point_vector.size() << " points in dimension " << point_vector[0].size() << std::endl;
}
