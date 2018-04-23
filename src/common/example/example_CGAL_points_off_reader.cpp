#include <gudhi/Points_off_io.h>

// For CGAL points type in dimension d
// cf. http://doc.cgal.org/latest/Kernel_d/classCGAL_1_1Point__d.html
#include <CGAL/Epick_d.h>

#include <iostream>
#include <string>
#include <vector>

using Kernel = CGAL::Epick_d< CGAL::Dynamic_dimension_tag >;
using Point_d = Kernel::Point_d;

void usage(char * const progName) {
  std::cerr << "Usage: " << progName << " inputFile.off" << std::endl;
  exit(-1);
}

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "Error: Number of arguments (" << argc << ") is not correct" << std::endl;
    usage(argv[0]);
  }

  std::string off_input_file(argv[1]);
  // Read the OFF file (input file name given as parameter) and triangulate points
  Gudhi::Points_off_reader<Point_d> off_reader(off_input_file);
  // Check the read operation was correct
  if (!off_reader.is_valid()) {
    std::cerr << "Unable to read file " << off_input_file << std::endl;
    usage(argv[0]);
  }

  // Retrieve the triangulation
  std::vector<Point_d> point_cloud = off_reader.get_point_cloud();

  int n {};
  for (auto point : point_cloud) {
    std::cout << "Point[" << n << "] = ";
    for (std::size_t i {0}; i < point.size(); i++)
      std::cout << point[i] << " ";
    std::cout << "\n";
    ++n;
  }
  return 0;
}
