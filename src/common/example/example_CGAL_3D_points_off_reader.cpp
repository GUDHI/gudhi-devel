#include <gudhi/Points_3D_off_io.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <iostream>
#include <string>
#include <vector>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = Kernel::Point_3;

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
  Gudhi::Points_3D_off_reader<Point_3> off_reader(off_input_file);
  // Check the read operation was correct
  if (!off_reader.is_valid()) {
    std::cerr << "Unable to read file " << off_input_file << std::endl;
    usage(argv[0]);
  }

  // Retrieve the triangulation
  std::vector<Point_3> point_cloud = off_reader.get_point_cloud();

  int n {};
  for (auto point : point_cloud) {
    ++n;
    std::clog << "Point[" << n << "] = (" << point[0] << ", " << point[1] << ", " << point[2] << ")\n";
  }
  return 0;
}
