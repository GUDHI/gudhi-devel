#include <gudhi/Points_off_io.h>

#include <iostream>
#include <string>
#include <vector>

using Point_d = std::vector<double>;

void usage(char * const progName) {
  std::cerr << "Usage: " << progName << " inputFile.off" << std::endl;
  exit(-1);
}

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "Error: Number of arguments (" << argc << ") is not correct" << std::endl;
    usage(argv[0]);
  }

  std::string offInputFile(argv[1]);
  // Read the OFF file (input file name given as parameter) and triangulate points
  Gudhi::Points_off_reader<Point_d> off_reader(offInputFile);
  // Check the read operation was correct
  if (!off_reader.is_valid()) {
    std::cerr << "Unable to read file " << offInputFile << std::endl;
    usage(argv[0]);
  }

  // Retrieve the triangulation
  std::vector<Point_d> point_cloud = off_reader.get_point_cloud();

  int n {0};
  for (auto point : point_cloud) {
    std::cout << "Point[" << n << "] = ";
    for (std::size_t i {0}; i < point.size(); i++)
      std::cout << point[i] << " ";
    std::cout << "\n";
    ++n;
  }
  return 0;
}
