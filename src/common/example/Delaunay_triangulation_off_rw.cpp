// to construct a Delaunay_triangulation from a OFF file
#include <gudhi/Delaunay_triangulation_off_io.h>

#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/Epick_d.h>

#include <iostream>
#include <string>

// Use dynamic_dimension_tag for the user to be able to set dimension
typedef CGAL::Epick_d< CGAL::Dynamic_dimension_tag > K;
typedef CGAL::Delaunay_triangulation<K> T;
// The triangulation uses the default instantiation of the 
// TriangulationDataStructure template parameter

void usage(char * const progName) {
  std::cerr << "Usage: " << progName << " inputFile.off outputFile.off" << std::endl;
  exit(-1); // ----- >>
}

int main(int argc, char **argv) {
  if (argc != 3) {
    std::cerr << "Error: Number of arguments (" << argc << ") is not correct" << std::endl;
    usage(argv[0]);
  }


#ifdef GUDHI_NDEBUG
  std::cout << "pouet pouet !!" << std::endl;
#endif

  std::string offInputFile(argv[1]);
  // Read the OFF file (input file name given as parameter) and triangulates points
  Gudhi::Delaunay_triangulation_off_reader<T> off_reader(offInputFile);
  // Check the read operation was correct
  if (!off_reader.is_valid()) {
    std::cerr << "Unable to read file " << offInputFile << std::endl;
    exit(-1); // ----- >>
  }
  
  // Retrieve the triangulation
  T* triangulation = off_reader.get_complex();
  // Operations on triangulation
  std::cout << "Number of vertices= " << triangulation->number_of_vertices() << std::endl;
  std::cout << "Number of finite full cells= " << triangulation->number_of_finite_full_cells() << std::endl;

  std::string outFileName(argv[2]);
  std::string offOutputFile(outFileName);
  // Write the OFF file (output file name given as parameter) with the points and triangulated cells as faces
  Gudhi::Delaunay_triangulation_off_writer<T> off_writer(offOutputFile, triangulation);

  // Check the write operation was correct
  if (!off_writer.is_valid()) {
    std::cerr << "Unable to write file " << offOutputFile << std::endl;
    exit(-1); // ----- >>
  }
  
  return 0;
}