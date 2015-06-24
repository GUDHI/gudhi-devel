// to construct a Delaunay_triangulation from a OFF file
#include "gudhi/Delaunay_triangulation_off_io.h"
#include "gudhi/Alpha_complex.h"

#include <stdio.h>
#include <stdlib.h>
#include <string>

void usage(char * const progName) {
  std::cerr << "Usage: " << progName << " filename.off" << std::endl;
  exit(-1); // ----- >>
}

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "Error: Number of arguments (" << argc << ") is not correct" << std::endl;
    usage(argv[0]);
  }

  std::string off_file_name(argv[1]);

  // ----------------------------------------------------------------------------
  // Init of an alpha complex from an OFF file
  // ----------------------------------------------------------------------------
  Gudhi::alphacomplex::Alpha_complex alpha_complex_from_file(off_file_name);
  
  // ----------------------------------------------------------------------------
  // Display information about the alpha complex
  // ----------------------------------------------------------------------------
  std::cout << "Alpha complex is of dimension " << alpha_complex_from_file.dimension() <<
      " - " << alpha_complex_from_file.num_simplices() << " simplices - " <<
      alpha_complex_from_file.num_vertices() << " vertices." << std::endl;
  
  std::cout << "Iterator on alpha complex simplices in the filtration order, with [filtration value]:" << std::endl;
  for (auto f_simplex : alpha_complex_from_file.filtration_simplex_range()) {
    std::cout << "   ( ";
    for (auto vertex : alpha_complex_from_file.simplex_vertex_range(f_simplex)) {
      std::cout << vertex << " ";
    }
    std::cout << ") -> " << "[" << alpha_complex_from_file.filtration(f_simplex) << "] ";
    std::cout << std::endl;
  }
  return 0;
}