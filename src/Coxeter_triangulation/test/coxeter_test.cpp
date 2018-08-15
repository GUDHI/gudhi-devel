#define CC_STAR_COMPLETION
// #define CC_A_V_VISITORS

#include <iostream>
#include <vector>

#include <gudhi/Points_off_io.h>
#include <gudhi/Coxeter_system.h>
#include <gudhi/Coxeter_complex.h>
#include <gudhi/Coxeter_complex/Off_point_range.h>
#include <gudhi/Clock.h>

#include <CGAL/Epick_d.h>

//#include <Eigen/Dense>

// #include "memory_usage.h"
// #include "cxx-prettyprint/prettyprint.hpp"
// #include "output_points_to_medit.h"

using K = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
using FT = K::FT;
using Point_d = K::Point_d;
using Point_vector = std::vector< Point_d >;
using Coxeter_complex = Gudhi::Coxeter_complex<Point_vector, Coxeter_system>;

/** Recursive procedure that checks test1 for all products of triangulations ~A_i at a given dimension */
void rec_test1(std::vector<unsigned>& decomposition, Coxeter_system& cs, unsigned dimension) {
  if (dimension == 0) {
    std::cout << std::endl << cs;
    typedef typename Coxeter_system::Alcove_id A_id;
    A_id a_id(1, cs.dimension());
    for (unsigned i = 0; i < cs.pos_root_count(); ++i)
      a_id.push_back(0);
    std::cout << "Cell " << a_id << ":\n";
    for (unsigned f_d = 0; f_d <= cs.dimension(); ++f_d) {
      unsigned total_faces_count = 0;
      std::cout << "Faces of dimension " << f_d << ":\n";
      for (auto f_id: cs.face_range(a_id, f_d)) {
        std::cout << " " << f_id  << "\n";
        // std::cout << " Faces of " << f_id << ":\n";
        // for (unsigned ff_d = 0; ff_d <= f_d; ++ff_d)
        //   for (auto ff_id: cs.face_range(f_id, ff_d))
        //     std::cout << "  " << ff_id  << "\n";
        std::cout << " Cofaces of " << f_id << ":\n";
        for (unsigned ff_d = f_d; ff_d <= cs.dimension(); ++ff_d)
          for (auto ff_id: cs.coface_range(f_id, ff_d))
            std::cout << "  " << ff_id  << "\n";
        total_faces_count++;
      }
      std::cout << "Total number of faces of dimension " << f_d << " is " << total_faces_count << ".\n";
    }
    return;
  }
  unsigned i = decomposition.back();
  if (decomposition.back() == 0)
    i = 1;
  for (; i <= dimension; ++i) {
    decomposition.push_back(i);
    cs.emplace_back('A', i);
    rec_test1(decomposition, cs, dimension-i);
    cs.pop_back();
    decomposition.pop_back();
  }
}

/** Test suite to check features of the code */ 
int main(int argc, char * const argv[]) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0]
              << " dimension\n";
    return 0;
  }
  unsigned dimension = atoi(argv[1]);
  /* Test1: print faces and cofaces of the simplex [0,...,0] */
  std::vector<unsigned> decomposition; // first coordinate is the sum
  decomposition.reserve(dimension);
  Coxeter_system cs;
  rec_test1(decomposition, cs, dimension);
  // cs.emplace_back('A', 2);
  // typename Coxeter_system::Alcove_id a_id(1, 0);
  // a_id.push_back(0, true);
  // a_id.push_back(0, true);
  // a_id.push_back(0, true);
  // std::cout << a_id << "\n";
  // for (auto f_id: cs.face_range(a_id, 0))
  //   std::cout <<  " " << f_id << "\n";
  // for (auto f_id: cs.simple_coxeter_system_begin()->coface_range(a_id, 0))
  //   std::cout <<  " " << f_id << "\n";
}
