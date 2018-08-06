#define CC_STAR_COMPLETION
// #define CC_A_V_VISITORS

#include <iostream>
#include <vector>
#include <fstream>

// #include <gudhi/Points_off_io.h>
// #include <gudhi/Coxeter_system.h>
// #include <gudhi/Coxeter_complex.h>
// #include <gudhi/Coxeter_complex/Off_point_range.h>
// #include <gudhi/Clock.h>
#include <gudhi/Hasse_diagram_persistence.h>

#include <CGAL/Epick_d.h>

// #include "memory_usage.h"
// #include "cxx-prettyprint/prettyprint.hpp"
// #include "output_points_to_medit.h"

using K = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
// using FT = K::FT;
using Point_d = K::Point_d;
// using Point_vector = std::vector< Point_d >;
// using Coxeter_complex = Gudhi::Coxeter_complex<Point_vector, Coxeter_system>;

void seed_expansion(const Point_d& p) {

}

int main(int argc, char * const argv[]) {
  std::vector<Point_d> seed_points;
  for (const auto& p: seed_points)
    seed_expansion(p);
}
