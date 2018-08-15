#define CC_STAR_COMPLETION
// #define CC_A_V_VISITORS

#include <iostream>
#include <vector>
#include <algorithm>

#include <gudhi/Points_off_io.h>
#include <gudhi/Coxeter_system.h>
#include <gudhi/Coxeter_complex.h>
#include <gudhi/Coxeter_complex/Off_point_range.h>
#include <gudhi/Clock.h>

#include <CGAL/Epick_d.h>
#include <CGAL/point_generators_d.h>

//#include <Eigen/Dense>

// #include "memory_usage.h"
// #include "cxx-prettyprint/prettyprint.hpp"
// #include "output_points_to_medit.h"

using K = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
using FT = K::FT;
using Point_d = K::Point_d;
using Point_vector = std::vector< Point_d >;
using Coxeter_complex = Gudhi::Coxeter_complex<Point_vector, Coxeter_system>;
using A_id = Coxeter_system::Alcove_id;

struct Alcove_compare {
  bool operator()(const A_id& a1, const A_id& a2) const {
    if (a1.size() < a2.size())
      return true;
    if (a1.size() > a2.size())
      return false;
    for (std::size_t k = 0; k < a1.size(); ++k)
      if (a1[k] < a2[k])
        return true;
      else if (a1[k] > a2[k])
        return false;
      else if (a1.is_fixed(k) && !a2.is_fixed(k))
        return true;
      else if (!a1.is_fixed(k) && a2.is_fixed(k))
        return false;
    return false;
  }
};

/** Recursive procedure that checks test1 for all products of triangulations ~A_i at a given dimension */
void rec_test1(std::vector<unsigned>& decomposition, Coxeter_system& cs, unsigned dimension) {
  if (dimension == 0) {
    std::cout << std::endl << cs;
    A_id a_id(1, cs.dimension());
    for (unsigned i = 0; i < cs.pos_root_count(); ++i)
      a_id.push_back(0);
    std::cout << "Cell " << a_id << ":\n";
    for (unsigned f_d = 0; f_d <= cs.dimension(); ++f_d) {
      std::vector<A_id> faces, faces2;
      Gudhi::Clock t;
      for (auto f_id: cs.face_range(a_id, f_d))
        faces.push_back(f_id);
      t.end();
      std::cout << "Computation time(old): " <<  t.num_seconds() << "s\n";
      t.begin();
      for (auto f_id: cs.face2_range(a_id, f_d))
        faces2.push_back(f_id);
      t.end();
      std::cout << "Computation time(new): " <<  t.num_seconds() << "s\n";
      std::sort(faces.begin(), faces.end(), Alcove_compare());
      std::sort(faces2.begin(), faces2.end(), Alcove_compare());
      if (faces != faces2)
        std::cerr << "The computed faces are not the same.";
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
  Coxeter_system cs('A', dimension);
  std::vector<double> total1(dimension + 1, 0), total2(dimension + 1, 0);
  CGAL::Random_points_on_sphere_d<Point_d> rp(dimension + 1, 5);
  unsigned num_tests = 5000;
  for (unsigned f_d = 0; f_d <= cs.dimension(); ++f_d) {
    Gudhi::Clock t;
    for (unsigned i = 0; i < num_tests; ++i) {
      A_id a_id = cs.query_point_location(*rp++, 1);
      for (auto f_id: cs.face_range(a_id, f_d)) {}
    }
    t.end();
    total1[f_d] += t.num_seconds() / num_tests * 1000;
  }
  for (unsigned f_d = 0; f_d <= cs.dimension(); ++f_d) {
    Gudhi::Clock t;
    for (unsigned i = 0; i < num_tests; ++i) {
      A_id a_id = cs.query_point_location(*rp++, 1);
      for (auto f_id: cs.face2_range(a_id, f_d)) {}
    }
    t.end();
    total2[f_d] += t.num_seconds() / num_tests * 1000;
  }

  std::cout << "\\hline\n Face dimension ";
  for (auto i = 0; i <= dimension; ++i)
    std::cout << "& " << i << " ";
  std::cout << "\\\\\n";
  std::cout << "\\hline\\hline\nOld algorithm\n";
  for (auto t: total1)
    std::cout << "& " << t << " ";
  std::cout << "\\\\\n\\hline\n";
  std::cout << "New algorithm\n";
  for (auto t: total2)
    std::cout << "& " << t << " ";
  std::cout << "\\\\\n\\hline\n";
  // rec_test1(decomposition, cs, dimension);
  // cs.emplace_back('A', 2);
  // typename Coxeter_system::Alcove_id a_id(1, 0);
  // a_id.push_back(0, false);
  // a_id.push_back(0, false);
  // a_id.push_back(0, false);
  // std::cout << a_id << "\n";
  // for (auto f_id: cs.face2_range(a_id, 0))
  //   std::cout <<  " " << f_id << "\n";
  // for (auto f_id: cs.simple_coxeter_system_begin()->coface_range(a_id, 0))
  //   std::cout <<  " " << f_id << "\n";
}
