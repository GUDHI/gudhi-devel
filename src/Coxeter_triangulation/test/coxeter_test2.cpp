#define CC_STAR_COMPLETION
// #define CC_A_V_VISITORS

#include <iostream>
#include <vector>
#include <algorithm>

#include <gudhi/Points_off_io.h>
#include <gudhi/Coxeter_system.h>
#include <gudhi/Coxeter_system2.h>
#include <gudhi/Coxeter_complex.h>
#include <gudhi/Coxeter_complex/Off_point_range.h>
#include <gudhi/Clock.h>

#include <boost/math/special_functions/binomial.hpp>

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

std::size_t BinomialCoefficient(std::size_t n, std::size_t k) {
  if (k == 0) { return 1; }
  else { return (n * BinomialCoefficient(n - 1, k - 1)) / k; }
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
  std::vector<std::vector<double> >
    total1(dimension + 1, std::vector<double>(dimension + 1, 0)),
    total2(dimension + 1, std::vector<double>(dimension + 1, 0));
  CGAL::Random_points_on_sphere_d<Point_d> rp(dimension, 5);
  unsigned num_tests = 1000;
  std::vector<std::vector<A_id> > faces(dimension + 1);
  for (unsigned i = 0; i < num_tests; ++i) {
    A_id a_id = cs.simple_coxeter_system_begin()->query_point_location(*rp++, 1);
    for (unsigned f_d = 0; f_d <= dimension; ++f_d) 
      for (auto f_id: Coxeter_system2('A', dimension).face_range(a_id, f_d))
        faces[f_d].push_back(f_id);
  }
  for (unsigned f_d = 0; f_d <= dimension; ++f_d) {
    for (unsigned ff_d = 0; ff_d <= f_d; ++ff_d) {
      Gudhi::Clock t;
      for (auto f_id: faces[f_d])
        for (auto ff_id: cs.face2_range(f_id, ff_d)) {}
      t.end();
      total2[f_d][ff_d] += t.num_seconds() / faces[f_d].size() * 1000;
    }
  }
  for (unsigned f_d = 0; f_d <= dimension; ++f_d) {
    for (unsigned ff_d = 0; ff_d <= f_d; ++ff_d) {
      Gudhi::Clock t;
      for (auto f_id: faces[f_d])
        for (auto ff_id: cs.face_range(f_id, ff_d)) {}
      t.end();
      total1[f_d][ff_d] += t.num_seconds() / faces[f_d].size() * 1000;
    }
  }

  std::cout << "\\begin{table}\n\\begin{tabular}{";
  for (unsigned i = 0; i <= std::min(dimension, (unsigned)8)+4; ++i)
    std::cout << "|l";
  std::cout << "|}\n\\hline\n\\multicolumn{3}{|c|}{Face dimension} ";
  for (unsigned i = 0; i <= std::min(dimension, (unsigned)8); ++i)
    std::cout << "& " << i << " ";
  std::cout << "\\\\ \\hline\n\\parbox[t]{2mm}{\\multirow{" << 2*(dimension + 1) << "}{*}{\\rotatebox[origin=c]{90}{\\centering Simplex dimension}}}";
  for (unsigned f_d = 0; f_d <= dimension; ++f_d) {
    // std::cout << "& \\multirow{2}{*}{" << f_d << "} & old ";
    std::cout << "& \\multirow{2}{*}{" << f_d << "} & first ";
    for (unsigned ff_d = 0; ff_d <= std::min(f_d, (unsigned)8); ++ff_d) {
      std::cout << "& ";
      if (std::round(total1[f_d][ff_d]*1000)/1000 > std::round(total2[f_d][ff_d]*1000)/1000)
        std::cout << "\\cellcolor{green!20} ";
      std::cout  << std::round(total1[f_d][ff_d]*1000)/1000. << " ";
    }
    for (unsigned ff_d = f_d + 1; ff_d <= std::min(dimension, (unsigned)8); ++ff_d) {
      std::cout << "& ";
    }
    // std::cout << "\\\\ \\cline{3-" << f_d + 4 << "}\n & & new "; 
    std::cout << "\\\\ \\cline{3-" << std::min(f_d, (unsigned)8) + 4 << "}\n & & second "; 
    for (unsigned ff_d = 0; ff_d <= std::min(f_d, (unsigned)8); ++ff_d) {
      std::cout << "& \\bf ";
      if (std::round(total1[f_d][ff_d]*1000)/1000 > std::round(total2[f_d][ff_d]*1000)/1000)
        std::cout << "\\cellcolor{green!20} ";
      std::cout << std::round(total2[f_d][ff_d]*1000)/1000. << " ";
    }
    for (unsigned ff_d = f_d + 1; ff_d <= std::min(dimension, (unsigned)8); ++ff_d) {
      std::cout << "& ";
    }
    std::cout << "\\\\ \\cline{2-" << std::min(dimension, (unsigned)8) + 4 << "}\n"; 
  }
  std::cout << "\\hline\n\\end{tabular}\n\\caption{The table of average running times in milliseconds of the old and the new algorithms to compute all faces of the simplices of various dimensions in a $" << dimension << "$-dimensional triangulation.}\n\\label{tab:compar-" << dimension << "}\n\\end{table}\n";

  // std::vector<double> total1(dimension + 1, 0), total2(dimension + 1, 0);
  // for (unsigned f_d = 0; f_d <= cs.dimension(); ++f_d) {
  //   Gudhi::Clock t;
  //   for (unsigned i = 0; i < num_tests; ++i) {
  //     A_id a_id = cs.query_point_location(*rp++, 1);
  //     for (auto f_id: cs.face_range(a_id, f_d)) {}
  //   }
  //   t.end();
  //   total1[f_d] += t.num_seconds() / num_tests * 1000;
  // }
  // for (unsigned f_d = 0; f_d <= cs.dimension(); ++f_d) {
  //   Gudhi::Clock t;
  //   for (unsigned i = 0; i < num_tests; ++i) {
  //     A_id a_id = cs.query_point_location(*rp++, 1);
  //     for (auto f_id: cs.face2_range(a_id, f_d)) {}
  //   }
  //   t.end();
  //   total2[f_d] += t.num_seconds() / num_tests * 1000;
  // }

  // std::cout << "\\begin{table}\n\\begin{tabular}{";
  // for (unsigned i = 0; i <= dimension+1; ++i)
  //   std::cout << "|l";
  // std::cout << "|}\n";
  // std::cout << "\\hline\n Face dimension ";
  // for (unsigned i = 0; i <= dimension; ++i)
  //   std::cout << "& " << i << " ";
  // std::cout << "\\\\\n";
  // std::cout << "\\hline\\hline\nOld algorithm\n";
  // for (auto t: total1)
  //   std::cout << "& " << t << " ";
  // std::cout << "\\\\\n\\hline\n";
  // std::cout << "New algorithm\n";
  // for (auto t: total2)
  //   std::cout << "& " << t << " ";
  // std::cout << "\\\\\n\\hline\n";
  // std::cout << "\\end{tabular}\n\\caption{The table of average running times in milliseconds of the two algorithms to compute all faces of the full-dimensional simplex in a $" << dimension << "$-dimensional triangulation.}\n\\label{tab:compar-" << dimension << "}\\end{table}\n";


  // std::vector<unsigned> decomposition(1, 0); // first coordinate is the sum
  // decomposition.reserve(dimension);
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
