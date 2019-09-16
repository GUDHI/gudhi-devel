// workaround for the annoying boost message in boost 1.69
#define BOOST_PENDING_INTEGER_LOG2_HPP
#include <boost/integer/integer_log2.hpp>
// end workaround 

#include <iostream>
#include <vector>
#include <algorithm>

#include <gudhi/Points_off_io.h>
#include <gudhi/Permutahedral_representation.h>
#include <gudhi/Coxeter_triangulation.h>
#include <gudhi/Manifold_tracing.h>
#include <gudhi/Implicit_manifold_intersection_oracle.h>
#include <gudhi/Functions/Function_Sm_in_Rd.h>
#include <gudhi/Functions/Embed_in_Rd.h>
#include <gudhi/Functions/random_orthogonal_matrix.h>
#include <gudhi/Clock.h>

#include <boost/math/special_functions/binomial.hpp>

#include <CGAL/Epick_d.h>
#include <CGAL/point_generators_d.h>

using namespace Gudhi::coxeter_triangulation;
using K = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
using FT = K::FT;
using Point_d = Eigen::VectorXd;
using Point_vector = std::vector< Point_d >;

using MT = Manifold_tracing<Coxeter_triangulation<> >;

/** Test suite to find an empiric dependency on the intrinsic dimension. */ 
int main(int argc, char * const argv[]) {
  if (argc != 7) {
    std::cerr << "Usage: " << argv[0]
              << " low_d high_d low_m high_m step level\n";
    return 0;
  }

  std::size_t low_d = atoi(argv[1]);
  std::size_t high_d = atoi(argv[2]);
  std::size_t low_m = atoi(argv[3]);
  std::size_t high_m = atoi(argv[4]);
  std::size_t step = atoi(argv[5]);
  double level = atof(argv[6]);

  MT mt;
  Gudhi::Clock t;
 
  std::vector<std::vector<double> > total1(high_d+1,
					   std::vector<double>(high_d+1, 0));
  std::vector<std::vector<double> > size1(high_d+1,
					  std::vector<double>(high_d+1, 0));
  for (std::size_t d = low_d; d <= high_d; d += step) {
    Coxeter_triangulation<> cox_tr(d);
    cox_tr.change_matrix(level * random_orthogonal_matrix(d) * cox_tr.matrix());
    cox_tr.change_offset(level*Eigen::VectorXd::Random(d));

    for (std::size_t m = low_m; m <= std::min(high_m, d-1); ++m) {
      Function_Sm_in_Rd sm(2.1, m);
      // std::cout << "sm(sm.seed) =\n" << sm(sm.seed()) << "\n";
      auto fun = make_embedding(sm, d);
      Point_d seed;
      fun.seed(seed);
      std::vector<Point_d> seed_points = {seed};
      // std::cout << "fun(fun.seed()) =\n" << fun(fun.seed()) << "\n\n";
      auto oracle = make_oracle(fun);
      typename MT::Out_simplex_map output;
      while (output.empty()) {
	t.begin();
        manifold_tracing_algorithm(seed_points, cox_tr, oracle, output);
	t.end();
	cox_tr.change_matrix(random_orthogonal_matrix(d) * cox_tr.matrix());
      }
      total1[d][m] = t.num_seconds();
      size1[d][m] = output.size();
    }
  }

  std::cout << "\\begin{table}\n\\centering\n\\begin{tabular}{|l|l";
  for (std::size_t i = low_m; i <= high_m; i += step)
    std::cout << "|c";
  std::cout << "|}\n\\hline\n \\multicolumn{2}{|c|}{\\backslashbox{d}{m}}\n";
  for (std::size_t i = low_m; i <= high_m; i += step)
    std::cout << "& " << i << " ";

  for (std::size_t d = low_d; d <= high_d; ++d) {
    std::cout << "\\\\\n \\hline\\hline\n\\multirow{" << 2 << "}{*}{" << d << "}\n";
    std::cout << "& time, s ";
    std::size_t m = low_m;
    for (; m < d && m <= high_m; m += step)
      std::cout << "& " << total1[d][m] << " ";
    for (; m <= high_m; m += step)
      std::cout << "& \\multirow{2}{*}{}\n";
    std::size_t col1 = 3;
    std::size_t col2 = std::floor((high_m-low_m+1)/(double)step) + 2;
    std::cout << "\\\\\n \\cline{2-2} \\cline{" << col1 << "-" << col2 << "}\n";
    std::cout << "& $|\\mathcal{S}|$ ";

    for (m = low_m; m < d && m <= high_m; m += step)
      std::cout << "& " << size1[d][m] << " ";
    for (; m <= high_m; m += step)
      std::cout << "& \\multirow{2}{*}{}\n";
  }
  std::cout << "\\\\ \\hline\n\\end{tabular}\n\\caption{The effect of the intrinsic dimension $m$ and the ambient dimension $d$ on the execution time and the output size of the manifold tracing algorithm.\nThe reconstructed manifold in the tests is the $m$-dimensional sphere embedded in $\\R^d$.\nThe ambient triangulation used in the reconstruction is a Coxeter triangulation of type $\\tessell{A}{d}$.\nThe time shown in the table is the computation time of the set $\\mathcal{S}$ in milliseconds.}\n\\end{table}\n";
  
}
