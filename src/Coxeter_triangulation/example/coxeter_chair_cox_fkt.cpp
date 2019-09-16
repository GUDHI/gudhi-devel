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
#include <gudhi/Functions/Function_chair_in_R3.h>
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

using MT_CX = Manifold_tracing<Coxeter_triangulation<> >;
using MT_FR = Manifold_tracing<Freudenthal_triangulation<> >;

/** Test suite to find an empiric dependency on the intrinsic dimension. */ 
int main(int argc, char * const argv[]) {
  if (argc != 4) {
    std::cerr << "Usage: " << argv[0]
              << " low_d high_d level\n";
    return 0;
  }

  std::size_t low_d = atoi(argv[1]);
  std::size_t high_d = atoi(argv[2]);
  double level = atof(argv[3]);

  Gudhi::Clock t;
 
  std::vector<double> total1(high_d+1, 0);
  std::vector<std::size_t> size1(high_d+1, 0);
  std::vector<double> total2(high_d+1, 0);
  std::vector<std::size_t> size2(high_d+1, 0);
  for (std::size_t d = low_d; d <= high_d; d ++) {
    Function_chair_in_R3 sm;
    // std::cout << "sm(sm.seed) =\n" << sm(sm.seed()) << "\n";
    auto fun = make_embedding(sm, d);
    Point_d seed;
    fun.seed(seed);
    std::vector<Point_d> seed_points = {seed};
    // std::cout << "fun(fun.seed()) =\n" << fun(fun.seed()) << "\n\n";
    auto oracle = make_oracle(fun);

    Coxeter_triangulation<> cox_tr(d);
    Freudenthal_triangulation<> fr_tr(d);
    fr_tr.change_matrix(level * random_orthogonal_matrix(d) * fr_tr.matrix());
    fr_tr.change_offset(level*Eigen::VectorXd::Random(d));
    cox_tr.change_matrix(random_orthogonal_matrix(d) * cox_tr.matrix());
    cox_tr.change_offset(Eigen::VectorXd::Random(d));
    double diam = 0;
    auto s = cox_tr.locate_point(seed);
    for (auto v1: s.vertex_range())
      for (auto v2: s.vertex_range()) {
	Point_d p1 = cox_tr.cartesian_coordinates(v1);
	Point_d p2 = cox_tr.cartesian_coordinates(v2);
	double dist = (p1-p2).norm();
	if (dist > diam)
	  diam = dist;
      }
    cox_tr.change_matrix((sqrt(d) * level/diam) * cox_tr.matrix());

    {
      typename MT_CX::Out_simplex_map output;
      while (output.empty()) {
	t.begin();
	manifold_tracing_algorithm(seed_points, cox_tr, oracle, output);
	t.end();
	cox_tr.change_matrix(random_orthogonal_matrix(d) * cox_tr.matrix());
      }
      total1[d] = t.num_seconds();
      size1[d] = output.size();
    }
    {
      typename MT_FR::Out_simplex_map output;
      while (output.empty()) {
	t.begin();
	manifold_tracing_algorithm(seed_points, fr_tr, oracle, output);
	t.end();
	fr_tr.change_matrix(random_orthogonal_matrix(d) * fr_tr.matrix());
      }
      total2[d] = t.num_seconds();
      size2[d] = output.size();
    }
  }

  std::cout << "\\begin{table}\n\\centering\n\\begin{tabular}{|l|l";
  for (std::size_t d = low_d; d <= high_d; d++)
    std::cout << "|c";
  std::cout << "|}\n\\hline\n \\multicolumn{2}{|c|}{Ambient dimension}\n";
  for (std::size_t d = low_d; d <= high_d; d++)
    std::cout << "& " << d << " ";

  std::cout << "\\\\\n \\hline\\hline\n\\parbox[t]{2mm}{\\multirow{" << 2 << "}{*}{\\rotatebox[origin=c]{90}{\\centering CT}}}\n";
  std::cout << "& time, s ";
  for (std::size_t d = low_d; d <= high_d; d++)
    std::cout << "& " << total1[d] << " ";
  std::cout << "\\\\\n \\cline{2-" << high_d-low_d + 3 << "}\n";
  std::cout << "& size $|\\mathcal{S}|$ ";
  for (std::size_t d = low_d; d <= high_d; d++)
    std::cout << "& " << size1[d] << " ";
  std::cout << "\\\\\n \\cline{2-" << high_d-low_d + 3 << "}\n";
  std::cout << "& av. time, ms ";
  for (std::size_t d = low_d; d <= high_d; d++)
    std::cout << "& " << 1000*total1[d]/size1[d] << " ";
  
  std::cout << "\\\\\n \\hline\\hline\n\\parbox[t]{2mm}{\\multirow{" << 2 << "}{*}{\\rotatebox[origin=c]{90}{\\centering FKT}}}\n";
  std::cout << "& time, s ";
  for (std::size_t d = low_d; d <= high_d; d++)
    std::cout << "& " << total2[d] << " ";
  std::cout << "\\\\\n \\cline{2-" << high_d-low_d + 3 << "}\n";
  std::cout << "& size $|\\mathcal{S}|$ ";
  for (std::size_t d = low_d; d <= high_d; d++)
      std::cout << "& " << size2[d] << " ";
  std::cout << "\\\\\n \\cline{2-" << high_d-low_d + 3 << "}\n";
  std::cout << "& av. time, ms ";
  for (std::size_t d = low_d; d <= high_d; d++)
    std::cout << "& " << 1000*total2[d]/size2[d] << " ";
  
  std::cout << "\\\\ \\hline\n\\end{tabular}\n\\caption{The comparison of the performance of the manifold tracing algorithm using two types of the ambient triangulation: a Coxeter triangulation of type $\\tessell{A}{d}$ (CT) and the Freudenthal-Kuhn triangulation of $\\R^d$ (FKT) with the same diameter $" << level << "\\sqrt{d}$ of $d$-dimensional simplices.\n"
	    << "The reconstructed manifold is the two-dimensional implicit surface ``Chair'' embedded in $\\R^d$ given by the equation: $f(x) = (x_1^2+x_2^2+x_3^2 - 0.8)^2 - 0.4\\left((x_3-1)^2 - 2x_1^2\\right)\\left((x_3+1)^2 - 2x_2^2\\right)$.}\n\\label{tab:ct-fkt}\n\\end{table}\n";
  
}
