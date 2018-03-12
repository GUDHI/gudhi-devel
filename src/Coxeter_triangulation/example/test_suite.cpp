#include <iostream>
#include <vector>
#include <queue>

#include <gudhi/Points_off_io.h>
#include <gudhi/Coxeter_system.h>
#include <gudhi/Coxeter_complex.h>
#include <gudhi/Coxeter_complex/Off_point_range.h>

#include <CGAL/Epick_d.h>

#include <boost/program_options.hpp>

//#include <Eigen/Dense>

#include "memory_usage.h"
#include "cxx-prettyprint/prettyprint.hpp"


using K = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
using FT = K::FT;
using Point_d = K::Point_d;
using Point_vector = std::vector< Point_d >;
using Coxeter_complex = Gudhi::Coxeter_complex<Point_vector, Coxeter_system>;

void program_options(int argc, char * const argv[]
                     , std::string & in_file_name
                     , double & lambda_step
                     , double & gamma_step
                     , double & time_limit
                     , int & max_abs_euler
                     , std::vector<int> & betti
                     ) {
  namespace po = boost::program_options;

  po::options_description visible("Allowed options", 100);
  visible.add_options()
    ("help,h", "produce help message")
    ("input-file,i", po::value<std::string>(&in_file_name),
     "The input file.")
    ("lambda-step,l", po::value<double>(&lambda_step)->default_value(0.01),
       "The step size for the lambda parameter.")
    ("gamma-step,g", po::value<double>(&gamma_step)->default_value(0.001),
       "The step size for the gamma parameter.")
    ("time-limit,t", po::value<double>(&time_limit)->default_value(std::numeric_limits<double>::infinity()),
       "The time limit for the execution.")
    ("max-abs-euler,e", po::value<int>(&max_abs_euler)->default_value(100),
       "The maximal Euler characteristic in absolute value.")
    ("betti,B", po::value<std::vector<int>>(&betti)->multitoken(),
     "The list of desired Betti numbers.");

  po::options_description all;
  all.add(visible);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).
            options(all).run(), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << std::endl;
    std::cout << "Find the parameters that are the quickest for the right Betti numbers. \n";
    std::cout << std::endl << std::endl;

    std::cout << "Usage: " << argv[0] << " [options]" << std::endl << std::endl;
    std::cout << visible << std::endl;
    std::abort();
  }
}


/** Current state of the algorithm.
 *  Input: a point cloud 'point_vector'
 *  Output: a reconstruction (a simplicial complex?, a Czech-like complex?)
 */

int main(int argc, char * const argv[]) {
  double init_level = 1, eps = 0;
  double lambda_step = 0.01, gamma_step = 0.01, time_limit = std::numeric_limits<double>::infinity();
  int max_abs_euler = 100;
  bool store_in_ram = true;
  std::vector<int> betti;
  std::string in_file_name;
  // if (argc >= 3)
  //   lambda_step = atof(argv[2]);
  // if (argc == 4)
  //   gamma_step = atof(argv[3]);
  program_options(argc, argv, in_file_name, lambda_step, gamma_step, time_limit, max_abs_euler, betti);
  std::cout << "Coxeter complex computation for " << in_file_name << ", lambda_step = " << lambda_step << ", gamma_step = " << gamma_step << ", time_limit = " << time_limit << ".\n";
  
  int d = 0;
  using Test = std::tuple<double, int, int, double, std::size_t, std::size_t, std::size_t>;
  using Conf = std::tuple<double, int, int>;
  auto test_comp =
    [](const Conf& t1, const Conf& t2)
    {return std::get<0>(t1) > std::get<0>(t2);};
  std::priority_queue<Conf, std::vector<Conf>, decltype(test_comp)> min_queue(test_comp); 
  if (store_in_ram) {
    Gudhi::Points_off_reader<Point_d> off_reader(in_file_name);
    if (!off_reader.is_valid()) {
      std::cerr << "Coxeter triangulations - Unable to read file " << in_file_name << "\n";
      exit(-1);
    }
    Point_vector* point_vector = new Point_vector(off_reader.get_point_cloud());
    int N = point_vector->size();
    d = (*point_vector)[0].size();
    std::cout << "Successfully read " << N << " points in dimension " << d << std::endl;
    bool too_coarse = true, found = false;
    std::vector<int> betti_;
    double a_time_, col_time_;
    int chi_;
    std::size_t max_simplices_, vertices_, coll_max_simplices_;
    Test best;
    int max_line = 0;
    int k = 0;
    // First step: Find the first parameter couple (k,l) such that the Betti numbers are at least the desired ones
    while (too_coarse) {
      Coxeter_system cs_A('A', d);
      Coxeter_complex cc(*point_vector, cs_A, init_level + k*lambda_step, eps);
      cc.collapse();
      a_time_ = cc.a_time_;
      col_time_ = cc.col_time_;
      chi_ = cc.chi_;
      betti_ = cc.betti_;
      max_simplices_ = cc.max_simplices_;
      coll_max_simplices_ = cc.coll_max_simplices_;
      vertices_ = cc.vertices_;
      if (betti_ == betti && a_time_ < time_limit) {
        found = true;
        best = std::make_tuple(a_time_, k, 0, col_time_, max_simplices_, vertices_, coll_max_simplices_);
        time_limit = a_time_;
      }
      std::cout << "Parameters lambda=" << init_level + k*lambda_step << " and gamma=" << eps << " a_time=" << a_time_ << " and col_time=" << col_time_ << ". Max simplices=" << max_simplices_ << ", vertices=" << vertices_ <<", collapse max simplices=" << coll_max_simplices_ << ", Betti numbers = " << betti_ << ", chi=" << chi_ << ".\n";
      too_coarse = (betti_ < betti);
      k++;
    }
    if (a_time_ < time_limit)    
      min_queue.push(std::make_tuple(a_time_, k-1, 0));
    // min_queue.push(std::make_tuple(a_time_, k, 0));
    max_line = k-1;
    // Main loop: Each element in min_queue corresponds to some maximal l for the given k
    while (!min_queue.empty()) {
      int k = std::get<1>(min_queue.top());
      int l = std::get<2>(min_queue.top())+1;
      Coxeter_system cs_A('A', d);
      Coxeter_complex cc(*point_vector, cs_A, init_level + k*lambda_step, eps + l*gamma_step);
      cc.collapse();
      a_time_ = cc.a_time_;
      col_time_ = cc.col_time_;
      betti_ = cc.betti_;
      chi_ = cc.chi_;
      max_simplices_ = cc.max_simplices_;
      vertices_ = cc.vertices_;
      if (betti_ == betti && a_time_ < time_limit) {
        found = true;
        best = std::make_tuple(a_time_, k, l, col_time_, max_simplices_, vertices_, coll_max_simplices_);
        time_limit = a_time_;
      }      
      std::cout << "Parameters lambda=" << init_level + k*lambda_step << " and gamma=" << eps + l*gamma_step << " a_time=" << a_time_ << " and col_time=" << col_time_ << ". Max simplices=" << max_simplices_ << ", vertices=" << vertices_ << ", collapse max simplices=" << coll_max_simplices_ << ", Betti numbers = " << betti_ << ", chi=" << chi_ << ".\n";
      min_queue.pop();
      if (betti_ >= betti && a_time_ < time_limit && std::abs(chi_) < max_abs_euler)
        min_queue.push(std::make_tuple(a_time_, k, l));
      if (k == max_line) {
        Coxeter_system cs_A('A', d);
        Coxeter_complex cc(*point_vector, cs_A, init_level + (k+1)*lambda_step, eps);
        cc.collapse();
        a_time_ = cc.a_time_;
        col_time_ = cc.col_time_;
        betti_ = cc.betti_;
        chi_ = cc.chi_;
        max_simplices_ = cc.max_simplices_;
        vertices_ = cc.vertices_;
        if (betti_ == betti) {
          found = true;
          best = std::make_tuple(a_time_, k+1, 0, col_time_, max_simplices_, vertices_, coll_max_simplices_);
          time_limit = a_time_;
        }
        std::cout << "Parameters lambda=" << init_level + (k+1)*lambda_step << " and gamma=" << eps << " a_time=" << a_time_ << " and col_time=" << col_time_ << ". Max simplices=" << max_simplices_ << ", vertices=" << vertices_ << ", collapse max simplices=" << coll_max_simplices_ << ", Betti numbers = " << betti_ << ", chi=" << chi_ << ".\n";
        if (a_time_ < time_limit && std::abs(chi_) < max_abs_euler)
          min_queue.push(std::make_tuple(a_time_, k+1, 0));
        max_line = k+1;
      }
    }
    if (found)
      std::cout << "The best time is a_time=" << std::get<0>(best) << ", col_time=" << std::get<3>(best) << " for values lambda=" << init_level + std::get<1>(best)*lambda_step << " and gamma=" << eps + std::get<2>(best)*gamma_step << ". Max simplices=" << std::get<4>(best) << ", vertices=" << std::get<5>(best) << ", collapse max simplices=" << std::get<6>(best) << ".\n";
    else
      std::cout << "No value found within the time and Euler characteristic limit.\n";
  }
  else {
    Gudhi::Off_point_range<Point_d> off_range(argv[1]);
    d = off_range.dimension();
    std::cout << "Successfully opened the file of points in dimension " << d << std::endl;
    using Coxeter_complex_off = Gudhi::Coxeter_complex<Gudhi::Off_point_range<Point_d>, Coxeter_system>;
    
    Coxeter_system cs_A('A', d);
    Coxeter_complex_off cc(off_range, cs_A, init_level, eps);  
    cc.write_mesh("sphere_coxeter_complex_A.mesh");
    std::cout << "Memory usage (Physical) before collapses: " << (float)getPhysicalValue()/1000 << "MB.\n";
    cc.collapse();
  }    
  std::cout << "Memory usage (Virtual): " << (float)getVirtualValue()/1000. << "MB.\n";
  std::cout << "Memory usage (Physical): " << (float)getPhysicalValue()/1000 << "MB.\n";
  // {
  //   Coxeter_system cs_B('B', d);
  //   Coxeter_complex cc(point_vector, cs_B, init_level);
  //   cc.write_mesh("sphere_coxeter_complex_B.mesh");
  // }
  // {
  //   Coxeter_system cs_C('C', d);
  //   Coxeter_complex cc(point_vector, cs_C, init_level); 
  //   cc.write_mesh("sphere_coxeter_complex_C.mesh");
  // }
 // Coxeter_system cs_D('D', d);
  // Coxeter_complex(point_vector, cs_D);  
  // Coxeter_system cs_E6('E', 6);
  // cs_E6.emplace_back('A', d-6);
  // Coxeter_complex(point_vector, cs_E6);  

  
}
