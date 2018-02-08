#include <iostream>
#include <vector>

#include <gudhi/Points_off_io.h>
#include <gudhi/Coxeter_system.h>
#include <gudhi/Coxeter_complex.h>
#include <gudhi/Coxeter_complex/Off_point_range.h>

#include <CGAL/Epick_d.h>

//#include <Eigen/Dense>

#include "memory_usage.h"
#include "cxx-prettyprint/prettyprint.hpp"

using K = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
using FT = K::FT;
using Point_d = K::Point_d;
using Point_vector = std::vector< Point_d >;
using Coxeter_complex = Gudhi::Coxeter_complex<Point_vector, Coxeter_system>;

std::vector<FT> bounding_box_dimensions(Point_vector& points) {
  std::vector<FT> lower, upper, difference;
  for (auto x: points[0]) {
    lower.push_back(x);
    upper.push_back(x);
  }
  for (auto p: points)
    for (unsigned i = 0; i < p.size(); i++) {
      if (p[i] < lower[i])
        lower[i] = p[i];
      if (p[i] > upper[i])
        upper[i] = p[i];
    }
  for (unsigned i = 0; i < lower.size(); i++)
    difference.push_back(upper[i]-lower[i]);
  return difference;
}


/** Current state of the algorithm.
 *  Input: a point cloud 'point_vector'
 *  Output: a reconstruction (a simplicial complex?, a Czech-like complex?)
 */

int main(int argc, char * const argv[]) {
  if (argc > 4 || argc < 2) {
    std::cerr << "Usage: " << argv[0]
        << " path_to_off_point_file [initial level] [epsilon]\n";
    return 0;
  }
  std::cout << "Coxeter complex computation for " << argv[1] << ", lev = " << argv[2] << ", eps = " << argv[3] << ".\n";
  
  double init_level = 1, eps = 0;
  bool store_in_ram = false;
  if (argc >= 3)
    init_level = atof(argv[2]);
  if (argc == 4)
    eps = atof(argv[3]);
  int d = 0;
  if (store_in_ram) {
    Gudhi::Points_off_reader<Point_d> off_reader(argv[1]);
    if (!off_reader.is_valid()) {
      std::cerr << "Coxeter triangulations - Unable to read file " << argv[1] << "\n";
      exit(-1);
    }
    Point_vector* point_vector = new Point_vector(off_reader.get_point_cloud());
    int N = point_vector->size();
    d = (*point_vector)[0].size();
    std::cout << "Successfully read " << N << " points in dimension " << d << std::endl;
    Coxeter_system cs_A('A', d);
    Coxeter_complex cc(*point_vector, cs_A, init_level, eps);
    delete point_vector;
    cc.write_mesh("sphere_coxeter_complex_A.mesh");
    cc.collapse();
    // cc.construct_clique_complex();
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
    // cc.construct_clique_complex();
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
