#include <sys/types.h>
#include <sys/stat.h>

#include <gudhi/Points_off_io.h>
#include <gudhi/choose_n_farthest_points.h>
#include "../example/generators.h"

#include <CGAL/Epick_d.h>

#include <boost/program_options.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> Kernel;
typedef typename Kernel::Point_d Point_d;
typedef std::vector<Point_d> Point_range;

void program_options(int argc, char * const argv[]
                     , std::size_t & number_of_points
                     , double & in_radius
                     , std::string & out_file_name
                     ) {
  namespace po = boost::program_options;

  po::options_description visible("Allowed options", 100);
  visible.add_options()
    ("help,h", "produce help message")
    ("max-number,n", po::value<std::size_t>(&number_of_points)->default_value(0),
       "Maximal number of generated points.")
    ("in-radius,r", po::value<double>(&in_radius)->default_value(0.5),
       "The inner circle radius.")
    ("output-file,o", po::value<std::string>(&out_file_name),
     "The output file.");

  po::options_description all;
  all.add(visible);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).
            options(all).run(), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << std::endl;
    std::cout << "Compute n random uniform points on a torus embedded in R^3. \n";
    std::cout << std::endl << std::endl;

    std::cout << "Usage: " << argv[0] << " [options]" << std::endl << std::endl;
    std::cout << visible << std::endl;
    std::abort();
  }
}

int main(int argc, char * const argv[]) {
  
  std::string out_name = "default";
  std::size_t number_of_points;

  double r = 0.5, R = 1;
  double pi = 3.14159265358;
  
  program_options(argc, argv, number_of_points, r, out_name);
  
  std::ofstream ofs(out_name, std::ofstream::out);
  ofs << "OFF\n" << number_of_points << " 0 0\n";
  
  CGAL::Random rand;
  for (unsigned i = 0; i < number_of_points; ++i) {
    double x, y, f;
    do {
      x = rand.get_double(0, 2.*pi);
      y = rand.get_double(0, 1./pi);
      f = (1. + (r/R)*cos(x))/(2*pi);    
    }
    while (y >= f);
    // we just need the value of x
    y = rand.get_double(0, 2.*pi);
    double x_coord = (R + r*cos(x))*cos(y);
    double y_coord = (R + r*cos(x))*sin(y);
    double z_coord = r*sin(x);
    ofs << x_coord << " " << y_coord << " " << z_coord << "\n";
  }
  ofs.close();

  std::cout << "Successfully generated " << number_of_points << " points.\n";
}
