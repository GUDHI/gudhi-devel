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
                     , std::string & out_file_name
                     ) {
  namespace po = boost::program_options;

  po::options_description visible("Allowed options", 100);
  visible.add_options()
    ("help,h", "produce help message")
    ("number,n", po::value<std::size_t>(&number_of_points)->default_value(0),
       "Number of generated points.")
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
    std::cout << "Compute n random uniform points on a product of spheres. \n";
    std::cout << std::endl << std::endl;

    std::cout << "Usage: " << argv[0] << " [options]" << std::endl << std::endl;
    std::cout << visible << std::endl;
    std::abort();
  }
}

int main(int argc, char * const argv[]) {
  
  std::string out_name = "default";
  std::size_t number_of_points;

  program_options(argc, argv, number_of_points, out_name);
  
  std::ofstream ofs(out_name, std::ofstream::out);
  ofs << "OFF\n" << number_of_points << " 0 0\n";
  
  // Generate uniformly random points on the (dim-1)-sphere
  for (unsigned i = 0; i < number_of_points; ++i) {
    CGAL::Random_points_on_sphere_d<Point_d> rp(3, 1);
    auto p_it = rp->cartesian_begin();
    double x = *p_it++, y = *p_it++, z = *p_it;
    ofs << x*y << " " << x*z << " " << y*y-z*z << " " << 2*y*z;
    ofs << "\n";
  }
  ofs.close();

  std::cout << "Successfully generated " << number_of_points << " points.\n";
}
