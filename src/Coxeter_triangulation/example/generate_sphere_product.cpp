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
                     , std::vector<int> & dimensions
                     , std::string & out_file_name
                     ) {
  namespace po = boost::program_options;

  po::options_description visible("Allowed options", 100);
  visible.add_options()
    ("help,h", "produce help message")
    ("number,n", po::value<std::size_t>(&number_of_points)->default_value(0),
       "Number of generated points.")
    ("dimensions,D", po::value<std::vector<int>>(&dimensions)->multitoken(),
     "The list of dimensions of spheres.")
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
    std::cout << "Compute n random points uniformly on a product of spheres. \n";
    std::cout << std::endl << std::endl;

    std::cout << "Usage: " << argv[0] << " [options]" << std::endl << std::endl;
    std::cout << visible << std::endl;
    std::abort();
  }
}

int main(int argc, char * const argv[]) {
  
  std::string out_name = "default";
  std::size_t number_of_points;
  std::vector<int> dimensions = {1};

  program_options(argc, argv, number_of_points, dimensions, out_name);
  
  std::ofstream ofs(out_name, std::ofstream::out);
  ofs << "OFF\n" << number_of_points << " 0 0\n";
  
  // Generate uniformly random points on the (dim-1)-sphere
  for (std::size_t i = 1; i <= number_of_points; i++) {
    for (int d: dimensions) {
      CGAL::Random_points_on_sphere_d<Point_d> rp(d+1, 1);
      for (auto p_it = rp->cartesian_begin(); p_it != rp->cartesian_end(); ++p_it)
        ofs << *p_it << " ";
    }
    // ofs << "0.01 ";
    ofs << "\n";
  }
  ofs.close();

  std::cout << "Successfully generated " << number_of_points << " points.\n";
}
