#include <gudhi/Points_off_io.h>
#include <gudhi/pick_n_random_points.h>

#include <CGAL/Epick_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Random.h>

#include <boost/program_options.hpp>

using Point = CGAL::Epick_d< CGAL::Dynamic_dimension_tag>::Point_d;
using Vector_of_points = std::vector<Point>;

const double PI  = 3.141592653589793238463;
#define _USE_MATH_DEFINES

class PointSetGen {
  public:
    void program_options(int argc, char * const argv[]
                         , double & steps
                         , double & end_thresold
                         , int 	  & repetetions
                         , char   & manifold
                         , int 	  & dimension
                         , int    & dim_max
                         , std::string & in_file_name
                         , std::string & out_file_name
                         ) {
      namespace po = boost::program_options;

      po::options_description visible("Allowed options", 100);
      visible.add_options()
        ("help,h", "produce help message")
      	("steps,s", po::value<double>(&steps)->default_value(0.1),
       		"Steps of the threshold")
        ("end_thresold,e", po::value<double>(&end_thresold)->default_value(1),
      		"Final threshold for rips complex.")
          
        ("repetetions,r", po::value<int>(&repetetions)->default_value(1),
        	"Num of repetetions of the experiments.")
      	("manifold,m", po::value<char>(&manifold)->default_value('s'),
       		"Type of manifold")
           
        ("dimensions,D", po::value<int>(&dimension)->default_value(2),
         "Dimension of the manifold.")

         ("dim_max,k ", po::value<int>(&dim_max)->default_value(2),
         "Maximum allowed dimension of the Rips complex.")

        ("input_file_name,i", po::value<std::string>(&in_file_name),
         "The input file.")
        ("out_file_name,o", po::value<std::string>(&out_file_name),
         "The output file.");

      po::options_description all;
      all.add(visible);

      po::variables_map vm;
      po::store(po::command_line_parser(argc, argv).
                options(all).run(), vm);
      po::notify(vm);

      if (vm.count("help")) {
        std::cout << std::endl;
        std::cout << "Computes rips complexes of different threshold values, to 'end_thresold', with priodic steps of 'steps' from a n random uniform point_vector on a selected manifold, . \n";
        std::cout << "Strongly collapses all the rips complexes and output the results in out_file. \n";
        std::cout << "The experiments are repeted 'repete' num of times for each threshold value. \n";
        std::cout << "type -m for manifold options, 's' for uni sphere, 'b' for unit ball, 'f' for file. \n";
        std::cout << "type -i 'filename' for Input file option for exported point sample. \n";
        std::cout << std::endl << std::endl;

        std::cout << "Usage: " << argv[0] << " [options]" << std::endl << std::endl;
        std::cout << visible << std::endl;
        std::abort();
      }
    }
};