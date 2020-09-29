#include <iostream>
#include <fstream>
#include <chrono>
#include <gudhi/Simplex_tree.h>
#include "gudhi/reader_utils.h"
#include <gudhi/distance_functions.h>
#include <gudhi/Points_off_io.h>
#include <boost/program_options.hpp>
#include <CGAL/Epick_d.h>
// #include <gudhi/choose_n_farthest_points.h>
// #include <gudhi/pick_n_random_points.h>

// Types definition
using Simplex_tree      = 
                Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_zigzag_persistence>;
using Zz_edge           = Zigzag_edge<Simplex_tree>;
using Filtration_value  = Simplex_tree::Filtration_value;
using K                 = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
using Point_d           = typename K::Point_d;
using Points_off_reader = Gudhi::Points_off_reader<Point_d>;

void program_options( int argc, char* argv[]
                    , std::string& off_file_points
                    , Filtration_value& nu
                    , Filtration_value& mu
                    , int& dim_max);

int main(int argc, char* argv[])
{
  std::string off_file_points;
  Filtration_value nu, mu;
  int dim_max;

  program_options(argc, argv, off_file_points, nu, mu, dim_max);
//sequence of insertion and deletions of vertices and edges
  std::vector< Zz_edge >        edge_filtration;
//epsilon_i values, size() == #points
  std::vector<Filtration_value> filtration_values;
//kernel
  K k_d; 
//extract points from file
  Points_off_reader off_reader(off_file_points); //read points
  std::cout << "Point cloud of size "<< off_reader.get_point_cloud().size() << "\n";
//remove duplicate points
  off_reader.no_duplicate();
  std::cout << "Remove duplicates: point cloud of size " << 
                                        off_reader.get_point_cloud().size() << "\n";

  // traverse the entire oscillating Rips zigzag filtration 
  Simplex_tree st;
  //initialize the zigzag filtration ; this is mandatory. Use the squared Eucliean distance for efficiency. Note that we must nu*nu and mu*mu
  st.initialize_filtration(nu*nu, mu*mu, dim_max, off_reader.get_point_cloud(),
                           k_d.squared_distance_d_object(),
                           farthest_point_ordering());//sort points by furthest pt order 
  //access the zigzag filtration range of simplices
  auto zz_rg(st.filtration_simplex_range());
  
  size_t num_arrows        = 0;//total number of insertion and deletion of simplices
  size_t max_size_complex  = 0;//max size of a complex in the filtration
  size_t curr_size_complex = 0;//size of the current complex
  std::cout << "Simplex filtration: \n";
  std::cout << " ins/del  simplex  [fil_value , key]    (size cpx, max size cpx)\n";
  for(auto it = zz_rg.begin(); it != zz_rg.end(); ++it ) {
    ++num_arrows;
    if(it.arrow_direction()) {//insertion
      std::cout << "+ ";
      if(++curr_size_complex > max_size_complex) 
      {  max_size_complex = curr_size_complex;  }
    }
    else { std::cout << "- "; --curr_size_complex; }//deletion
    //print list of vertices
    for(auto u : st.simplex_vertex_range(*it)) { std::cout << u << " "; }
    std::cout << "  [" << sqrt(st.filtration(*it)) << "," << st.key(*it) << "]\n";
    std::cout << "  (" << curr_size_complex << "," << max_size_complex << ")\n";
  }
  std::cout << std::endl << std::endl;
  std::cout << "Total number of arrows:        " << num_arrows << std::endl;
  std::cout << "Maximal size complexes:        " << max_size_complex << std::endl;
  std::cout << "Total number of 1-skel arrows: " << edge_filtration.size() << "\n";

  return 0;
}
 
//program options
void program_options(int argc, char* argv[], std::string& off_file_points, Filtration_value& nu, Filtration_value &mu, int& dim_max) {
  namespace po = boost::program_options;
  po::options_description hidden("Hidden options");
  hidden.add_options()("input-file", po::value<std::string>(&off_file_points),
                       "Name of an OFF file containing a point set.\n");

  po::options_description visible("Allowed options", 100);
  visible.add_options()("help,h", "produce help message")
  ( "nu,n",
    po::value<Filtration_value>(&nu)->default_value(3.0),
    "Lower multiplicative factor in the oscillating Rips zigzag filtration.")
  ( "mu,m",
    po::value<Filtration_value>(&mu)->default_value(3.2),
    "Upper multiplicative factor in the oscillating Rips zigzag filtration.")
  ( "cpx-dimension,d", po::value<int>(&dim_max)->default_value(1),
    "Maximal dimension of the oscillating Rips complexes in the filtration.");

  po::positional_options_description pos;
  pos.add("input-file", 1);

  po::options_description all;
  all.add(visible).add(hidden);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(all).positional(pos).run(), vm);
  po::notify(vm);

  if (vm.count("help") || !vm.count("input-file")) {
    std::cout << std::endl;
    std::cout << "Compute the oscillating Rips zigzag filtration based on a point cloud, with Euclidean metric.\n\n";
    std::cout << "Usage: " << argv[0] << " [options] input-file" << std::endl << std::endl;
    std::cout << visible << std::endl;
    std::abort();
  }
}
