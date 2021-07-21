#include <iostream>
#include <fstream>
#include <chrono>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Zigzag_persistence.h>
#include "gudhi/reader_utils.h"
#include <gudhi/distance_functions.h>
#include <gudhi/Points_off_io.h>
#include <boost/program_options.hpp>
#include <CGAL/Epick_d.h>
#include <gudhi/choose_n_farthest_points.h>
#include <gudhi/pick_n_random_points.h>

// Types definition
using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_zigzag_persistence>;
using Zz_edge = Gudhi::Zigzag_edge<Simplex_tree>;
using Filtration_value = Simplex_tree::Filtration_value;

using Zz_persistence = Gudhi::zigzag_persistence::Zigzag_persistence<Simplex_tree>;
using Zz_persistence_fast = Gudhi::zigzag_persistence::Zigzag_persistence<Simplex_tree,Gudhi::zigzag_persistence::Zigzag_persistence_fast>;

using K = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
using Point_d = typename K::Point_d;
using Points_off_reader = Gudhi::Points_off_reader<Point_d>;

void program_options( int argc, char* argv[]
                    , std::string& off_file_points
                    , std::string& output_diagram
                    , Filtration_value& nu
                    , Filtration_value& mu
                    , int& dim_max);

#define _VERBATIM_ 0

int main(int argc, char* argv[])
{
  std::string off_file_points, output_diagram;
  Filtration_value nu, mu;
  int dim_max;

  program_options(argc, argv, off_file_points, output_diagram, nu, mu, dim_max);

  std::chrono::time_point<std::chrono::high_resolution_clock> end_mem1, start_mem1;
  std::chrono::time_point<std::chrono::high_resolution_clock> end_mem2, start_mem2;

//sequence of insertion and deletions of vertices and edges
//epsilon_i values, size() == #points
  std::vector<Filtration_value> filtration_values;
//CGAL geometry kernel
  K k_d; 
//extract points from file
  Points_off_reader off_reader(off_file_points); //read points
  std::cout << "Point cloud of size "<< off_reader.get_point_cloud().size() << "\n";
//remove duplicate points
  off_reader.no_duplicate();
  // std::cout << "!!!!Should remove duplicates!!!!\n";
  std::cout << "Remove duplicates: point cloud of size " << 
                                        off_reader.get_point_cloud().size() << "\n";
{
  std::cout << "\n***** ZIGZAG PERSISTENT HOMOLOGY (columns as lists) *****\n";
  // traverse the entire oscillating Rips zigzag filtration 
  Simplex_tree st;
  //initialize the zigzag filtration ; this is mandatory. Use the squared Euclidean distance for efficiency. Note that we must use nu*nu and mu*mu. Reorder point by
  //farthest point ordering.
  auto start_init = std::chrono::high_resolution_clock::now();
  //  
  st.initialize_filtration( nu//*nu //use Euclidean distance squared
                          , mu//*mu //use Euclidean distance squared
                          , dim_max
            , off_reader.get_point_cloud()
            , k_d.squared_distance_d_object()//Euclidean distance squared
            , Gudhi::farthest_point_ordering()//sort points by furthest pt order 
            , Gudhi::sqrt_filtration<Simplex_tree>());//apply sqrt to all edg length
//
  auto end_init = std::chrono::high_resolution_clock::now();

  std::cout << "Initialize filtration in " << std::chrono::duration_cast<std::chrono::milliseconds>(end_init-start_init).count()/(double)1000 << " s\n";
  std::cout << "   (warning: this includes a farthest point sort of the data, which is slow)\n";

  std::cout << "Begin zigzag persistence:\n";
  auto start_zz = std::chrono::high_resolution_clock::now();
  Zz_persistence zz(st);
  zz.compute_zigzag_persistence();
  auto end_zz = std::chrono::high_resolution_clock::now();

  start_mem1 = start_zz;  end_mem1 = end_zz;

  // std::cout << "Persistence diagram (dionysus): \n";
  // zz.output_diagram_dionysus(std::cout);//, 0.0000001);
  
  // std::cout << "\n\n\n";
  // zz.output_diagram_log2(std::cout, 0.0000001);

  // std::cout << "Index persistence diagram:\n";
  // zz.output_index_diagram();

  std::cout << "\n";
  std::cout << "Compute zigzag persistence in: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_zz-start_zz).count()/(double)1000 << " s\n";
  std::cout << "\n";
}


{
  std::cout << "\n***** ZIGZAG PERSISTENT HOMOLOGY (columns as sets) *****\n";
  // traverse the entire oscillating Rips zigzag filtration 
  Simplex_tree st_opt;
  //initialize the zigzag filtration ; this is mandatory. Use the squared Euclidean distance for efficiency. Note that we must use nu*nu and mu*mu. Reorder point by
  //farthest point ordering.
  auto start_init_opt = std::chrono::high_resolution_clock::now();
  //  
  st_opt.initialize_filtration( nu//*nu //use Euclidean distance squared
                          , mu//*mu //use Euclidean distance squared
                          , dim_max
            , off_reader.get_point_cloud()
            , k_d.squared_distance_d_object()//Euclidean distance squared
            , Gudhi::farthest_point_ordering()//sort points by furthest pt order 
            , Gudhi::sqrt_filtration<Simplex_tree>());//apply sqrt to all edg length
  //
  auto end_init_opt = std::chrono::high_resolution_clock::now();

  std::cout << "Initialize filtration in " << std::chrono::duration_cast<std::chrono::milliseconds>(end_init_opt-start_init_opt).count()/(double)1000 << " s\n";
  std::cout << "   (warning: this includes a farthest point sort of the data, which is slow)\n";

  std::cout << "Begin zigzag persistence:\n";
  auto start_zz_opt = std::chrono::high_resolution_clock::now();
  Zz_persistence_fast zz_opt(st_opt);
  zz_opt.compute_zigzag_persistence();
  auto end_zz_opt = std::chrono::high_resolution_clock::now();

  std::cout << "Persistence diagram (dionysus): \n";
  // zz.output_diagram_dionysus(std::cout);//, 0.0000001);
  
  // std::cout << "\n\n\n";
  // zz.output_diagram_log2(std::cout, 0.0000001);

  std::cout << "Index persistence diagram:\n";
  // zz.output_index_diagram();

  std::cout << "\n";
  std::cout << "Compute optimized zigzag persistence in: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_zz_opt-start_zz_opt).count()/(double)1000 << " s\n";

  start_mem2 = start_zz_opt; end_mem2 = end_zz_opt;

  // std::cout << "Forward arrow: " << timing_forward_arrows_ << "\n";
  // std::cout << "Backward arrow: " << timing_backward_arrows_ << "\n";
  std::cout << "\n";

}







{
  long num_vertices = 0;
  long num_edges = 0;
  long num_insertions = 0;
  long num_removals = 0;

  num_vertices = off_reader.get_point_cloud().size();

  //time spent constructing the complex on the fly
  Simplex_tree st_timing;
  // st_timing.initialize_filtration(edge_filtration, dim_max); 
  st_timing.initialize_filtration( nu//*nu //use Euclidean distance squared
                                 , mu//*mu //use Euclidean distance squared
                                 , dim_max
            , off_reader.get_point_cloud()
            , k_d.squared_distance_d_object()//Euclidean distance squared
            , Gudhi::farthest_point_ordering()//sort points by furthest pt order 
            , Gudhi::sqrt_filtration<Simplex_tree>());//apply sqrt to all edg length
//

  if(_VERBATIM_) {
    std::cout << "Unordered list of data points, no duplicate: x1 ... xd\n";
    for(auto p : off_reader.get_point_cloud()) {
      std::cout << p << std::endl;
    }
    std::cout << std::endl;
  }

  auto start_stonly = std::chrono::high_resolution_clock::now();

  auto zzrg = st_timing.filtration_simplex_range();
  auto zzit = zzrg.begin();

  if(_VERBATIM_) {
    std::cout << "List of insertions and deletions of edges of the flag filtration (points ordered with farthest point): fil_value <-> u v \n";
    for(auto e : zzit.zigzag_edge_filtration_) {
      std::cout << " " << e << "\n";
    }
    std::cout << std::endl;
  }

  num_edges = zzit.zigzag_edge_filtration_.size();

  if(_VERBATIM_) { std::cout << "The whole zigzag filtration:\n"; }

  size_t num_arrow    = 0;//total number of insertions and deletions
  size_t max_complex  = 0;//maximal size of a complex
  size_t size_complex = 0;//current size of the complex
  while( zzit != zzrg.end() )
  { //new arrow

    //print filtration
    if(_VERBATIM_) {
      auto curr_simp = *zzit;
      // std::cout << st_timing.filtration(curr_simp) << " ";
      std::cout << zzit.filtration() << "   ";
      if(zzit.arrow_direction()) { std::cout << "-> "; }
      else { std::cout << "<- "; }
      for(auto v : st_timing.simplex_vertex_range(curr_simp)) {
        std::cout << v << " ";
      }
      std::cout << std::endl;
    }
    //end print filtration


    ++num_arrow;
    if(zzit.arrow_direction()) //insertion of a simplex
    { 
      ++num_insertions;
      ++size_complex; 
      if(size_complex > max_complex) { ++max_complex; }
    } 
    else { ++num_removals; --size_complex; } //removal of a simplex
    
    ++zzit;
  }
  auto end_stonly = std::chrono::high_resolution_clock::now();

  if(_VERBATIM_) { std::cout << "\n"; }


  std::cout << "Time spent on constructing the complex: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_stonly-start_stonly).count()/(double)1000 << " s\n";

  std::cout << "Total number of insertions and deletions: " << num_arrow << std::endl;

  std::cout << "Maximal size of a complex in the filtration: " << max_complex << std::endl;

    std::cout << "\n\n";

    std::cout << "Statistics:\n";
    std::cout << "Number of arrows:     " << num_removals + num_insertions << "\n";
    std::cout << "Number of insertions: " << num_insertions << "\n";
    std::cout << "Number of removals:   " << num_removals << "\n";
    std::cout << "Number of vertices:   " << num_vertices << "\n";
    // std::cout << "Number of edges:      " << num_edges << "\n";
    std::cout << "Max size complex:     " << max_complex << "\n";

    std::cout << "\n\n**********Compute zigzag persistence in: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_mem1-start_mem1).count()/(double)1000 << " s\n";

    std::cout << "      -> time per arrow = " << std::chrono::duration_cast<std::chrono::milliseconds>(end_mem1-start_mem1).count()/((double)1000 * (double)(num_removals + num_insertions)) << "\n";

    std::cout << "\n\n**********Compute optimized zigzag persistence in: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_mem2-start_mem2).count()/(double)1000 << " s\n";

    std::cout << "      -> time per arrow = " << std::chrono::duration_cast<std::chrono::milliseconds>(end_mem2-start_mem2).count()/((double)1000 * (double)(num_removals + num_insertions)) << "\n";

    std::cout << "\n";
    std::cout << " --- STANDARD / OPTIMIZED = " << (double)std::chrono::duration_cast<std::chrono::milliseconds>(end_mem1-start_mem1).count()/(double)std::chrono::duration_cast<std::chrono::milliseconds>(end_mem2-start_mem2).count() << "\n\n";

#ifdef _PROFILING_ZIGZAG_PERSISTENCE_
    std::cout << " ------- with:\n";
    std::cout << "   calls to forward_arrow:                     " << _num_forward_arrow_ << "\n";
    std::cout << "   calls to backward_arrow:                    " << _num_backward_arrow_ << "\n";
    std::cout << "   calls to arrow transposition case study:    " << _num_arrow_trans_case_study_ << "\n";
    std::cout << "     with in average arrow_trans/backward_arrow: " << (float)_num_arrow_trans_case_study_/(float)_num_backward_arrow_ << "\n";
    std::cout << "   calls to plus_equal_col in total:           " << _num_total_plus_equal_col_ << "\n";
    std::cout << "   calls to plus_equal_col within arrow_trans: " << _num_arrow_trans_case_study_plus_equal_col_ << "\n";
    std::cout << "     which represents (in %):                  " << (float)_num_arrow_trans_case_study_plus_equal_col_/(float)_num_total_plus_equal_col_ * 100. << "\n";
    std::cout << "\n";
    std::cout << "   average length manipulated column:          " << (double) _total_length_columns_/(double)_num_total_plus_equal_col_ << "\n";
#endif
    
}
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
  ( "nu-factor,n",
    po::value<Filtration_value>(&nu)->default_value(3.0),
    "Lower multiplicative factor in the oscillating Rips zigzag filtration.")
  ( "mu-factor,m",
    po::value<Filtration_value>(&mu)->default_value(3.2),
    "Upper multiplicative factor in the oscillating Rips zigzag filtration.")
  ( "cpx-dimension,d", po::value<int>(&dim_max)->default_value(1),
    "Maximal dimension of the oscillating Rips complexes in the filtration.")
  // ( "field-charac,p", po::value<int>(&p)->default_value(11),
  //   "Characteristic p of the coefficient field Z/pZ for computing homology.")
  // ( "min-persistence,m", po::value<Filtration_value>(&min_persistence),
  //   "Minimal lifetime of homology feature to be recorded. Default is 0. Enter a negative value to see zero length "
  //     "intervals")
  ;

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
    std::cout << "Usage: " << argv[0] << " [options] input-file-points output-file-persistence-diagram" << std::endl << std::endl;
    std::cout << visible << std::endl;
    std::abort();
  }
}