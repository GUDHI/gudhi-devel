#include <iostream>
#include <fstream>
#include <chrono>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Zigzag_persistence.h>
#include "gudhi/reader_utils.h"
#include <gudhi/distance_functions.h>
// #include <gudhi/Zigzag_filtration.h>
#include <gudhi/Points_off_io.h>
#include <boost/program_options.hpp>
#include <CGAL/Epick_d.h>
#include <gudhi/choose_n_farthest_points.h>
#include <gudhi/pick_n_random_points.h>

// Types definition
using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_zigzag_persistence>;
using Zz_persistence = Gudhi::zigzag_persistence::Zigzag_persistence<Simplex_tree> ;
using Zz_edge = Gudhi::Zigzag_edge<Simplex_tree>;
using Filtration_value = Simplex_tree::Filtration_value;
using K = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
using Point_d = typename K::Point_d;
using Points_off_reader = Gudhi::Points_off_reader<Point_d>;

void program_options( int argc, char* argv[]
                    , std::string& off_file_points
                    , Filtration_value& nu
                    , Filtration_value& mu
                    , int& dim_max);

// bool lex_cmp( Point_d p, Point_d q ) {
//   auto itp = p.begin(); auto itq = q.begin();
//   while(itp != p.end() && itq != q.end()) {
//     if(*itp != *itq) { return *itp < *itq; }
//     ++itp; ++itq;
//   }
//   return false;
// }

#define _VERBATIM_ 0

int main(int argc, char* argv[])
{
  // std::string off_file_points = "/Users/cmaria/Projects/gudhi-devel/data/points/tore3D_1307.off";
  // Filtration_value nu = 3.5;
  // Filtration_value mu = 3.7;
  // int dim_max = 3;
  std::string off_file_points = "/Users/cmaria/Projects/gudhi-devel/data/points/sphere3D_2646.off";
  Filtration_value nu = 3.2;
  Filtration_value mu = 3.2;
  int dim_max = 4;

  std::cout << off_file_points << " -n " << nu << " -m " << mu << " -d " << dim_max << "\n";

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

  std::cout << "Compute zigzag persistence in: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_zz-start_zz).count()/(double)1000 << " s\n";

#ifdef _PROFILING_ZIGZAG_PERSISTENCE_
    std::cout << " ------- with:\n";
    std::cout << "   calls to forward_arrow:                     " << _num_forward_arrow_ << "\n";
    std::cout << "   calls to backward_arrow:                    " << _num_backward_arrow_ << "\n";
    std::cout << "   calls to arrow transposition case study:    " << _num_arrow_trans_case_study_ << "\n";
    std::cout << "   calls to plus_equal_col in total:           " << _num_total_plus_equal_col_ << "\n";
    std::cout << "   calls to plus_equal_col within arrow_trans: " << _num_arrow_trans_case_study_plus_equal_col_ << "\n";
    std::cout << "     which represents (in %):                  " << (float)_num_arrow_trans_case_study_plus_equal_col_/(float)_num_total_plus_equal_col_ * 100. << "\n";
#endif

  long num_vertices = 0;
  long num_edges = 0;
  long num_insertions = 0;
  long num_removals = 0;

  num_vertices = off_reader.get_point_cloud().size();

  //time spent constructing the complex on the fly
  Simplex_tree st_timing;
  // st_timing.initialize_filtration(edge_filtration, dim_max); 
  st_timing.initialize_filtration( nu//*nu//*nu //use Euclidean distance squared
                                 , mu//*mu//*mu //use Euclidean distance squared
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

  return 0;
}
