#include <iostream>
#include <fstream>
#include <gudhi/Simplex_tree.h>
#include <gudhi/distance_functions.h>
#include <gudhi/choose_n_farthest_points.h>

#define PRINT_ZIGZAG_FILTRATION

// Types definition
using Simplex_tree      = 
                Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_zigzag_persistence>;
using Zz_edge           = Zigzag_edge<Simplex_tree>;
using Filtration_value  = Simplex_tree::Filtration_value;

using Point_d = std::vector<Filtration_value>;


int main(int argc, char* argv[])
{
  Filtration_value nu = 2.;
  Filtration_value mu = 3.;
  int dim_max         = 3;

//sequence of insertion and deletions of vertices and edges
  std::vector< Zz_edge >        edge_filtration;
//epsilon_i values, size() == #points
  std::vector<Filtration_value> filtration_values;
//extract points from file
  std::vector<Point_d> point_cloud;
  point_cloud.push_back({0.,0.});
  point_cloud.push_back({1.,0.});
  point_cloud.push_back({2.,1.});
  point_cloud.push_back({0.25,1.});

  std::cout << "-----Point cloud, original: \n";
  for(auto point : point_cloud) {
    for(auto x : point) { std::cout << x << " "; }
    std::cout << std::endl;
  }
  std::cout << std::endl; 

  //Compute edge filtration, after ordering points by furthest point ordering.
	zigzag_filtration_one_skeleton( point_cloud, 
                                  Gudhi::Euclidean_distance(), 
                                  nu, 
                                  mu, 
                                  edge_filtration,
                                  farthest_point_ordering() );

  //Print the zigzag filtration of the 1-skeleton
  std::cout << "-----Edge filtration: \n";
  for(auto edg : edge_filtration) 
  { //insertion or deletion
    if(edg.type()) { std::cout << "+ "; } else { std::cout << "- "; }
    //check whether this is a vertex or an edge ( u() == v() )
    if(edg.u() == edg.v()) { std::cout << edg.u() << " "; }//vertex
    else { std::cout <<  edg.u() << " " << edg.v() << " "; }//edge
    std::cout << "  [" << edg.fil() << "]" << std::endl;//filtration value
  }
  std::cout << std::endl;
 
  // traverse the entire oscillating Rips zigzag filtration 
  Simplex_tree st;
  //initialize the zigzag filtration ; this is mandatory
  st.initialize_filtration(edge_filtration, dim_max); 
  //access the zigzag filtration range of simplices
  auto zz_rg = std::move(st.filtration_simplex_range());
  
  size_t num_arrows        = 0;//total number of insertion and deletion of simplices
  size_t max_size_complex  = 0;//max size of a complex in the filtration
  size_t curr_size_complex = 0;//size of the current complex
  size_t curr_morse_complex = 0;//number of critical faces in current complex
  size_t max_morse_complex  = 0;//max number of critical faces in a complex

#ifdef PRINT_ZIGZAG_FILTRATION
  std::cout << "Simplex filtration: \n";
  std::cout << " ins/del  simplex  [fil_value , key]    (size cpx, max size cpx)\n";
#endif 

  for(auto it = zz_rg.begin(); it != zz_rg.end(); ++it ) {
    ++num_arrows;
    if(it.arrow_direction()) {//insertion

#ifdef PRINT_ZIGZAG_FILTRATION
      std::cout << "+ ";
#endif
      if(++curr_size_complex > max_size_complex) 
      {  max_size_complex = curr_size_complex;  }
      if(st.critical(*it)) {
        if(++curr_morse_complex > max_morse_complex) 
        { max_morse_complex = curr_morse_complex; }
      }
    }
    else { 
#ifdef PRINT_ZIGZAG_FILTRATION
      std::cout << "- "; 
#endif
      --curr_size_complex; 

      if(it.break_morse_pair()) {
        if(++curr_morse_complex > max_morse_complex) 
        { max_morse_complex = curr_morse_complex; }
      }
      else{ --curr_morse_complex; }

    }//deletion
#ifdef PRINT_ZIGZAG_FILTRATION
    // print list of vertices
    for(auto u : st.simplex_vertex_range(*it)) { std::cout << u << " "; }
    std::cout << "  [" << sqrt(st.filtration(*it)) << "," << st.key(*it) << "]";
    std::cout << "     (" << curr_size_complex << "," << max_size_complex << ")  --- ";

    if(it.arrow_direction()) {//forward arrow
      // std::cout << "  ";
      if(st.critical(*it)) { std::cout << "critical"; }
      else {
        auto p_sh = st.paired_with(*it);
        std::cout << "paired with:   ";
        for(auto u : st.simplex_vertex_range(p_sh)) { std::cout << u << " "; }
        std::cout << "  of key " << st.key(p_sh);
      }
    }
    else {//backward arrow
      if(it.break_morse_pair()) {
        auto p_sh = st.paired_with(*it);
        std::cout << "was paired with:   ";
        for(auto u : st.simplex_vertex_range(p_sh)) { std::cout << u << " "; }
        std::cout << "  of key " << st.key(p_sh) << "--> which becomes critical";
      
        it.make_critical();//it paired with a simplex that becomes critical
      }
      else {
        if(st.critical(*it)) { std::cout << "critical"; }
        else {
          auto p_sh = st.paired_with(*it);
          std::cout << "paired with:   ";
          for(auto u : st.simplex_vertex_range(p_sh)) { std::cout << u << " "; }
          std::cout << "  of key " << st.key(p_sh) << "--> which is also removed";
        }
      }
    }
    std::cout << std::endl;
#endif 
  }
#ifdef PRINT_ZIGZAG_FILTRATION
  std::cout << std::endl << std::endl;
#endif

  std::cout << std::endl << std::endl;
  std::cout << "Total number of arrows:        " << num_arrows << std::endl;
  std::cout << "Maximal size complexes:        " << max_size_complex << std::endl;
  std::cout << "Total number of 1-skel arrows: " << edge_filtration.size() << "\n";

  return 0;
}
 
