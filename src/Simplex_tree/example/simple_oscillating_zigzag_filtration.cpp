#include <iostream>
#include <fstream>
#include <gudhi/Simplex_tree.h>
#include <gudhi/distance_functions.h>
#include <gudhi/choose_n_farthest_points.h>

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
  int dim_max         = 10;

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
    std::cout << "  [" << st.filtration(*it) << "," << st.key(*it) << "]";
    std::cout << "  (" << curr_size_complex << "," << max_size_complex << ")\n";
  }
  std::cout << std::endl << std::endl;
  std::cout << "Total number of arrows:        " << num_arrows << std::endl;
  std::cout << "Maximal size complexes:        " << max_size_complex << std::endl;
  std::cout << "Total number of 1-skel arrows: " << edge_filtration.size() << "\n";

  return 0;
}
 
