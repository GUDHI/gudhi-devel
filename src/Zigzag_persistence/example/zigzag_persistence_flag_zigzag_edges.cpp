/*    This file is a prototype for the Gudhi Library.
 *    Author(s):       Clément Maria
 *    Copyright (C) 2021 Inria
 *    This version is under developement, please do not redistribute this software. 
 *    This program is for academic research use only. 
 */

#include <iostream>
#include <fstream>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Zigzag_persistence.h>

// Types 
using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_zigzag_persistence>;
using Zz_persistence = Gudhi::zigzag_persistence::Zigzag_persistence<Simplex_tree> ;
using Zz_edge = Gudhi::Zigzag_edge<Simplex_tree>;
using Filtration_value = Simplex_tree::Filtration_value;

int main(int argc, char* argv[])
{
  int dim_max = 2;

  std::vector<Zz_edge> edge_filtration;
  //add 6 vertices 1 ... 6
  edge_filtration.emplace_back(1,1,0.1,true);
  edge_filtration.emplace_back(2,2,0.2,true);
  edge_filtration.emplace_back(3,3,0.3,true);
  edge_filtration.emplace_back(4,4,0.4,true);
  edge_filtration.emplace_back(5,5,0.5,true);
  edge_filtration.emplace_back(6,6,0.6,true);
  //make a hexagon 1-2...-6-1
  edge_filtration.emplace_back(1,2,1.1,true);
  edge_filtration.emplace_back(2,3,1.2,true);
  edge_filtration.emplace_back(3,4,1.3,true);
  edge_filtration.emplace_back(4,5,1.4,true);
  edge_filtration.emplace_back(5,6,1.5,true);
  edge_filtration.emplace_back(6,1,1.6,true);
  //add a centre 0 connected to every second vertex
  edge_filtration.emplace_back(0,0,1.65,true);
  edge_filtration.emplace_back(0,2,1.7,true);
  edge_filtration.emplace_back(0,4,1.8,true);
  edge_filtration.emplace_back(0,6,1.9,true);
  //add a new vertex 7, an connected it to 1 .. 6
  edge_filtration.emplace_back(7,7,2.0,true);
  edge_filtration.emplace_back(7,1,2.1,true);
  edge_filtration.emplace_back(7,2,2.2,true);
  edge_filtration.emplace_back(7,3,2.3,true);
  edge_filtration.emplace_back(7,4,2.4,true);
  edge_filtration.emplace_back(7,5,2.5,true);
  edge_filtration.emplace_back(7,6,2.6,true);
  //remove all edges adjacent to 7, and 7 itself
  edge_filtration.emplace_back(7,1,3.1,false);
  edge_filtration.emplace_back(7,2,3.2,false);
  edge_filtration.emplace_back(7,3,3.3,false);
  edge_filtration.emplace_back(7,4,3.4,false);
  edge_filtration.emplace_back(7,5,3.5,false);
  edge_filtration.emplace_back(7,6,3.6,false);
  edge_filtration.emplace_back(7,7,4.0,false);
  //remove the three central edges in order of apparition
  edge_filtration.emplace_back(0,2,4.1,false);
  edge_filtration.emplace_back(0,4,4.2,false);
  edge_filtration.emplace_back(0,6,4.3,false);  
  edge_filtration.emplace_back(0,0,4.4,false);

 // traverse the entire oscillating Rips zigzag filtration 
  Simplex_tree st;
  //initialize the zigzag filtration with the vertices and edges; this is mandatory. 
  st.initialize_filtration( edge_filtration, dim_max );
  //compute zigzag persistence
  Zz_persistence zz(st);
  zz.zigzag_persistent_homology();

  std::cout << "Persistence diagram (dionysus): \n";
  zz.persistence_diagram(std::cout, 0.5);

  //informations on the whole filtration
  Simplex_tree st_stats;
  st_stats.initialize_filtration( edge_filtration, dim_max );
  //zigzag filtration iterator
  auto zzrg = st_stats.filtration_simplex_range();
  auto zzit = zzrg.begin();

  std::cout << "List of insertions and deletions of edges of the flag filtration: [num_arrow] fil_value <-> u v \n";
  for(auto e : zzit.zigzag_edge_filtration_) {
    std::cout << " " << e << "\n";
  }
  std::cout << std::endl;

  std::cout << "The whole zigzag filtration:\n";

  size_t num_arrow    = 0;//total number of insertions and deletions
  size_t max_complex  = 0;//maximal size of a complex
  size_t size_complex = 0;//current size of the complex
  while( zzit != zzrg.end() )
  { //new arrow
    std::cout << "[" << num_arrow << "] ";
    auto curr_simp = *zzit;
    std::cout << zzit.filtration() << " ";
    if(zzit.arrow_direction()) { std::cout << "-> "; }
    else { std::cout << "<- "; }
    for(auto v : st_stats.simplex_vertex_range(curr_simp)) {
      std::cout << v << " ";
    }
    std::cout << std::endl;

    ++num_arrow;
    if(zzit.arrow_direction()) //insertion of a simplex
    { 
      ++size_complex; 
      if(size_complex > max_complex) { ++max_complex; }
    } 
    else { --size_complex; } //removal of a simplex
    
    ++zzit;
  }
  std::cout << "\n"; 
  std::cout << "Total number of insertions and deletions: " 
            << num_arrow << std::endl;
  std::cout << "Maximal size of a complex in the filtration: " 
            << max_complex << std::endl;

  return 0;
}
 