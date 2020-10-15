/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <iostream>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Discrete_morse_theory.h>

//allows Morse matching
using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_zigzag_persistence>;
using Vertex_handle = Simplex_tree::Vertex_handle;
using Filtration_value = Simplex_tree::Filtration_value;
using Simplex_handle = Simplex_tree::Simplex_handle;

int main(int argc, char* argv[])
{
  Simplex_tree st;

  //loop 0-1-2-3-4-5-0, plus triangle 3-4-6
  for(int v = 0; v < 7; ++v) {//insert vertices
    st.insert_simplex({v},(Filtration_value)(v));
  } 
  //edges 0-1, 1-2, 2-3, 3-4
  for(int v = 0; v < 4; ++v) {
    st.insert_simplex({v,v+1},(Filtration_value)(7+v));
  }
  //edges 3-6 and 4-6, triangle 3-4-6
  st.insert_simplex({3,6},(Filtration_value)(11));
  st.insert_simplex({4,6},(Filtration_value)(12));
  st.insert_simplex({3,4,6},(Filtration_value)(13));
  //edges 4-5 and 5-0
  st.insert_simplex({4,5},(Filtration_value)(14));
  st.insert_simplex({0,5},(Filtration_value)(15));

  std::vector<Simplex_handle> sh_range; 
  for(auto sh : st.filtration_simplex_range()) {  sh_range.push_back(sh); }
//sort by filtration value. Read from end to begin, only max simplices
  sort(sh_range.begin(),sh_range.end(), 
        [&](Simplex_handle sh1, Simplex_handle sh2)->bool {
              return st.filtration(sh1) < st.filtration(sh2);   });

  std::cout << "filtration_value:  v0  v1 ... vd \n";
  for(auto sh : sh_range) {
    std::cout << st.filtration(sh) << ":  ";
    for(auto v : st.simplex_vertex_range(sh)) { std::cout << v << " "; }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  std::cout << "Morse matching on the full complex by Benedetti-Lutz heuristic:\n";

  std::cout << "   Range of simplices to match:\n";
  for(auto it = sh_range.rbegin(); it != sh_range.rend(); ++it) {
    std::cout << "   ";
    for(auto v : st.simplex_vertex_range(*it)) { std::cout << v << " "; }
    std::cout << std::endl;  
  }

  std::cout << std::endl;
  std::cout << "   Matching:\n";

  Discrete_morse_theory<Simplex_tree> dmt;
  dmt.compute_matching(sh_range, &st);

  for(auto sh : st.complex_simplex_range()) {
    std::cout << "{ ";
    for(auto v : st.simplex_vertex_range(sh)) { std::cout << v << " "; }
    std::cout << "}  ";
    if(st.critical(sh)) { std::cout << "critical"; }
    else { 
      auto p_sh = st.paired_with(sh); 
      std::cout << "paired with  { ";
      for(auto v : st.simplex_vertex_range(p_sh)) { std::cout << v << " "; }
      std::cout << "}";      
    }
    std::cout << std::endl; 
  }

  std::cout << std::endl;
  std::cout << "Clear the Morse matching:\n";

  dmt.clear_matching(&st);

  for(auto sh : st.complex_simplex_range()) {
    std::cout << "{ ";
    for(auto v : st.simplex_vertex_range(sh)) { std::cout << v << " "; }
    std::cout << "}  ";
    if(st.critical(sh)) { std::cout << "critical\n"; }
    else { 
      auto p_sh = st.paired_with(sh); 
      std::cout << "paired with  { ";
      for(auto v : st.simplex_vertex_range(p_sh)) { std::cout << v << " "; }
      std::cout << "} \n";      
    } 
  }

  std::cout << std::endl;
  std::cout << "Morse matching excluding simplices {0}, {0,5}, {4,5}:\n";

  auto it_beg = sh_range.rbegin(); it_beg++; it_beg++; //exclude edges 0-5 and 4-5
  auto it_end = sh_range.rend(); it_end--; //exclude vertex 0

  std::cout << "   Range of simplices to match:\n";
  for(auto it = it_beg; it != it_end; ++it) {
    std::cout << "   ";
    for(auto v : st.simplex_vertex_range(*it)) { std::cout << v << " "; }
    std::cout << std::endl;  
  }

  std::cout << std::endl;
  std::cout << "   Matching:\n";

  dmt.compute_matching(it_beg, it_end, &st);

  for(auto sh : st.complex_simplex_range()) {
    std::cout << "{ ";
    for(auto v : st.simplex_vertex_range(sh)) { std::cout << v << " "; }
    std::cout << "}  ";
    if(st.critical(sh)) { std::cout << "critical\n"; }
    else { 
      auto p_sh = st.paired_with(sh); 
      std::cout << "paired with  { ";
      for(auto v : st.simplex_vertex_range(p_sh)) { std::cout << v << " "; }
      std::cout << "} \n";      
    } 
  }

  return 0;
}
