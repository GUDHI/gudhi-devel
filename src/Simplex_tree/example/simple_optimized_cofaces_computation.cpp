#include <iostream>
#include <gudhi/Simplex_tree.h>
#include <chrono>

// Types definition, allowing the fast cofaces computation
using Simplex_tree      = 
                Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_zigzag_persistence>;

int main(int argc, char* argv[])
{ 
  // traverse the entire oscillating Rips zigzag filtration 
  Simplex_tree st;
  std::vector<int> simplex = {0,1,2,3,4};
  std::vector<int> face = {1,3};

  st.insert_simplex_and_subfaces(simplex);

  size_t num_simp=0;
  std::cout << "Complex:\n";
  for(auto sh : st.complex_simplex_range()) {
    ++num_simp;
    for(auto v : st.simplex_vertex_range(sh)) {
      std::cout << v << " ";
    }
    std::cout << std::endl;
  }

  auto sh_face = st.find(face);
  size_t num_cof=0;

  std::cout << "Cofaces with optimized search:\n";
  for(auto cof_sh : st.star_simplex_range(sh_face)) {
    ++num_cof;
    for(auto v : st.simplex_vertex_range(cof_sh)) {
      std::cout << v << " ";
    }
    std::cout << std::endl;
  }
  std::cout << num_simp << " simplices, and " << num_cof << " cofaces.\n\n"; 

  return 0;
}
 
