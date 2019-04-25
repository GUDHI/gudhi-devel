#include <iostream>
#include <gudhi/Coxeter_triangulation/Permutation_iterator.h>
#include <gudhi/Coxeter_triangulation/Combination_iterator.h>
#include <gudhi/Coxeter_triangulation/Integer_combination_iterator.h>
#include <gudhi/Coxeter_triangulation/Set_partition_iterator.h>
#include <gudhi/Coxeter_triangulation/Ordered_set_partition_iterator.h>

#include <gudhi/Coxeter_triangulation/Freudenthal_representation.h>
#include <gudhi/Coxeter_triangulation/Face_iterator.h>
#include <gudhi/Coxeter_triangulation/Coface_iterator.h>
#include "../example/cxx-prettyprint/prettyprint.hpp"

int main() {
  Permutation_iterator p_it(5), p_end;
  unsigned counter = 0;
  // for (; p_it != p_end; ++p_it, ++counter)
  //   std::cout << *p_it << "\n";
  // std::cout << "counter = " << counter << "\n";

  // std::cout << "New cycle...\n";
  // counter = 0;
  // p_it.reinitialize();
  // for (; p_it != p_end; ++p_it, ++counter)
  //   std::cout << *p_it << "\n";
  // std::cout << "counter = " << counter << "\n";

  // Combination_iterator c_it(10,3), c_end; 
  // for (; c_it != c_end; ++c_it, ++counter)
  //   std::cout << *c_it << "\n";
  // std::cout << "counter = " << counter << "\n";

  // std::vector<unsigned> bounds = {2};
  // Integer_combination_iterator i_it(1, 1, bounds), i_end; 
  // for (; i_it != i_end; ++i_it, ++counter)
  //   std::cout << *i_it << "\n";
  // std::cout << "counter = " << counter << "\n";

  // Set_partition_iterator s_it(2,1), s_end; 
  // for (; s_it != s_end; ++s_it, ++counter) {
  //   for (uint i = 0; i < s_it->size(); i++)
  //     std::cout << (*s_it)[i] << " ";
  //   std::cout << "\n";
  // }
  // std::cout << "counter = " << counter << "\n";

  // Ordered_set_partition_iterator o_it(6,3), o_end; 
  // for (; o_it != o_end; ++o_it, ++counter) {
  //   for (uint i = 0; i < o_it->size(); i++)
  //     std::cout << (*o_it)[i] << " ";
  //   std::cout << "\n";
  // }
  // std::cout << "counter = " << counter << "\n";

  typedef std::vector<int> Vertex;
  typedef std::vector<std::vector<uint> > Ordered_partition;
  typedef Freudenthal_representation<Vertex, Ordered_partition> FR;
  typedef Face_iterator<FR> Face_iterator;
  typedef Coface_iterator<FR> Coface_iterator;
  FR fr = {Vertex(8,0), Ordered_partition({std::vector<uint>({1}),
  					   std::vector<uint>({3,6}),
  					   std::vector<uint>({4}),
  					   std::vector<uint>({2,5}),
  					   std::vector<uint>({0,7,8})})};
  // FR fr = {Vertex(3,0), Ordered_partition({std::vector<uint>({0,1}),
  // 					   std::vector<uint>({2,3})})};

  std::cout << fr.vertex << " ";
  for (uint i = 0; i < fr.partition.size(); i++)
    std::cout << fr.partition[i] << " ";
  std::cout << "\n\n";
  
  // Face_iterator f_it(fr, 2), f_end;
  // for (; f_it != f_end; ++f_it, ++counter) {
  //   std::cout << f_it->vertex << " ";
  //   for (uint i = 0; i < f_it->partition.size(); i++)
  //     std::cout << f_it->partition[i] << " ";
  //   std::cout << "\n";
  // }
  // std::cout << "counter = " << counter << "\n";

  Coface_iterator f_it(fr, 6), f_end;
  for (; f_it != f_end; ++f_it, ++counter) {
    std::cout << f_it->vertex << " ";
    for (uint i = 0; i < f_it->partition.size(); i++)
      std::cout << f_it->partition[i] << " ";
    std::cout << "\n";
  }
  std::cout << "counter = " << counter << "\n";

}
