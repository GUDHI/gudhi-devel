#ifndef COXETER_COMPLEX_
#define COXETER_COMPLEX_

#include <vector>
#include <map>
#include <ctime>
#include <Eigen/Sparse>
#include <algorithm>
#include <utility>
#include <CGAL/Epick_d.h>

namespace Gudhi {

template <class Point_range,
          class Coxeter_system>
class Coxeter_complex {
  
public:
  
  using Simplex_id = std::vector<int>;
  using Vertex_id = Simplex_id;    

  using Index_range = std::vector<std::size_t>;
  using Vertex_map = std::map<Vertex_id, Index_range>;

  using Point_pointers = std::vector<typename Point_range::const_iterator>;
  using Vertex_pointers = std::vector<typename Vertex_map::iterator>;
  using Alcove = std::tuple<std::size_t, Point_pointers, Vertex_pointers>;
  using Alcove_map = std::map<Simplex_id, Alcove>;
  
  const Point_range& point_vector_;
  const Coxeter_system& cs_;
  
  Alcove_map a_map;
  Vertex_map v_map;
  std::size_t max_id;
  
  template <class SPMap_iterator>
  void subdivide_cell(SPMap_iterator sp_it) {
  }

  Coxeter_complex(const Point_range& point_vector, const Coxeter_system& cs)
    : point_vector_(point_vector), cs_(cs), max_id(0) {
    clock_t start, end, global_start;
    double time;
    global_start = clock();
    start = clock();
    for (auto p_it = point_vector.begin(); p_it != point_vector.end(); ++p_it) {
      Simplex_id s_id = cs.alcove_coordinates(*p_it, 1); 
      auto find_it = a_map.find(s_id);
      if (find_it == a_map.end())
        a_map.emplace(s_id, std::make_tuple(max_id++, Point_pointers(1, p_it), Vertex_pointers()));
      else
        std::get<1>(find_it->second).push_back(p_it);
    }
    end = clock();
    time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    std::cout << "Computed alcove coordinate map in " << time << " s. \n";
            
    start = clock();
    for (auto a_it = a_map.begin(); a_it != a_map.end(); ++a_it) {
      std::vector<Vertex_id> vertices = cs.vertices_of_alcove(a_it->first);
      for (Vertex_id v: vertices) {
        auto find_it = v_map.find(v);
        if (find_it == v_map.end())
          v_map.emplace(v, Index_range(1, std::get<0>(a_it->second)));
        else
          find_it->second.push_back(std::get<0>(a_it->second));   
      }
    }
    end = clock();
    time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    std::cout << "Computed vertex alcove map in " << time << " s. \n";
    end = clock();
    time = static_cast<double>(end - global_start) / CLOCKS_PER_SEC;
    std::cout << "Total time: " << time << " s. \n";
      
    std::size_t max_dim = 0; 
    for (auto m: v_map) {
      if (m.second.size()-1 > max_dim)
        max_dim = m.second.size()-1;
    }
    std::cout << "Dimension of the complex is " << max_dim << ".\n";    
  }

};

}
#endif
