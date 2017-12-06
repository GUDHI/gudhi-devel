#ifndef COXETER_COMPLEX_
#define COXETER_COMPLEX_

#include <vector>
#include <map>
#include <ctime>
#include <Eigen/Sparse>
#include <algorithm>
#include <utility>
#include <CGAL/Epick_d.h>

#include "../../example/cxx-prettyprint/prettyprint.hpp"


namespace Gudhi {

template <class Point_range,
          class Coxeter_system>
class Coxeter_complex {
  
public:
  
  using Alcove_id = std::vector<int>;
  using Vertex_id = Alcove_id;    

  using Index_range = std::vector<std::size_t>;
  using Vertex_map = std::map<Vertex_id, Index_range>;

  using Point_pointers = std::vector<typename Point_range::const_iterator>;
  using Vertex_pointers = std::vector<typename Vertex_map::iterator>;
  using Alcove = std::tuple<std::size_t, Point_pointers, Vertex_pointers>;
  using Alcove_map = std::map<Alcove_id, Alcove>;
  
  const Point_range& point_vector_;
  const Coxeter_system& cs_;
  
  Alcove_map a_map;
  Vertex_map v_map;
  std::size_t max_id;
  
  template <class AMap_iterator>
  void subdivide_alcove(AMap_iterator a_it) {
    for (auto v_it: std::get<2>(a_it->second)) {
      auto find_it = std::find(v_it->second.begin(), v_it->second.end(), std::get<0>(a_it->second));
      v_it->second.erase(find_it);
      if (v_it->second.empty())
        v_map.erase(v_it);
    }
    for (auto p_it: std::get<1>(a_it->second)) {
      Alcove_id s_id = cs_.alcove_coordinates(*p_it, 2 * a_it->first[0]); 
      auto new_a_it = a_map.find(s_id);
      if (new_a_it == a_map.end()) {
        auto success_pair = a_map.emplace(s_id, std::make_tuple(max_id++, Point_pointers(1, p_it), Vertex_pointers()));
        new_a_it = success_pair.first;
        std::vector<Vertex_id> vertices = cs_.vertices_of_alcove(s_id);
        for (Vertex_id v: vertices) {
          auto new_v_it = v_map.find(v);
          if (new_v_it == v_map.end()) {
            auto success_pair = v_map.emplace(v, Index_range(1, std::get<0>(new_a_it->second)));
            new_v_it = success_pair.first;
            // add the coarser alcoves to the adjacency of the new vertex on the new level
            if (v[0] == new_a_it->first[0]) {
              auto m_it = a_map.begin();
              while (m_it != a_map.end() && m_it->first[0] < new_a_it->first[0]) {
                if (cs_.is_adjacent(new_v_it->first, m_it->first)) {
                  std::get<2>(m_it->second).push_back(new_v_it);
                  new_v_it->second.push_back(std::get<0>(m_it->second));
                }
                m_it++;
              }
            }
          }
          else
            new_v_it->second.push_back(std::get<0>(new_a_it->second));
          std::get<2>(new_a_it->second).push_back(new_v_it);
        }
      }
      else
        std::get<1>(new_a_it->second).push_back(p_it);      
    }
    a_map.erase(a_it);
  }

  Coxeter_complex(const Point_range& point_vector, const Coxeter_system& cs)
    : point_vector_(point_vector), cs_(cs), max_id(0) {
    clock_t start, end, global_start;
    double time;
    global_start = clock();
    start = clock();
    for (auto p_it = point_vector.begin(); p_it != point_vector.end(); ++p_it) {
      Alcove_id s_id = cs.alcove_coordinates(*p_it, 1); 
      auto a_it = a_map.find(s_id);
      if (a_it == a_map.end())
        a_map.emplace(s_id, std::make_tuple(max_id++, Point_pointers(1, p_it), Vertex_pointers()));
      else
        std::get<1>(a_it->second).push_back(p_it);
    }
    end = clock();
    time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    std::cout << "Computed alcove coordinate map in " << time << " s. \n";
            
    start = clock();
    for (auto a_it = a_map.begin(); a_it != a_map.end(); ++a_it) {
      std::vector<Vertex_id> vertices = cs_.vertices_of_alcove(a_it->first);
      for (Vertex_id v: vertices) {
        auto v_it = v_map.find(v);
        if (v_it == v_map.end()) {
          auto success_pair = v_map.emplace(v, Index_range(1, std::get<0>(a_it->second)));
          v_it = success_pair.first;
        }
        else
          v_it->second.push_back(std::get<0>(a_it->second));
        std::get<2>(a_it->second).push_back(v_it);
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
    std::cout << "Dimension of the complex is " << max_dim << ".\n\n";    

    // subdivision part
    
    // std::cout << "AMap:\n";
    // for (auto m: a_map) 
    //   std::cout << m.first << ": " << std::get<0>(m.second) << ", "
    //             << "size=" << std::get<1>(m.second).size() << std::endl;    
    // std::cout << "\n";
    
    // std::cout << "VMap:\n";
    // for (auto m: v_map) 
    //   std::cout << m.first << ": " << m.second << std::endl;
    // std::cout << "\n";
    
    // bool subdivision_needed = true;
    // int current_level = 1;
    // while (subdivision_needed) {
    //   std::cout << "Subdivision level " << 2*current_level << std::endl;
    //   subdivision_needed = false;
    //   auto a_it = a_map.begin();
    //   while (a_it != a_map.end() && a_it->first[0] <= current_level) {
    //     if (std::get<1>(a_it->second).size() > 2) {
    //       subdivision_needed = true;
    //       subdivide_alcove(a_it++);
    //     }
    //     else
    //       a_it++;
    //   }
    //   current_level *= 2;
    //   std::cout << "AMap:\n";
    //   for (auto m: a_map) 
    //     std::cout << m.first << ": " << std::get<0>(m.second) << ", "
    //               << "size=" << std::get<1>(m.second).size() << std::endl;    
    //   std::cout << "\n";
    
    //   std::cout << "VMap:\n";
    //   for (auto m: v_map) 
    //     std::cout << m.first << ": " << m.second << std::endl;
    //   std::cout << "\n";

    // }
  }
 
};

}
#endif
