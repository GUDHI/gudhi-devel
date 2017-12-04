#ifndef COXETER_COMPLEX_
#define COXETER_COMPLEX_

#include <vector>
#include <map>
#include <ctime>
#include <Eigen/Sparse>
#include <algorithm>
#include <CGAL/Epick_d.h>


class Test {
  
public:
  
  template <class Point_range,
            class Coxeter_system>
  Test(Point_range& point_vector, Coxeter_system& cs) {

    using Simplex_id = std::vector<int>;
    using Vertex_id = Simplex_id;    
    using Pointer_range = std::vector<typename Point_range::iterator>;
    using SPMap = std::map<Simplex_id, Pointer_range>;
    using VSMap = std::map<Vertex_id, std::vector<int>>;

    struct Lexicographic_ptr {
      bool operator() (const typename SPMap::iterator &lhs, const typename SPMap::iterator &rhs) {
        return std::lexicographical_compare(lhs->first.begin(), lhs->first.end(), rhs->first.begin(), rhs->first.end());
      }
    };
    
    using SiMap = std::map<typename SPMap::iterator, int, Lexicographic_ptr>;

    clock_t start, end, global_start;
    global_start = clock();
    start = clock();
    SPMap sp_map;
    for (auto p_it = point_vector.begin(); p_it != point_vector.end(); ++p_it) {
      Simplex_id s_id = cs.alcove_coordinates(*p_it, 1); 
      auto find_it = sp_map.find(s_id);
      if (find_it == sp_map.end())
        sp_map.emplace(s_id, Pointer_range(1, p_it));
      else
        find_it->second.push_back(p_it);
    }
    end = clock();
    double time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    std::cout << "Computed alcove coordinate map in " << time << " s. \n";
      
    start = clock();
    SiMap si_map;
    int si_index = 0;
    for (auto m_it = sp_map.begin(); m_it != sp_map.end(); ++m_it, si_index++)
      si_map.emplace(m_it, si_index);
    end = clock();
    time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    std::cout << "Computed alcove index  map in " << time << " s. \n";
      
    start = clock();
    VSMap vs_map;
    for (auto si_it = si_map.begin(); si_it != si_map.end(); ++si_it) {
      std::vector<Vertex_id> vertices = cs.vertices_of_alcove(si_it->first->first);
      for (Vertex_id v: vertices) {
        auto find_it = vs_map.find(v);
        if (find_it == vs_map.end())
          vs_map.emplace(v, std::vector<int>(1, si_it->second));
        else
          find_it->second.push_back(si_it->second);    
      }
    }
    end = clock();
    time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    std::cout << "Computed vertex alcove map in " << time << " s. \n";
    end = clock();
    time = static_cast<double>(end - global_start) / CLOCKS_PER_SEC;
    std::cout << "Total time: " << time << " s. \n";
      
    std::size_t max_dim = 0; 
    for (auto m: vs_map) {
      if (m.second.size()-1 > max_dim)
        max_dim = m.second.size()-1;
    }
    std::cout << "Dimension of the complex is " << max_dim << ".\n";    
  }

};
#endif
