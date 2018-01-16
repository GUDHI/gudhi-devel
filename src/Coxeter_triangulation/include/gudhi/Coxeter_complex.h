#ifndef COXETER_COMPLEX_
#define COXETER_COMPLEX_

#include <vector>
#include <map>
#include <ctime>
#include <Eigen/Sparse>
#include <algorithm>
#include <utility>
#include <CGAL/Epick_d.h>

#include <gudhi/SparseMsMatrix.h>
#include <gudhi/Fake_simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Dim_lists.h>
#include <gudhi/Coxeter_complex/Collapse.h>
// #include <gudhi/Coxeter_complex/Simplex_with_cofaces.h>

#include "../../example/cxx-prettyprint/prettyprint.hpp"

#include <boost/iterator/iterator_facade.hpp>
// #include <boost/graph/graph_traits.hpp>
// #include <boost/graph/adjacency_list.hpp>
// #include <boost/graph/connected_components.hpp>

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

  struct Pointer_compare {
    typedef typename Vertex_map::iterator Pointer;
    bool operator()(const Pointer& lhs, const Pointer& rhs) const { 
      return lhs->first < rhs->first;
    }
  };
  using Vertex_index_map = std::map<typename Vertex_map::iterator, int, Pointer_compare>;
  
  const Point_range& point_vector_;
  const Coxeter_system& cs_;
  
  Alcove_map a_map;
  Vertex_map v_map;
  std::size_t max_id;
  Vertex_index_map vi_map;

  
  class Alcove_iterator : public boost::iterator_facade< Alcove_iterator,
                                                         std::vector<std::size_t> const,
                                                         boost::forward_traversal_tag> {
  private:
    typename Alcove_map::iterator it_;
    Alcove_map& a_map_;
    Vertex_index_map& vi_map_;
    std::vector<std::size_t> value_; 
    friend class boost::iterator_core_access;

    void update_value() {
      value_.clear();
      if (it_ != a_map_.end()) {
        for (auto m_it: std::get<2>(it_->second))
          value_.push_back(vi_map_[m_it]);
        std::sort(value_.begin(), value_.end());
      }
    }
    
    bool equal(Alcove_iterator const& other) const {
      return it_ == other.it_;
    }

    std::vector<std::size_t> const& dereference() const {
      return value_;
    }

    void increment() {
      it_++;
      update_value();
    }
    
  public:
    Alcove_iterator(typename Alcove_map::iterator it,
                    Alcove_map& a_map,
                    Vertex_index_map& vi_map)
      : it_(it), a_map_(a_map), vi_map_(vi_map) {
      update_value();
    }
  };
  
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

  Coxeter_complex(const Point_range& point_vector, const Coxeter_system& cs, int init_level=1)
    : point_vector_(point_vector), cs_(cs), max_id(0) {
    clock_t start, end, global_start;
    double time;
    global_start = clock();
    start = clock();
    for (auto p_it = point_vector.begin(); p_it != point_vector.end(); ++p_it) {
      Alcove_id s_id = cs.alcove_coordinates(*p_it, init_level); 
      auto a_it = a_map.find(s_id);
      if (a_it == a_map.end())
        a_map.emplace(s_id, std::make_tuple(max_id++, Point_pointers(1, p_it), Vertex_pointers()));
      else
        std::get<1>(a_it->second).push_back(p_it);
    }
    end = clock();
    time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    std::cout << "Computed alcove coordinate map in " << time << " s. \n";

    // std::cout << "AMap:\n";
    // for (auto m: a_map) 
    //   std::cout << m.first << ": " << std::get<0>(m.second) << ", "
    //             << "size=" << std::get<1>(m.second).size() << std::endl;    
    // std::cout << "\n";
    
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
    std::cout << "Dimension of the nerve is " << max_dim << ".\n\n";    

    // graph part to test the proximity

    // typedef typename boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> Adj_graph;
    // typedef boost::graph_traits<Adj_graph>::vertex_descriptor Vertex_t;
    // typedef boost::graph_traits<Adj_graph>::edge_descriptor Edge_t;
    // typedef boost::graph_traits<Adj_graph>::adjacency_iterator Adj_it;
    // typedef std::pair<Adj_it, Adj_it> Out_edge_it;

    // typedef boost::container::flat_map<std::size_t, Vertex_t> Graph_map;
    // typedef boost::container::flat_map<Vertex_t, std::size_t> Inv_graph_map;    

    // Adj_graph adj_graph;
    // Graph_map g_map;
    // Inv_graph_map ig_map;
    // Vertex_t vert;
    // for (auto a_it = a_map.begin(); a_it != a_map.end(); a_it++) {
    //   vert = boost::add_vertex(adj_graph);
    //   g_map.emplace(std::get<0>(a_it->second), vert);
    //   ig_map.emplace(vert, std::get<0>(a_it->second));
    // }
    // for (auto v_it = v_map.begin(); v_it != v_map.end(); v_it++)
    //   for (auto al_it1 = v_it->second.begin(); al_it1 != v_it->second.end(); ++al_it1)
    //     for (auto al_it2 = al_it1+1; al_it2 != v_it->second.end(); ++al_it2)
    //       boost::add_edge(g_map[*al_it1], g_map[*al_it2], adj_graph);

    /* It was uncommented  \/ */
    
    // std::map<typename Vertex_map::iterator, typename Point_range::value_type> vp_map;
    // for (auto v_it = v_map.begin(); v_it != v_map.end(); v_it++) {
      
    // }
    
    // std::vector<std::vector<bool>> adj_table(a_map.size(), std::vector<bool>(a_map.size(), false));
    // for (auto a_it = a_map.begin(); a_it != a_map.end(); a_it++) {
    //   Vertex_pointers& v_list = std::get<2>(a_it->second);
    //   for (auto v_it = v_list.begin(); v_it != v_list.end(); v_it++)
    //     for (auto aa: (*v_it)->second)
    //       adj_table[std::get<0>(a_it->second)][aa] = true;
    // }

    // for (unsigned i = 0; i != adj_table.size(); i++)
    //   for (unsigned j = 0; j != adj_table[i].size(); j++)
    //     if (!adj_table[i][j]) {
          
    //     }
          
    
    // subdivision part
    
    
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
    
      // std::cout << "VMap:\n";
      // for (auto m: v_map) 
      //   std::cout << m.first << ": " << m.second << std::endl;
      // std::cout << "\n";

    // }
  }

  void write_coxeter_mesh(std::string file_name = "coxeter.mesh") {
    cs_.write_coxeter_mesh(v_map, a_map, file_name);
  }

  template <class Fake_simplex_tree>
  void write_toplex_mesh(Fake_simplex_tree& stree, std::string file_name = "toplex.mesh") const {
    cs_.write_toplex_mesh(v_map, stree, file_name);
  }

  
  void collapse(bool pers_out = true) {
    Fake_simplex_tree stree;

    clock_t start, end;
    double time;
    
    unsigned index = 0;
    for (auto v_it = v_map.begin(); v_it != v_map.end(); ++v_it, ++index)
      vi_map.emplace(v_it, index);
    
    start = clock();
    // for (auto a: a_map) {
    //   std::vector<int> vertices;
    //   for (auto v_it: std::get<2>(a.second))
    //     vertices.push_back(vi_map[v_it]);
    //   stree.insert_simplex_and_subfaces(vertices, 0);
    // }
    // SparseMsMatrix mat(stree);


    typedef boost::iterator_range<Alcove_iterator> Max_simplex_range;
    Max_simplex_range max_simplex_range(Alcove_iterator(a_map.begin(), a_map, vi_map),
                                        Alcove_iterator(a_map.end(), a_map, vi_map));
    end = clock();
    // time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    // std::cout << "Computed the ToplexMap in time " << time << " s. \n";
    // SparseMsMatrix mat(a_map.size(), max_simplex_range);
    
    // auto matrix_formed  = std::chrono::high_resolution_clock::now();
    // std::cout << "Start strong collapse..." << std::endl;

    // Fake_simplex_tree coll_tree = mat.collapsed_tree();
    // auto collapse_done = std::chrono::high_resolution_clock::now();
    // std::cout << "Strong collapse done." << std::endl;
    // auto collapseTime = std::chrono::duration<double, std::milli>(collapse_done- matrix_formed).count();
    // std::cout << "Time for Collapse : " << collapseTime << " ms\n" << std::endl;

    // int originalDim, collDim;
    // int originalNumVert, colNumVert;
    // long originalNumMxSimp, colNumMxSimp;

    // originalDim = stree.dimension();
    // collDim 	= coll_tree.dimension();
		
    // originalNumVert = stree.num_vertices();
    // colNumVert		= coll_tree.num_vertices();

    // originalNumMxSimp 	= stree.num_simplices();
    // colNumMxSimp 		= coll_tree.num_simplices();

    // // std::cout << "Dimension of the collapsed complex is " << collDim << std::endl; 
    // // std::cout << "Number of maximal simplices is " << colNumSimp << std::endl;
    
    // std::cout << "Coxeter complex is of dimension " << originalDim << " with " << originalNumMxSimp << " maximal simplices and " << originalNumVert << " vertices." << std::endl;
    // std::cout << "Collapsed Coxeter complex is of dimension " << collDim << " with " <<  colNumMxSimp << " maximal simplices and " << colNumVert << " vertices." << std::endl;
    
    // std::cout << "** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** " << std::endl;

    Simplex_tree<> coll_stree;
    Collapse coll(max_simplex_range, coll_stree);
    // std::cout << coll_stree << "\n";
    unsigned dim_complex = 0;
    for (auto sh: coll_stree.complex_simplex_range()) {
      if (coll_stree.dimension(sh) > dim_complex)
        dim_complex = coll_stree.dimension(sh);
    }
    std::cout << "Dimension of the collapsed complex is " << dim_complex << "\n";
    
    // if (pers_out) {
    //   std::cout << "Start persistence calculation..." << std::endl;
    //   using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
    //   using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree<>, Field_Zp>;
    //   double min_persistence = 0;
    //   Simplex_tree<> real_stree;
    //   for (auto s: coll_tree.max_simplices()) {
    //     real_stree.insert_simplex_and_subfaces(s);
    //   }
    //   // Sort the simplices in the order of the filtration
    //   real_stree.initialize_filtration();
    //   // Compute the persistence diagram of the complex
    //   Persistent_cohomology pcoh(real_stree);
    //   // initializes the coefficient field for homology
    //   pcoh.init_coefficients(3);
      
    //   pcoh.compute_persistent_cohomology(min_persistence);
    //   std::ofstream out("persdiag_cox.out");
    //   pcoh.output_diagram(out);
    //   out.close();
    //   std::cout << "Persistence complete." << std::endl;

      
    //   witness_complex::Dim_lists<Simplex_tree<>> simplices(real_stree, collDim, 1);
    //   start = clock();
    //   simplices.collapse();
    //   end = clock();
    //   double time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    //   std::cout << "Collapses took " << time << " s. \n";
    //   // Simplicial_complex collapsed_tree;
    //   // for (auto sh: simplices) {
    //   //   std::vector<int> vertices;
    //   //   for (int v: collapsed_tree.simplex_vertex_range(sh))
    //   //     vertices.push_back(v);
    //   //   collapsed_tree.insert_simplex(vertices, simplex_tree.filtration(sh));
    //   // }
    //   std::cout << "The dimension after the simple collapses is " << simplices.dimension() << ".\n";
    // }

    
    // write_toplex_mesh(coll_tree);
    
  }
  
};

}
#endif
