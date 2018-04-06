#ifndef COXETER_COMPLEX_
#define COXETER_COMPLEX_

#define DEBUG_TRACES

#include <string>
#include <map>
#include <utility>
#include <boost/graph/adjacency_list.hpp>

#include "../../example/cxx-prettyprint/prettyprint.hpp"

double num_error = 1.0e-10;

namespace Gudhi {
  
template <class Point_range,
          class Coxeter_system>
class Coxeter_complex {

  // Definition of a filtered alcove
  using Id = typename Coxeter_system::Alcove_id;
  using Filtration = double;
  using Filtered_alcove = typename Coxeter_system::Filtered_alcove;

  // Graph definitions
  using Graph = boost::adjacency_list< boost::listS,          
                                       boost::listS,            
                                       boost::bidirectionalS >;      
  using Graph_v = typename Graph::vertex_descriptor;
  struct Fields {
    Graph_v gv;
    Filtration f;
    Fields(Graph_v gv_in, Filtration f_in) : gv(gv_in), f(f_in) {}
  };
  
  // The data structure that contains the alcoves and the vertices
  struct Alcove_vertex_graph {
    using Id_v_map = std::map<Id, Fields >;
    using V_id_map = std::map<Graph_v, typename Id_v_map::iterator>;
    Id_v_map a_map, v_map;
    V_id_map inv_map;
    Graph graph;
  } av_graph_;
  // The underlying Coxeter system
  const Coxeter_system& cs_;

  // The visitor that inserts the alcoves and vertices to the Alcove_vertex_graph
  struct Alcove_vertex_visitor {
    Alcove_vertex_visitor(typename Point_range::const_iterator& p_it,
                          Alcove_vertex_graph& av_graph,
                          bool& store_points)
      : p_it_(p_it), av_graph_(av_graph), store_points_(store_points) {}
    
    void operator() (const Filtered_alcove& a,
                     const std::vector<Id>& vertices) {
      auto& a_map = av_graph_.a_map;
      auto& v_map = av_graph_.v_map;
      auto& inv_map = av_graph_.inv_map;
      auto& graph = av_graph_.graph;
      auto a_it = a_map.find(a.id);
      if (a_it == a_map.end()) {
        Graph_v a_v = boost::add_vertex(graph);
        a_it = a_map.emplace(a.id, Fields(a_v, a.f)).first;
        inv_map.emplace(a_v, a_it);
        for (Id v_id: vertices) {
          auto v_it = v_map.find(v_id);
          if (v_it == v_map.end()) {
            Graph_v v_v = boost::add_vertex(graph);
            v_it = v_map.emplace(v_id, Fields(v_v, a.f)).first;
            inv_map.emplace(v_v, v_it);
            boost::add_edge(a_v, v_v, graph);
          }
          else {
            v_it->second.f = std::min(v_it->second.f, a.f);
            boost::add_edge(a_v, v_it->second.gv, graph);
          }
        }
      }
      else {
        a_it->second.f = std::min(a_it->second.f, a.f);
        Graph_v a_v = a_it->second.gv;
        for (Id v_id: vertices) {
          auto v_it = v_map.find(v_id);
          if (v_it == v_map.end()) {
            Graph_v v_v = boost::add_vertex(graph);
            v_it = v_map.emplace(v_id, Fields(v_v, a.f)).first;
            inv_map.emplace(v_v, v_it);
            boost::add_edge(a_v, v_v, graph);
          }
          else
            v_it->second.f = std::min(v_it->second.f, a.f);
        }
      }
    }
  private :
    typename Point_range::const_iterator& p_it_;
    Alcove_vertex_graph& av_graph_;
    bool& store_points_;
  };

  // The procedure that fills av_graph_
  void compute_graph(const Point_range& point_vector, double init_level, double eps, bool store_points) {
    for (auto p_it = point_vector.begin(); p_it != point_vector.end(); ++p_it)
      cs_.alcoves_of_ball(*p_it,
                          init_level,
                          eps,
                          Alcove_vertex_visitor(p_it, av_graph_, store_points));
    for (auto m: av_graph_.a_map)
      m.second.f = std::sqrt(m.second.f);
  }
  
public:

  Coxeter_complex(const Point_range& point_vector,
                  const Coxeter_system& cs,
                  double init_level=1,
                  double eps=0,
                  bool store_points = false) : cs_(cs)
  {
    compute_graph(point_vector, init_level, eps, store_points);
#ifdef DEBUG_TRACES
    auto& a_map = av_graph_.a_map;
    auto& v_map = av_graph_.v_map;
    auto& inv_map = av_graph_.inv_map;
    auto& graph = av_graph_.graph;
    std::cout << "Alcove map size is " << a_map.size() << ".\n";    
    std::cout << "AMap:\n";
    for (auto m: a_map) {
      std::cout << m.first << ": filt=" << m.second.f;
      // std::cout << " vertices: [ ";
      // typename Graph::out_edge_iterator e_it, e_end;
      // for (std::tie(e_it, e_end) = boost::out_edges(m.second.gv, graph); e_it != e_end; ++e_it) {
      //   std::cout << inv_map[boost::target(*e_it, graph)]->first << " ";
      // }
      // std::cout << "]";
      std::cout << "\n";
    }
    std::cout << "Vertex map size is " << v_map.size() << ".\n";    
    std::cout << "VMap:\n";
    for (auto m: v_map) {
      std::cout << m.first << ": filt=" << m.second.f;    
      // std::cout << " alcoves: [ ";
      // typename Graph::in_edge_iterator e_it, e_end;
      // for (std::tie(e_it, e_end) = boost::in_edges(m.second.gv, graph); e_it != e_end; ++e_it) {
      //   std::cout << inv_map[boost::source(*e_it, graph)]->first << " ";
      // }
      // std::cout << "]";
      std::cout << "\n";
    }
#endif

    // compute_v_map();
  }

  template <class Simplex_range>
  void write_mesh(Simplex_range& range, std::string file_name) const {
  }
  
  void write_mesh(std::string file_name = "toplex.mesh") const {
  }

    template <class Filtered_simplex_range>
  void write_bb(Filtered_simplex_range& range, std::string file_name) const {
  }
  
  void write_bb(std::string file_name = "toplex.bb") const {
  }

  void collapse(bool pers_out = true) {
  }
  
};

}
#endif
