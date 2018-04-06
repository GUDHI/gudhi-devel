#ifndef COXETER_COMPLEX_
#define COXETER_COMPLEX_

#define DEBUG_TRACES

#include <string>
#include <map>
#include <utility>

#include <gudhi/Hasse_diagram_persistence.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/dynamic_bitset.hpp>

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
    typedef std::map<Id, Fields > Id_v_map;
    typedef std::map<Graph_v, typename Id_v_map::iterator> V_id_map;
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

  // Compute the ordered partition of an alcove with respect to a vertex
  using Part = boost::dynamic_bitset<>;
  using Ordered_partition = std::vector<Part>;
  Ordered_partition partition(const Id& a_id, const Id& v_id) {
    int d = v_id.size();
    std::vector<int> point;
    for (int l = 0; l < d+1; ++l)
      point.push_back(l);

    auto a_it = a_id.begin();
    for (int j = 1; j < d+1; ++j) {
      int v_coord = 0;
      for (int i = j-1; i >= 0; --i) {
        v_coord += v_id[i];
        point[i] -= (*a_it - v_coord);
        point[j] += (*a_it - v_coord);
        a_it++;
      }
    }
    Ordered_partition op(d+1, Part(d+1));
    for (int l = 0; l < d+1; ++l)
      op[l][point[l]] = 1;
    return op;
  }

  template <class PCMap>
  bool rec_subfaces_are_present(Ordered_partition& op_face,
                                const Ordered_partition& op,
                                int i, int j, int d,
                                const PCMap& pc_map) {
    while (j < d+1 && op[i][j] != 1)
      j++;
    if (j == d+1)
      return (op_face[i].none() || op_face[i+1].none() ||
              pc_map.find(op_face) != pc_map.end());
    bool present = true;
    op_face[i][j] = 1;
    op_face[i+1][j] = 0;
    present = rec_subfaces_are_present(op_face, op, i, j+1, d, pc_map);
    if (present) {
      op_face[i][j] = 0;
      op_face[i+1][j] = 1;
      present = rec_subfaces_are_present(op_face, op, i, j+1, d, pc_map);
    }
    return present;
  }
  
  template <class PCMap>
  bool subparts_are_present(const Ordered_partition& op, const PCMap& pc_map) {
    int d = av_graph_.v_map.begin()->first.size();
    bool present = true;
    for (unsigned i = 0; i < op.size() && present; ++i) {
      Ordered_partition op_face(op.size()+1, Part(d+1));
      for (unsigned j = 0; j < i; ++j)
        op_face[j] = op[j];
      for (unsigned j = i+2; j < op.size()+1; ++j)
        op_face[j] = op[j-1];
      present = rec_subfaces_are_present(op_face, op, i, 0, d, pc_map);
    }
    return present;
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
      std::cout << " vertices: [ ";
      typename Graph::out_edge_iterator e_it, e_end;
      for (std::tie(e_it, e_end) = boost::out_edges(m.second.gv, graph); e_it != e_end; ++e_it) {
        std::cout << inv_map[boost::target(*e_it, graph)]->first << " ";
      }
      std::cout << "]";
      std::cout << "\n";
    }
    std::cout << "Vertex map size is " << v_map.size() << ".\n";    
    std::cout << "VMap:\n";
    for (auto m: v_map) {
      std::cout << m.first << ": filt=" << m.second.f;    
      std::cout << " alcoves: [ ";
      typename Graph::in_edge_iterator e_it, e_end;
      for (std::tie(e_it, e_end) = boost::in_edges(m.second.gv, graph); e_it != e_end; ++e_it) {
        std::cout << inv_map[boost::source(*e_it, graph)]->first << " ";
      }
      std::cout << "]";
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

  // Computing Voronoi skeleton as in the paper
  // k is the desired dimension
  void voronoi_skeleton(int k) {
    using Hasse_cell = Gudhi::Hasse_diagram::Hasse_diagram_cell<int, double, double>;
    using PCMap = std::map<Ordered_partition, Hasse_cell*>;
    using CPMap = std::map<Hasse_cell*, typename PCMap::iterator>;
    using IVMap_it = typename Alcove_vertex_graph::Id_v_map::iterator;
    struct ACMap_comparator {
      bool operator() (const IVMap_it& l_it, const IVMap_it& r_it) {
        return l_it->first < r_it->first;
      }
    };
    using ACMap = std::map<IVMap_it, Hasse_cell*, ACMap_comparator>;

    std::vector<std::vector<Hasse_cell*>> hasse_diagram(k+1);
    ACMap ac_map;
    auto& v_map = av_graph_.v_map;
    auto& inv_map = av_graph_.inv_map;
    auto& graph = av_graph_.graph;
    for (auto v_pair: v_map) {
      PCMap* pc_map_faces = new PCMap(), *pc_map_cofaces;
      CPMap* cp_map_faces = new CPMap(), *cp_map_cofaces;
      typename Graph::in_edge_iterator e_it, e_end;
      // Dual vertex insertion
      for (std::tie(e_it, e_end) = boost::in_edges(v_pair.second.gv, graph); e_it != e_end; ++e_it) {
        auto a_pair_it = inv_map[boost::source(*e_it, graph)];
        Hasse_cell* cell = new Hasse_cell(0, a_pair_it->second.f);
        auto ac_emplace_result = ac_map.emplace(a_pair_it, cell);
        typename PCMap::iterator pc_map_it;
        if (ac_emplace_result.second) {
          pc_map_it = pc_map_faces->emplace(partition(a_pair_it->first, v_pair.first), cell).first;
          hasse_diagram[0].push_back(cell);
        }
        else {
          pc_map_it = pc_map_faces->emplace(partition(a_pair_it->first, v_pair.first),
                                            ac_emplace_result.first->second).first;
          delete cell;
        }
        cp_map_faces->emplace(pc_map_it->second, pc_map_it);
      }
      // Dual face insertion
      for (int curr_dim = 1; curr_dim <= k; ++curr_dim) {
        pc_map_cofaces = new PCMap();
        cp_map_cofaces = new CPMap();
        PCMap candidates; // Some are repetitions with other Voronoi cells
        for (auto pc_pair: *pc_map_faces) {
          int d = v_pair.first.size();
          for (int i = 0; i < d+1-curr_dim; ++i) {
            Ordered_partition op(d+1-curr_dim, Part(d+1));
            for (int j = 0; j <= i; ++j)
              op[j] |= pc_pair.first[j];
            for (int j = i; j < d+1-curr_dim; ++j)
              op[j] |= pc_pair.first[j+1];
            subparts_are_present(op, *pc_map_faces);
          }
        }
        delete pc_map_faces;
        delete cp_map_faces;
        pc_map_faces = pc_map_cofaces;
        cp_map_faces = cp_map_cofaces;
      }
      delete pc_map_cofaces;
      delete cp_map_cofaces;      
    }
  }
  
};

}
#endif
