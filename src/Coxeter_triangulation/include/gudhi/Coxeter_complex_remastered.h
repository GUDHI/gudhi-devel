#ifndef COXETER_COMPLEX_
#define COXETER_COMPLEX_

#define DEBUG_TRACES

#include <string>
#include <map>
#include <utility>

#include <gudhi/Hasse_diagram_persistence.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Simplex_tree.h>

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
#ifndef CC_A_V_VISITORS
      if (!vertices.empty())
        std::cout << "Vertices are not empty!\n";
#endif
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
#ifndef CC_A_V_VISITORS
    auto& inv_map = av_graph_.inv_map;
    auto& v_map = av_graph_.v_map;
    auto& graph = av_graph_.graph;
    for (auto a_pair: av_graph_.a_map) {
      std::vector<Id> vertices = cs_.vertices_of_alcove(a_pair.first);
      Graph_v a_v = a_pair.second.gv;
      for (Id v_id: vertices) {
        auto v_it = v_map.find(v_id);
        if (v_it == v_map.end()) {
          Graph_v v_v = boost::add_vertex(graph);
          v_it = v_map.emplace(v_id, Fields(v_v, a_pair.second.f)).first;
          inv_map.emplace(v_v, v_it);
          boost::add_edge(a_v, v_v, graph);
        }
        else {
          v_it->second.f = std::min(v_it->second.f, a_pair.second.f);
          boost::add_edge(a_v, v_it->second.gv, graph);
        }
      }
    }
#endif
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
      std::cout << ", vertices: [ ";
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
      std::cout << ", alcoves: [ ";
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
    auto& v_map = av_graph_.v_map;
    cs_.write_mesh(v_map, range, file_name);
  }

private:

  // The map from vertices to indices
  using Id_it = typename Alcove_vertex_graph::Id_v_map::const_iterator;
  struct Pointer_compare {
    typedef Id_it Pointer;
    bool operator()(const Pointer& lhs, const Pointer& rhs) const { 
      return lhs->first < rhs->first;
    }
  };
  using Vertex_index_map = std::map<Id_it, std::size_t, Pointer_compare>;

  /* An iterator over Simplex tree simplices */
  template <class Simplex_tree>
  class Simplex_tree_simplex_iterator : public boost::iterator_facade< Simplex_tree_simplex_iterator<Simplex_tree>,
                                                                       std::vector<std::size_t> const,
                                                                       boost::forward_traversal_tag> {
  private:
    typename Simplex_tree::Complex_simplex_iterator it_;
    Simplex_tree& stree_;
    std::vector<std::size_t> value_; 
    friend class boost::iterator_core_access;
    void update_value() {
      if (it_ != stree_.complex_simplex_range().end()) {
        auto range = stree_.simplex_vertex_range(*it_);
        value_ = std::vector<size_t>(range.begin(), range.end());
      }
    }
    bool equal(Simplex_tree_simplex_iterator<Simplex_tree> const& other) const {
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
    Simplex_tree_simplex_iterator(typename Simplex_tree::Complex_simplex_iterator it,
                                  Simplex_tree& stree)
      : it_(it), stree_(stree) {
      update_value();
    }
  };

  /* Alcove iterator */
  class No_filtration_alcove_iterator : public boost::iterator_facade< No_filtration_alcove_iterator,
                                                                       std::vector<std::size_t> const,
                                                                       boost::forward_traversal_tag> {
  private:
    using Alcove_map = typename Alcove_vertex_graph::Id_v_map;
    typename Alcove_map::const_iterator it_;
    Alcove_vertex_graph const& av_graph_;
    Vertex_index_map const& vi_map_;
    std::vector<std::size_t> value_;
    friend class boost::iterator_core_access;
    void update_value() {
      auto& a_map = av_graph_.a_map;
      auto& inv_map = av_graph_.inv_map;
      auto& graph = av_graph_.graph;    
      value_.clear();
      if (it_ != a_map.end()) {
        typename Graph::out_edge_iterator e_it, e_end;
        for (std::tie(e_it, e_end) = boost::out_edges(it_->second.gv, graph); e_it != e_end; ++e_it)
          value_.push_back(vi_map_.at(inv_map.at(boost::target(*e_it, graph))));
        std::sort(value_.begin(), value_.end());
      }
    }
    bool equal(No_filtration_alcove_iterator const& other) const {
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
    No_filtration_alcove_iterator(typename Alcove_map::const_iterator it,
                                  Alcove_vertex_graph const& av_graph,
                                  Vertex_index_map const& vi_map)
      : it_(it), av_graph_(av_graph), vi_map_(vi_map) {
      update_value();
    }
    int dimension() {
      return value_.size()-1;
    }
  };
  
public:
  void write_mesh(std::string file_name = "toplex.mesh") const {
    Vertex_index_map vi_map;
    std::size_t index = 0;
    auto& v_map = av_graph_.v_map;
    auto& a_map = av_graph_.a_map;
    for (Id_it v_it = v_map.begin(); v_it != v_map.end(); ++v_it, ++index)
      vi_map.emplace(v_it, index);
    
    typedef boost::iterator_range<No_filtration_alcove_iterator> Max_simplex_range;
    Max_simplex_range range(No_filtration_alcove_iterator(a_map.begin(), av_graph_, vi_map),
                            No_filtration_alcove_iterator(a_map.end(), av_graph_, vi_map));
    using Simplex_tree = Gudhi::Simplex_tree<>;
    Simplex_tree stree;
    for (auto m: range)
      stree.insert_simplex_and_subfaces(m);
    typedef Simplex_tree_simplex_iterator<Simplex_tree> Simplex_tree_iterator;
    typedef boost::iterator_range<Simplex_tree_iterator> Simplex_tree_range;
    Simplex_tree_range simplex_tree_range(Simplex_tree_iterator(stree.complex_simplex_range().begin(), stree),
                                          Simplex_tree_iterator(stree.complex_simplex_range().end(), stree));
    cs_.write_mesh(v_map, simplex_tree_range, file_name);
  }

    template <class Filtered_simplex_range>
  void write_bb(Filtered_simplex_range& range, std::string file_name) const {
  }
  
  void write_bb(std::string file_name = "toplex.bb") const {
  }

  void collapse(bool pers_out = true) {
  }

private:

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

  // Check if all the subfaces are present for the face given by a partition
  using Hasse_cell = Gudhi::Hasse_diagram::Hasse_diagram_cell<int, double, double>;
  using Hasse_boundary = std::vector<std::pair<Hasse_cell*, int> >;
  template <class PCMap>
  bool rec_subfaces_are_present(Ordered_partition& op_face,
                                const Ordered_partition& op,
                                int i, int j, int d, double& f,
                                Hasse_boundary& boundary,
                                const PCMap& pc_map) {
    while (j < d+1 && op[i][j] != 1)
      j++;
    if (j == d+1) {
      if (op_face[i].none() || op_face[i+1].none())
        return true;
      auto pc_pair_it = pc_map.find(op_face);
      if (pc_pair_it != pc_map.end()) {
        f = std::max(f, pc_pair_it->second->get_filtration());
        boundary.emplace_back(std::make_pair(pc_pair_it->second, 1));
        return true;
      }
      else
        return false;
    }
    bool present = true;
    op_face[i][j] = 1;
    op_face[i+1][j] = 0;
    present = rec_subfaces_are_present(op_face, op, i, j+1, d, f, boundary, pc_map);
    if (present) {
      op_face[i][j] = 0;
      op_face[i+1][j] = 1;
      present = rec_subfaces_are_present(op_face, op, i, j+1, d, f, boundary, pc_map);
    }
    return present;
  }
  
  template <class PCMap>
  bool subparts_are_present(const Ordered_partition& op,
                            const PCMap& pc_map,
                            double& f,
                            Hasse_boundary& boundary) {
    int d = av_graph_.v_map.begin()->first.size();
    bool present = true;
    for (unsigned i = 0; i < op.size() && present; ++i) {
      Ordered_partition op_face(op.size()+1, Part(d+1));
      for (unsigned j = 0; j < i; ++j)
        op_face[j] = op[j];
      for (unsigned j = i+2; j < op.size()+1; ++j)
        op_face[j] = op[j-1];
      present = rec_subfaces_are_present(op_face, op, i, 0, d, f, boundary, pc_map);
    }
    return present;
  }

public:  
  // Computing Voronoi skeleton as in the paper
  // k is the desired dimension
  void voronoi_skeleton(int k) {
    using PCMap = std::map<Ordered_partition, Hasse_cell*>; 
    using IVMap_it = typename Alcove_vertex_graph::Id_v_map::iterator;
    struct ACMap_comparator {
      bool operator() (const IVMap_it& l_it, const IVMap_it& r_it) {
        return l_it->first < r_it->first;
      }
    };
    using ACMap = std::map<IVMap_it, Hasse_cell*, ACMap_comparator>;

    struct Hasse_cell_comparator {
      bool operator() (Hasse_cell* l_it, Hasse_cell* r_it) {
        if (l_it->get_dimension() == 0)
          return l_it < r_it; 
        else
          return l_it->get_boundary() < r_it->get_boundary();
      }
    };
    using Hasse_diagram = std::set<Hasse_cell*, Hasse_cell_comparator>;
    Hasse_diagram hasse_diagram;
    ACMap ac_map;
    auto& v_map = av_graph_.v_map;
    auto& inv_map = av_graph_.inv_map;
    auto& graph = av_graph_.graph;
    for (auto v_pair: v_map) {
#ifdef DEBUG_TRACES
      std::cout << "Current vertex: " << v_pair.first << std::endl;
#endif
      PCMap* pc_map_faces = new PCMap(), *pc_map_cofaces;
      typename Graph::in_edge_iterator e_it, e_end;
      // Dual vertex insertion
      for (std::tie(e_it, e_end) = boost::in_edges(v_pair.second.gv, graph); e_it != e_end; ++e_it) {
        auto a_pair_it = inv_map[boost::source(*e_it, graph)];
        Hasse_cell* cell = new Hasse_cell(0, a_pair_it->second.f);
        auto ac_emplace_result = ac_map.emplace(a_pair_it, cell);
        typename PCMap::iterator pc_map_it;
        if (ac_emplace_result.second) {
          pc_map_it = pc_map_faces->emplace(partition(a_pair_it->first, v_pair.first), cell).first;
          hasse_diagram.emplace(cell);
        }
        else {
          pc_map_it = pc_map_faces->emplace(partition(a_pair_it->first, v_pair.first),
                                            ac_emplace_result.first->second).first;
          delete cell;
        }
#ifdef DEBUG_TRACES
        std::cout << "dual vertex for " << a_pair_it->first
                  << ", partition=" << partition(a_pair_it->first, v_pair.first) << std::endl;
#endif
      }
      // Dual face insertion
      for (int curr_dim = 1; curr_dim <= k; ++curr_dim) {
        pc_map_cofaces = new PCMap();
        for (auto pc_pair: *pc_map_faces) {
          int d = v_pair.first.size();
          for (int i = 0; i < d+1-curr_dim; ++i) {
            Ordered_partition op(d+1-curr_dim, Part(d+1));
            for (int j = 0; j <= i; ++j)
              op[j] |= pc_pair.first[j];
            for (int j = i; j < d+1-curr_dim; ++j)
              op[j] |= pc_pair.first[j+1];
            double f = 0;
            Hasse_boundary boundary;
            if (pc_map_cofaces->find(op) == pc_map_cofaces->end() &&
                subparts_are_present(op, *pc_map_faces, f, boundary)) {
              Hasse_cell* cell = new Hasse_cell(i, f);
              std::sort(boundary.begin(), boundary.end());
              cell->get_boundary() = boundary;
#ifdef DEBUG_TRACES
              if (pc_map_cofaces->emplace(op, cell).second)
                std::cout << "dual face " << op << std::endl;
#else
              pc_map_cofaces->emplace(op, cell);
#endif
            }
          }
        }
        for (auto pc_pair: *pc_map_cofaces)
          hasse_diagram.emplace(pc_pair.second);
        delete pc_map_faces;
        pc_map_faces = pc_map_cofaces;
      }
      delete pc_map_cofaces;
    }
    typedef Gudhi::Hasse_diagram::Hasse_diagram_persistence<Hasse_cell> Hasse_pers_vector;
    Hasse_pers_vector hdp(std::vector<Hasse_cell*>(hasse_diagram.begin(),
                                                   hasse_diagram.end()));  
    typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;
    typedef Gudhi::persistent_cohomology::Persistent_cohomology
                          <Hasse_pers_vector, Field_Zp> Persistent_cohomology;
    
    Persistent_cohomology pcoh(hdp, true);  
    unsigned field_characteristic = 11;
    double min_persistence = 0;
    
    pcoh.init_coefficients(field_characteristic);
    pcoh.increment_dimension();
    pcoh.compute_persistent_cohomology(min_persistence);
    std::cout << pcoh.persistent_betti_numbers(0,0) << std::endl;
    std::ofstream out("persdiag_vor.out");
    pcoh.output_diagram(out);
    out.close();
#ifdef VERBOSE_DEBUG_TRACES
    std::cout << "Hasse diagram:\n" << hdp << std::endl;
#endif
    for (auto c_ptr: hasse_diagram) {
      delete c_ptr;
    }
  }
  
};

}
#endif
