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

  // The data structure that contains the alcoves and the vertices
  struct Alcove_vertex_graph {
    using Graph = boost::adjacency_list< boost::listS,          
                                         boost::listS,            
                                         boost::bidirectionalS >;      
    std::map<Id, std::pair<typename Graph::vertex_descriptor, double> > a_map, v_map;    
    Graph graph;
    void emplace_alcove(const Filtered_alcove& a) {
      typename Graph::vertex_descriptor v = boost::add_vertex(graph);
      a_map.emplace(a.id, std::make_pair(v, a.f));
    }
  } av_graph_;
  // The underlying Coxeter system
  const Coxeter_system& cs_;

  // The visitor that inserts the alcoves to the Alcove_vertex_list
  struct Alcove_visitor {
    Alcove_visitor(typename Point_range::const_iterator& p_it,
                   Alcove_vertex_graph& av_graph,
                   bool& store_points)
      : p_it_(p_it), av_graph_(av_graph), store_points_(store_points) {}
    
    void operator() (const Filtered_alcove& a) {
      auto& a_map = av_graph_.a_map;
      auto a_it = a_map.find(a.id);
      if (a_it == a_map.end())
        av_graph_.emplace_alcove(a);
      else
        a_it->second.second = std::min(a_it->second.second, a.f);
    }
  private :
    typename Point_range::const_iterator& p_it_;
    Alcove_vertex_graph& av_graph_;
    bool& store_points_;
  };

  // The procedure that emplaces the alcoves in the map
  void compute_a_map(const Point_range& point_vector, double init_level, double eps, bool store_points) {
    for (auto p_it = point_vector.begin(); p_it != point_vector.end(); ++p_it)
      cs_.alcoves_of_ball(*p_it, init_level, eps, Alcove_visitor(p_it, av_graph_, store_points));
    for (auto m: av_graph_.a_map)
      m.second.second = std::sqrt(m.second.second);
    // typedef typename Alcove_map::iterator Alcove_iterator;
    // std::vector<Alcove_iterator> iterators;
    // for (auto m_it = a_map.begin(); m_it != a_map.end(); ++m_it)
    //   iterators.push_back(m_it);
    // struct Filt_compare {
    //   bool operator() (const Alcove_iterator& lhs, const Alcove_iterator& rhs) const {
    //     return std::get<3>(lhs->second) < std::get<3>(rhs->second);
    //   }
    // };
    // std::sort(iterators.begin(), iterators.end(), Filt_compare());
    // int k = 0;
    // for (auto it: iterators)
    //   std::get<0>(it->second) = k++;
  }

  
public:

  Coxeter_complex(const Point_range& point_vector,
                  const Coxeter_system& cs,
                  double init_level=1,
                  double eps=0,
                  bool store_points = false) : cs_(cs)
  {
    compute_a_map(point_vector, init_level, eps, store_points);
#ifdef DEBUG_TRACES
    auto& a_map = av_graph_.a_map;
    std::cout << "Alcove map size is " << a_map.size() << ".\n";
    
    std::cout << "AMap:\n";
    for (auto m: a_map) 
      std::cout << m.first << ": filt=" << m.second.second << std::endl;    
    std::cout << "\n";
#endif
  
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
