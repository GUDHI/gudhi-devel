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

#include <gudhi/SparseMsMatrix.h>
#include <gudhi/Fake_simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Coxeter_complex/Collapse.h>
#include <gudhi/Coxeter_complex/Hasse_diagram.h>

#include <gudhi/Simplex_tree.h>
#include <limits>
#include <iomanip>
#include <gudhi/Coxeter_complex/Cech_blocker.h>
#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/distance_functions.h>
// #include <gudhi/Coxeter_complex/Simplex_with_cofaces.h>


#include <boost/iterator/iterator_facade.hpp>

double num_error = 1.0e-10;

namespace Gudhi {

std::ostream& operator<<(std::ostream & os, const typename Simple_coxeter_system::Alcove_id& a_id) {
  int i = 0, j = 0;
  os << "[";
  if (a_id.empty())
    return os;
  auto a_it = a_id.begin();
  if (a_it == a_id.end() - 1)
    os << "\033[1;32m" << *a_it++ << "\033[0m";
  else if (j == i)
    os << "\033[1;31m" << *a_it++ << "\033[0m";
  else
    os << *a_it++;
  if (j == 0) {
    i++;
    j = i;
  }
  else
    j--;
  while (a_it != a_id.end()) {
    if (a_it == a_id.end() - 1)
      os << ", \033[1;32m" << *a_it++ << "\033[0m";
    else if (j == i)
      os << ", \033[1;31m" << *a_it++ << "\033[0m";
    else
      os << ", " << *a_it++;
    if (j == 0) {
      i++;
      j = i;
    }
    else
      j--;
  }
  std::cout << "]";
  return os;
}
  
template <class Point_range,
          class Coxeter_system>
class Coxeter_complex {
  
public:
  
  using Alcove_id = typename Coxeter_system::Alcove_id;
  using Vertex_id = Alcove_id;    

  using Index_range = std::vector<std::size_t>;
  using Vertex_map = std::map<Vertex_id, Index_range>;

  using Point_pointers = std::vector<typename Point_range::const_iterator>;
  using Vertex_pointers = std::vector<typename Vertex_map::iterator>;
  using Alcove = std::tuple<std::size_t, Point_pointers, Vertex_pointers>;
  using Alcove_map = std::map<Alcove_id, Alcove>;

  using Simplex_tree = Gudhi::Simplex_tree<>;
  using Proximity_graph = Gudhi::Proximity_graph<Simplex_tree>;
  
  struct Pointer_compare {
    typedef typename Vertex_map::iterator Pointer;
    bool operator()(const Pointer& lhs, const Pointer& rhs) const { 
      return lhs->first < rhs->first;
    }
  };
  using Vertex_index_map = std::map<typename Vertex_map::iterator, int, Pointer_compare>;

  struct Alcove_pair_compare {
    typedef std::pair<Alcove_id, double> Pair;
    bool operator()(const Pair& lhs, const Pair& rhs) const { 
      return lhs.second < rhs.second || lhs.first < lhs.first;
    }
  };
  using Clique_id = int;
  using Mask = std::map<Alcove_id, std::vector<Clique_id> >;
  using Mask_pointers = std::vector<typename Mask::iterator>;
  struct Double_compare {
    bool operator()(const double& lhs, const double& rhs) const { 
      return lhs + num_error < rhs;
    }
  };
  using Dist_mask_map = std::map<double, Mask_pointers, Double_compare>;
  // using Mask = std::pair<Mask_clique_map, Dist_mask_map> ;
  using Cech_nerve = std::vector<std::pair<std::size_t, Clique_id> >;

  using Filtered_alcove = std::tuple<std::size_t, double, Vertex_pointers>;
  using Filtered_alcove_map = std::map<Alcove_id, Filtered_alcove>;
  
  const Coxeter_system& cs_;
  
  Alcove_map a_map;
  Vertex_map v_map;
  std::size_t max_id;
  Vertex_index_map vi_map;

  class Alcove_iterator : public boost::iterator_facade< Alcove_iterator,
                                                         std::vector<std::size_t> const,
                                                         boost::forward_traversal_tag> {
  private:
    typename Alcove_map::const_iterator it_;
    Alcove_map const& a_map_;
    Vertex_index_map const& vi_map_;
    std::vector<std::size_t> value_; 
    friend class boost::iterator_core_access;

    void update_value() {
      value_.clear();
      if (it_ != a_map_.end()) {
        for (auto m_it: std::get<2>(it_->second))
          value_.push_back(vi_map_.at(m_it));
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
    Alcove_iterator(typename Alcove_map::const_iterator it,
                    Alcove_map const& a_map,
                    Vertex_index_map const& vi_map)
      : it_(it), a_map_(a_map), vi_map_(vi_map) {
      update_value();
    }

    int dimension() {
      return value_.size()-1;
    }
  };

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

private:

  struct Alcove_visitor {
    Alcove_visitor(typename Point_range::const_iterator& p_it,
                   Alcove_map& a_map,
                   bool& store_points)
      : p_it_(p_it), a_map_(a_map), store_points_(store_points) {}
    
    void operator() (const Alcove_id& a) {
      auto a_it = a_map_.find(a);
      if (a_it == a_map_.end()) {
        if (store_points_)
          a_map_.emplace(a, std::make_tuple(a_map_.size(), Point_pointers(1, p_it_), Vertex_pointers()));
        else
          a_map_.emplace(a, std::make_tuple(a_map_.size(), Point_pointers(), Vertex_pointers()));
      }
      else
        if (store_points_)
          std::get<1>(a_it->second).push_back(p_it_);
    }
  private :
    typename Point_range::const_iterator& p_it_;
    Alcove_map& a_map_;
    bool& store_points_;
  };
  
  void compute_a_map(const Point_range& point_vector, double init_level, double eps, bool store_points) {
    for (auto p_it = point_vector.begin(); p_it != point_vector.end(); ++p_it)
      cs_.alcoves_of_ball(*p_it, init_level, eps, Alcove_visitor(p_it, a_map, store_points));
  }
  
  void compute_v_map() {
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
  }
  
public:
  
  Coxeter_complex(const Point_range& point_vector, const Coxeter_system& cs, double init_level=1, double eps=0, bool store_points = false)
    : cs_(cs), max_id(0) {
    clock_t start, end, global_start;
    double time;
    global_start = clock();
    start = clock();
    compute_a_map(point_vector, init_level, eps, store_points);
    end = clock();
    time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    std::cout << "Computed alcove coordinate map in " << time << " s. \n";
    std::cout << "Alcove map size is " << a_map.size() << ".\n";
    
    // std::cout << "AMap:\n";
    // for (auto m: a_map) 
    //   std::cout << m.first << ": " << std::get<0>(m.second) << ", "
    //             << "size=" << std::get<1>(m.second).size() << std::endl;    
    // std::cout << "\n";
    
    start = clock();
    compute_v_map();
    end = clock();
    time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    std::cout << "Computed vertex map in " << time << " s. \n";
    std::cout << "Vertex map size is " << v_map.size() << ".\n";
    end = clock();
    time = static_cast<double>(end - global_start) / CLOCKS_PER_SEC;
    std::cout << "Total time: " << time << " s. \n";

    // std::cout << "VMap:\n";
    // for (auto m: v_map) 
    //   std::cout << m.first << ": " << m.second << std::endl;
    // std::cout << "\n";
    
    std::size_t max_dim = 0; 
    for (auto m: v_map) {
      if (m.second.size()-1 > max_dim)
        max_dim = m.second.size()-1;
    }
    std::cout << "Dimension of the nerve is " << max_dim << ".\n\n";    

    unsigned index = 0;
    for (auto v_it = v_map.begin(); v_it != v_map.end(); ++v_it, ++index)
      vi_map.emplace(v_it, index);
                  
    
    // subdivision part
        
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

  // Coxeter_complex(const std::string& input_file, const Coxeter_system& cs, int init_level=1, bool store_points = false)
  //   : cs_(cs), max_id(0) {
  //   clock_t start, end, global_start;
  //   double time;
  //   global_start = clock();
  //   start = clock();
  //   Off_point_range<typename Point_range::value_type> off_range(input_file);
  //   for (auto p: off_range) {
  //     Alcove_id s_id = cs.alcove_coordinates(p, init_level); 
  //     auto a_it = a_map.find(s_id);
  //     if (a_it == a_map.end())
  //       a_map.emplace(s_id, std::make_tuple(max_id++, Point_pointers(), Vertex_pointers()));
  //     else
  //       std::get<1>(a_it->second).push_back(p);
  //   }
  //   end = clock();
  //   time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
  //   std::cout << "Computed alcove coordinate map in " << time << " s. \n";
  //   std::cout << "Alcove map size is " << a_map.size() << ".\n";
    
  //   // std::cout << "AMap:\n";
  //   // for (auto m: a_map) 
  //   //   std::cout << m.first << ": " << std::get<0>(m.second) << ", "
  //   //             << "size=" << std::get<1>(m.second).size() << std::endl;    
  //   // std::cout << "\n";
    
  //   start = clock();
  //   compute_v_map();
  //   end = clock();
  //   time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
  //   std::cout << "Computed vertex map in " << time << " s. \n";
  //   std::cout << "Vertex map size is " << v_map.size() << ".\n";
  //   end = clock();
  //   time = static_cast<double>(end - global_start) / CLOCKS_PER_SEC;
  //   std::cout << "Total time: " << time << " s. \n";

  //   // std::cout << "VMap:\n";
  //   // for (auto m: v_map) 
  //   //   std::cout << m.first << ": " << m.second << std::endl;
  //   // std::cout << "\n";
    
  //   std::size_t max_dim = 0; 
  //   for (auto m: v_map) {
  //     if (m.second.size()-1 > max_dim)
  //       max_dim = m.second.size()-1;
  //   }
  //   std::cout << "Dimension of the nerve is " << max_dim << ".\n\n";    

  //   unsigned index = 0;
  //   for (auto v_it = v_map.begin(); v_it != v_map.end(); ++v_it, ++index)
  //     vi_map.emplace(v_it, index);
  // }

  /* Returns true if and only if the simplex is critical.
   */
  template <class SimplexTree>
  bool is_maximal(typename SimplexTree::Simplex_handle sh,
                  SimplexTree& st)
  {
    assert(st.num_simplices() != 0);
    return st.cofaces_simplex_range(sh,1).empty();
  }

private:
  struct Alcove_mask_visitor {
    Alcove_mask_visitor(Mask& mask,
                        Dist_mask_map& dmm,
                        const Coxeter_system& cs,
                        std::vector<double>& barycenter,
                        double rad)
      : mask_(mask), dmm_(dmm), cs_(cs), barycenter_(barycenter), rad_(rad) {}
    
    void operator() (const Alcove_id& a) {
      std::vector<double> bary_a = cs_.barycenter(a);
      double dist = Euclidean_distance()(barycenter_, bary_a);
      if (dist <= 2*rad_+1.0e-10) {
        auto r_pair = mask_.emplace(a, std::vector<Clique_id>());
        auto d_it = dmm_.find(dist/rad_);
        if (d_it == dmm_.end())
          dmm_.emplace(dist/rad_, Mask_pointers(1, r_pair.first));
        else
          d_it->second.push_back(r_pair.first);
      }
    }
  private :
    Mask& mask_;
    Dist_mask_map& dmm_;
    const Coxeter_system& cs_;
    std::vector<double>& barycenter_;
    double rad_;
  };
  
public:
  /* Should never be called if a_map or v_map are empty */
  void build_mask(Mask& mask) {
    const unsigned d = v_map.begin()->first.size();
    double init_level = a_map.begin()->first.level();
    Alcove_id base(init_level);
    for (unsigned i = 0; i < (d*(d+1))/2; i++)
      base.push_back(0);
    std::vector<double> barycenter = cs_.barycenter(base);
    std::cout << barycenter << "\n";
    double rad = std::sqrt(d*(d+2.)/(d+1)/12)/init_level;
    Dist_mask_map dmm;
    cs_.alcoves_of_ball(barycenter,
                        init_level,
                        2*rad,
                        Alcove_mask_visitor(mask, dmm, cs_, barycenter, rad));
    for (auto m: mask)
      std::cout << m.first << " " << m.second << std::endl;
    std::cout << mask.size() << std::endl;
    std::cout << std::endl;
    std::cout << std::fixed;
    std::cout << std::setprecision(6);
    for (auto d: dmm) {
      std::cout << d.first << " " << d.second.size() << " [";
      if (!d.second.empty()) {
        auto p_it = d.second.begin();
        std::cout << (*p_it++)->first;
        while (p_it != d.second.end())
          std::cout << ", " << (*p_it++)->first;
      }
      std::cout << "]\n";
    }
  }
  
  void construct_clique_complex() {
    Mask mask;
    build_mask(mask);
  }
  
  template <class AMap_iterator>
  void subdivide_alcove(AMap_iterator a_it) {
    for (auto v_it: std::get<2>(a_it->second)) {
      auto find_it = std::find(v_it->second.begin(), v_it->second.end(), std::get<0>(a_it->second));
      v_it->second.erase(find_it);
      if (v_it->second.empty())
        v_map.erase(v_it);
    }
    for (auto p_it: std::get<1>(a_it->second)) {
      Alcove_id s_id = cs_.alcove_coordinates(*p_it, 2 * a_it->first.level()); 
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
            if (v.level() == new_a_it->first.level()) {
              auto m_it = a_map.begin();
              while (m_it != a_map.end() && m_it->first.level() < new_a_it->first.level()) {
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

  template <class Simplex_range>
  void write_mesh(Simplex_range& range, std::string file_name) const {
    cs_.write_mesh(v_map, range, file_name);
  }

  void write_mesh(std::string file_name = "toplex.mesh") const {
    typedef boost::iterator_range<Alcove_iterator> Max_simplex_range;
    Max_simplex_range range(Alcove_iterator(a_map.begin(), a_map, vi_map),
                            Alcove_iterator(a_map.end(), a_map, vi_map));
    cs_.write_mesh(v_map, range, file_name);
  }

  void strong_collapse(bool pers_out = true) {
    struct Size_compare {
      bool operator() (const Vertex_pointers& lhs, const Vertex_pointers& rhs) const {
        return lhs->second.size() < rhs->second.size();
      }
    };
    bool stable = false;
    while (!stable) {
      stable = true;
      for (auto m: a_map) {
        std::sort(std::get<2>(m.second).begin(), std::get<2>(m.second).end(), Size_compare());
      }
    }
  }
  
  void collapse(bool pers_out = true) {
    Fake_simplex_tree stree;

    clock_t start, end;
    double time;
        
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

    /* Full simplicial complex computation */
    // if (pers_out) {
    //   Simplex_tree<> full_stree;
    //   for (auto s: max_simplex_range)
    //     full_stree.insert_simplex_and_subfaces(s);
    //   std::cout << "Start persistence calculation..." << std::endl;
    //   using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
    //   using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree<>, Field_Zp>;
    //   double min_persistence = 0;
    //   // Sort the simplices in the order of the filtration
    //   full_stree.initialize_filtration();
    //   // Compute the persistence diagram of the complex
    //   Persistent_cohomology pcoh(full_stree);
    //   // initializes the coefficient field for homology
    //   pcoh.init_coefficients(3);
      
    //   pcoh.compute_persistent_cohomology(min_persistence);
    //   std::ofstream out("persdiag_cox_full.out");
    //   pcoh.output_diagram(out);
    //   out.close();
    //   std::cout << "Persistence complete." << std::endl;
    //   int chi = 0;
    //   for (auto sh: full_stree.complex_simplex_range())
    //     chi += 1-2*(full_stree.dimension(sh)%2);
    //   std::cout << "Euler characteristic of full_stree is " << chi << std::endl;
    // }

    start = clock();

    // SparseMsMatrix mat(a_map.size(), max_simplex_range);
                                
    // // auto matrix_formed  = std::chrono::high_resolution_clock::now();
    // std::cout << "Start strong collapse..." << std::endl;
    // mat.strong_collapse();
    // std::vector<std::vector<std::size_t>> scoll_simplex_range;
    // mat.output_simplices(std::back_inserter(scoll_simplex_range));

    // struct Size_comparison {
    //   using Container = std::vector<std::size_t>;
    //   bool operator() (const Container& lhs, const Container& rhs) {
    //     return lhs.size() < rhs.size();
    //   }
    // };
    // std::sort(scoll_simplex_range.begin(), scoll_simplex_range.end(), Size_comparison());
    
    Simplex_tree coll_stree;
    Collapse coll(max_simplex_range, coll_stree);
    end = clock();
    time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    std::cout << "Number of all simplices after collapse: " << coll_stree.num_simplices() << "\n";
    std::cout << "Collapse took " << time << " s. \n";
    // std::cout << coll_stree << "\n";
    int dim_complex = 0;
    for (auto sh: coll_stree.complex_simplex_range()) {
      if (coll_stree.dimension(sh) > dim_complex)
        dim_complex = coll_stree.dimension(sh);
    }
    std::cout << "Dimension of the collapsed complex is " << dim_complex << "\n";
    // int chi = 0;
    // for (auto sh: coll_stree.complex_simplex_range())
    //   chi += 1-2*(coll_stree.dimension(sh)%2);
    // std::cout << "Euler characteristic of coll_stree is " << chi << std::endl;
    
    if (pers_out) {
      std::cout << "Start persistence calculation..." << std::endl;
      int chi = 0;
      for (auto sh: coll_stree.complex_simplex_range())
        chi += 1-2*(coll_stree.dimension(sh)%2);
      std::cout << "Euler characteristic of coll_stree is " << chi << std::endl;

      using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
      using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp>;
      double min_persistence = 0;
      // Sort the simplices in the order of the filtration
      coll_stree.set_dimension(dim_complex+1);
      coll_stree.initialize_filtration();
      // Compute the persistence diagram of the complex
      Persistent_cohomology pcoh(coll_stree);
      // initializes the coefficient field for homology
      pcoh.init_coefficients(3);
      
      pcoh.compute_persistent_cohomology(min_persistence);
      std::cout << pcoh.persistent_betti_numbers(0,0) << std::endl;
      std::ofstream out("persdiag_cox.out");
      pcoh.output_diagram(out);
      out.close();
      std::cout << "Persistence complete." << std::endl;
      typedef Simplex_tree_simplex_iterator<Simplex_tree> Simplex_tree_iterator;
      typedef boost::iterator_range<Simplex_tree_iterator> Simplex_tree_range;
      Simplex_tree_range simplex_tree_range(Simplex_tree_iterator(coll_stree.complex_simplex_range().begin(), coll_stree),
                                            Simplex_tree_iterator(coll_stree.complex_simplex_range().end(), coll_stree));
      write_mesh(simplex_tree_range, "coll_stree.mesh");
      
    //   witness_complex::Dim_lists<Simplex_tree<>> simplices(coll_stree, collDim, 1);
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
    }
    
    
  }

  // struct Filtered_alcove_visitor {
  //   Alcove_visitor(typename Point_range::const_iterator& p_it,
  //                  Filtered_alcove_map& a_map,
  //                  bool& store_points)
  //     : p_it_(p_it), a_map_(a_map), store_points_(store_points) {}
    
  //   void operator() (const Alcove_id& a) {
  //     auto a_it = a_map_.find(a);
  //     if (a_it == a_map_.end()) {
  //       if (store_points_)
  //         a_map_.emplace(a, std::make_tuple(a_map_.size(), Point_pointers(1, p_it_), Vertex_pointers()));
  //       else
  //         a_map_.emplace(a, std::make_tuple(a_map_.size(), Point_pointers(), Vertex_pointers()));
  //     }
  //     else
  //       if (store_points_)
  //         std::get<1>(a_it->second).push_back(p_it_);
  //   }
  // private :
  //   typename Point_range::const_iterator& p_it_;
  //   Filtered_alcove_map& a_map_;
  //   bool& store_points_;
  // };
  
  // void filtered_complex(const Point_range& point_vector, double init_level, double eps) {
  //   typedef typename Hasse_diagram::Node_iterator Node_iterator;
  //   Hasse_diagram hasse;
  //   std::map<std::size_t, Node_iterator> iv_map;
  //   Filtered_alcove_map filt_a_map;

  //   for (auto p_it = point_vector.begin(); p_it != point_vector.end(); ++p_it)
  //     cs_.alcoves_of_ball(*p_it, init_level, eps, Filtered_alcove_visitor(p_it, filt_a_map, false));
  // }
  
};

}
#endif
