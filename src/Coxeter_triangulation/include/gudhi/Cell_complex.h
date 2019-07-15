/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2019 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CELL_COMPLEX_H_
#define CELL_COMPLEX_H_

#include <Eigen/Dense>

#include <gudhi/Permutahedral_representation/Simplex_comparator.h>
#include <fstream> // for Hasse_diagram_persistence.h
#include <gudhi/Hasse_diagram_persistence.h> // for Hasse_cel

namespace Gudhi {

namespace coxeter_triangulation {


/** \class Cell_complex 
 *  \brief A class that assembles methods for manifold tracing algorithm.
 *
 *  \tparam Out_simplex_map_ The type of a map from a simplex type that is a
 *   model of SimplexInCoxeterTriangulation to Eigen::VectorXd.
 *
 *  \ingroup coxeter_triangulation
 */
template <class Out_simplex_map_>
class Cell_complex {
public:
  typedef typename Out_simplex_map_::key_type Simplex_handle;
  using Hasse_cell = Gudhi::Hasse_diagram::Hasse_diagram_cell<int, double, double>;
  typedef std::map<Simplex_handle,
		   Hasse_cell*,
		   Simplex_comparator<Simplex_handle> > Simplex_cell_map;
  typedef std::vector<Simplex_cell_map> Simplex_cell_maps;
  
  typedef std::map<Hasse_cell*, Simplex_handle> Cell_simplex_map;

private:
  Hasse_cell* insert_cell(const Simplex_handle& simplex, std::size_t cell_d) {
    Simplex_cell_map& simplex_cell_map = simplex_cell_maps[cell_d];
    auto curr_it = simplex_cell_map.lower_bound(simplex.smallest_coface());
    auto upper_bound = simplex_cell_map.upper_bound(simplex);
    for (; curr_it != upper_bound; ++curr_it) 
      if (simplex.is_face_of(curr_it->first)) {
	std::cout << "Current simplex \033[1;31m" << simplex
		  << "\033[0m is a face of \033[1;34m" << curr_it->first << "\033[0m\n";
	Hasse_cell* cell = curr_it->second;
	cell_simplex_map.at(cell) = simplex;
	simplex_cell_map.erase(curr_it);
	simplex_cell_map.emplace(std::make_pair(simplex, cell));
	for (auto& sc_pair: simplex_cell_map)
	  std::cout << " \033[1;31m" << sc_pair.first << "\033[0m\n";
	return cell;
      }
    upper_bound = simplex_cell_map.upper_bound(simplex.greatest_face());
    for (; curr_it != upper_bound; ++curr_it) 
      if (curr_it->first.is_face_of(simplex)) {
	std::cout << "Current simplex \033[1;34m" << simplex
		  << "\033[0m is a coface of \033[1;31m" << curr_it->first << "\033[0m\n";
	for (auto& sc_pair: simplex_cell_map)
	  std::cout << " \033[1;31m" << sc_pair.first << "\033[0m\n";
	return curr_it->second;
      }
    hasse_cells.push_back(new Hasse_cell(cell_d, 0.0));
    Hasse_cell* new_cell = hasse_cells.back();
    cell_simplex_map.emplace(std::make_pair(new_cell, simplex));
    simplex_cell_map.emplace(std::make_pair(simplex, new_cell));
    std::cout << "Current simplex \033[1;31m" << simplex << "\n";
    for (auto& sc_pair: simplex_cell_map)
      std::cout << " \033[1;31m" << sc_pair.first << "\033[0m\n";

    return new_cell;
  }

  void expand_level(std::size_t cell_d) {
  }

  void construct_complex_(const Out_simplex_map_& out_simplex_map) {
    for (auto& os_pair: out_simplex_map)
      insert_cell(os_pair.first, 0);    
  }

  void construct_complex_(const Out_simplex_map_& interior_simplex_map,
			  const Out_simplex_map_& boundary_simplex_map) {
  }
  
public:  

  Cell_complex(std::size_t intrinsic_dimension)
    : intr_d_(intrinsic_dimension) {}
  
  void construct_complex(const Out_simplex_map_& out_simplex_map) {
    simplex_cell_maps.resize(intr_d_ + 1);
    construct_complex_(out_simplex_map);
  }
  
  void construct_complex(const Out_simplex_map_& out_simplex_map,
			 std::size_t limit_dimension) {
    simplex_cell_maps.resize(limit_dimension + 1);
    construct_complex_(out_simplex_map);
  }

  void construct_complex(const Out_simplex_map_& interior_simplex_map,
			 const Out_simplex_map_& boundary_simplex_map) {
    construct_complex_(interior_simplex_map, boundary_simplex_map);
  }

  void construct_complex(const Out_simplex_map_& interior_simplex_map,
			 const Out_simplex_map_& boundary_simplex_map,
			 std::size_t limit_d) {
    simplex_cell_maps.resize(limit_d + 1);
    construct_complex_(interior_simplex_map, boundary_simplex_map);
  }

  const Simplex_cell_map& simplex_cell_map(std::size_t cell_d) const {
    return simplex_cell_maps[cell_d];
  }
  
private:
  std::size_t intr_d_;
  Simplex_cell_maps simplex_cell_maps;
  Cell_simplex_map cell_simplex_map;
  std::vector<Hasse_cell*> hasse_cells;
  
};

} // namespace coxeter_triangulation 

} // namespace Gudhi

#endif
