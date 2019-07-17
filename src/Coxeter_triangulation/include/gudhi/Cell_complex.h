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

#include <algorithm>
#include <gudhi/Permutahedral_representation/Simplex_comparator.h>
#include <fstream> // for Hasse_diagram_persistence.h

#include <gudhi/Hasse_diagram_persistence.h> // for Hasse_cell
#include <gudhi/Permutahedral_representation/join.h>

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
  typedef Gudhi::Hasse_diagram::Hasse_diagram_cell<int, double, double> Hasse_cell;
  typedef std::map<Simplex_handle,
		   Hasse_cell*,
		   Simplex_comparator<Simplex_handle> > Simplex_cell_map;
  typedef std::vector<Simplex_cell_map> Simplex_cell_maps;
  
  typedef std::map<Hasse_cell*, Simplex_handle> Cell_simplex_map;

  typedef std::map<Hasse_cell*, Eigen::VectorXd> Cell_point_map;

private:
  
  enum class Result_type {self, face, coface, inserted};
  
  struct Insert_result {
    Hasse_cell* cell;
    bool success;
  };
  
  std::pair<Hasse_cell*, Result_type> insert_cell(const Simplex_handle& simplex, std::size_t cell_d) {
    // std::cout << "Insert simplex for " << simplex << "\n";
    Simplex_cell_map& simplex_cell_map = simplex_cell_maps_[cell_d];
    auto curr_it = simplex_cell_map.lower_bound(simplex.smallest_coface());
    auto upper_bound = simplex_cell_map.upper_bound(simplex);
    for (; curr_it != upper_bound; ++curr_it) {
      // std::cout << " Checking if face: " << curr_it->first << "\n";
      if (simplex == curr_it->first) {
	// std::cout << "  Same simplex!\n";
	return std::make_pair(curr_it->second, Result_type::self);
      }
      if (simplex.is_face_of(curr_it->first)) {
	// std::cout << "  Face!\n";
	Hasse_cell* cell = curr_it->second;
	cell_simplex_map_.at(cell) = simplex;
	simplex_cell_map.erase(curr_it++);
	while (curr_it != upper_bound) {
	  if (simplex == curr_it->first)
	    continue;
	  if (simplex.is_face_of(curr_it->first)) {
	    // std::cout << "  Post-deleting " << curr_it->first << "\n";
	    simplex_cell_map.erase(curr_it++);
	  }
	  else
	    curr_it++;
	}
	simplex_cell_map.emplace(std::make_pair(simplex, cell));
        return std::make_pair(cell, Result_type::face);
      }
    }
    upper_bound = simplex_cell_map.upper_bound(simplex.greatest_face());
    for (; curr_it != upper_bound; ++curr_it) {
      // std::cout << " Checking if coface: " << curr_it->first << "\n";
      if (curr_it->first.is_face_of(simplex)) {
	// std::cout << "  Coface!\n";
	return std::make_pair(curr_it->second, Result_type::coface);
      }
    }
    hasse_cells_.push_back(new Hasse_cell(cell_d, 0.0));
    Hasse_cell* new_cell = hasse_cells_.back();
    cell_simplex_map_.emplace(std::make_pair(new_cell, simplex));
    simplex_cell_map.emplace(std::make_pair(simplex, new_cell));
    // std::cout << " OK for insertion!\n";
    return std::make_pair(new_cell, Result_type::inserted);
  }

  void expand_level(std::size_t cell_d) {
    std::size_t min_d = cell_d + intr_d_;
    std::size_t amb_d = simplex_cell_maps_[0].begin()->first.vertex().size();
    for (auto& sc_pair: simplex_cell_maps_[cell_d - 1])
      min_d = std::min(min_d, sc_pair.first.dimension());
    for (std::size_t i = min_d + 1; i <= std::min(cell_d + intr_d_, amb_d); ++i) {
      for (auto& sc_pair: simplex_cell_maps_[cell_d - 1]) {
	const Simplex_handle& simplex = sc_pair.first;
	if (simplex.dimension() >= i)
	  continue;
        Hasse_cell* cell = sc_pair.second;
	for (Simplex_handle coface: simplex.coface_range(i)) {
	  Hasse_cell* new_cell;
	  Result_type success;
	  std::tie(new_cell, success) = insert_cell(coface, cell_d);
	  if (success == Result_type::self || success == Result_type::inserted) {
	    new_cell->get_boundary().emplace_back(std::make_pair(cell, 1));
	    // std::cout << "The cell " << cell_simplex_map_[cell]
	    // 	      << " (" << cell <<  ") became a face of the cell "
	    // 	      << cell_simplex_map_[new_cell] << " (" << new_cell << ")\n";
	  }
	}
      }
      // std::cout << "\nThe layer " << cell_d << " before the cleanup at simplex dimension " << i << "\n";
      // for (auto& sc_pair: simplex_cell_maps_[cell_d]) {
      // 	std::cout << " \033[1;31m" << sc_pair.first << "\033[0m: [ ";
      // 	if (sc_pair.second->get_boundary().empty()) {
      // 	  std::cout << "]";
      // 	  continue;
      // 	}
      // 	Hasse_cell* cell = sc_pair.second;
      // 	auto b_it = cell->get_boundary().begin();
      // 	std::cout << "\033[1;32m" << cell_simplex_map_[(b_it++)->first] << "\033[0m";
      // 	for (; b_it != cell->get_boundary().end(); ++b_it) {
      // 	  std::cout << ",  \033[1;32m" << cell_simplex_map_[b_it->first] << "\033[0m";
      // 	}
      // 	std::cout << " ]";
      // 	std::vector<Simplex_handle> faces;
      // 	faces.reserve(cell->get_boundary().size());
      // 	for (auto& bi_pair: cell->get_boundary()) {
      // 	  const Simplex_handle& face = cell_simplex_map_.at(bi_pair.first);
      // 	  faces.push_back(face);
      // 	}
      // 	Simplex_handle join = Gudhi::coxeter_triangulation::join(faces);
      // 	std::cout << " join = \033[1;33m" << join << "\033[0m\n";
      // }
      Simplex_cell_map& new_layer = simplex_cell_maps_[cell_d];
      auto curr_it = new_layer.begin();
      while (curr_it != new_layer.end()) {
	std::vector<Simplex_handle> faces;
	const Simplex_handle& simplex = curr_it->first;
	Hasse_cell* cell = curr_it->second;
	if (cell->get_boundary().size() == 1) {
	  new_layer.erase(curr_it++);
	  continue;
	}
	faces.reserve(cell->get_boundary().size());
	for (auto& bi_pair: cell->get_boundary()) {
	  const Simplex_handle& face = cell_simplex_map_.at(bi_pair.first);
	  faces.push_back(face);
	}
	Simplex_handle join = Gudhi::coxeter_triangulation::join(faces);
	bool join_is_face = false;
	for (auto& bi_pair: cell->get_boundary()) {
	  const Simplex_handle& face = cell_simplex_map_.at(bi_pair.first);
	  if (face == join) {
	    join_is_face = true;
	    break;
	  }
	}
	if (join_is_face)
	  new_layer.erase(curr_it++);
	else if (simplex != join) {
	  new_layer.erase(curr_it++);
	  insert_cell(join, cell_d);
	}
	else
	  curr_it++;
      }
      // std::cout << "\nThe layer " << cell_d << " after the cleanup at simplex dimension " << i << "\n";
      // for (auto& sc_pair: simplex_cell_maps_[cell_d]) {
      // 	std::cout << " \033[1;31m" << sc_pair.first << "\033[0m: [ ";
      // 	if (sc_pair.second->get_boundary().empty()) {
      // 	  std::cout << "]";
      // 	  continue;
      // 	}
      // 	auto b_it = sc_pair.second->get_boundary().begin();
      // 	std::cout << "\033[1;32m" << cell_simplex_map_[(b_it++)->first] << "\033[0m";
      // 	for (; b_it != sc_pair.second->get_boundary().end(); ++b_it) {
      // 	  std::cout << ",  \033[1;32m" << cell_simplex_map_[b_it->first] << "\033[0m";
      // 	}
      // 	std::cout << " ]\n";
      // }
    }
  }

  void construct_complex_(const Out_simplex_map_& out_simplex_map) {
    for (auto& os_pair: out_simplex_map) {
      Hasse_cell* new_cell;
      Result_type success;	  
      std::tie(new_cell, success) = insert_cell(os_pair.first, 0);
      if (success == Result_type::inserted)
	cell_point_map_.emplace(std::make_pair(new_cell, os_pair.second));
      else if (success == Result_type::face)
	cell_point_map_.at(new_cell) = os_pair.second;	
    }
    // std::cout << "Finished building the layer 0. Simplices:\n";
    // std::size_t color = 31;
    // for (auto& sc_pair: simplex_cell_maps_[0])
    //   std::cout << "\033[1;" << color << "m" << sc_pair.first << "\033[0m\n";
    // std::cout << "Size: " << simplex_cell_maps_[0].size() << "\n";
    // std::cout << "\n";
    // color++;
    for (std::size_t cell_d = 1;
	 cell_d < simplex_cell_maps_.size() && !simplex_cell_maps_[cell_d - 1].empty();
	 ++cell_d) {
      expand_level(cell_d);
      // ++color;
      // std::cout << "Finished building the layer " << cell_d << ". Simplices:\n";
      // for (auto& sc_pair: simplex_cell_maps_[cell_d])
      // 	std::cout << "\033[1;" << color << "m" << sc_pair.first << "\033[0m\n";
      // std::cout << "Size: " << simplex_cell_maps_[cell_d].size() << "\n";
      // std::cout << "\n";
    }
  }
  
  void construct_complex_(const Out_simplex_map_& interior_simplex_map,
			  const Out_simplex_map_& boundary_simplex_map) {
  }
  
public:  

  Cell_complex(std::size_t intrinsic_dimension)
    : intr_d_(intrinsic_dimension) {}
  
  void construct_complex(const Out_simplex_map_& out_simplex_map) {
    simplex_cell_maps_.resize(intr_d_ + 1);
    construct_complex_(out_simplex_map);
  }
  
  void construct_complex(const Out_simplex_map_& out_simplex_map,
			 std::size_t limit_dimension) {
    simplex_cell_maps_.resize(limit_dimension + 1);
    construct_complex_(out_simplex_map);
  }

  void construct_complex(const Out_simplex_map_& interior_simplex_map,
			 const Out_simplex_map_& boundary_simplex_map) {
    construct_complex_(interior_simplex_map, boundary_simplex_map);
  }

  void construct_complex(const Out_simplex_map_& interior_simplex_map,
			 const Out_simplex_map_& boundary_simplex_map,
			 std::size_t limit_d) {
    simplex_cell_maps_.resize(limit_d + 1);
    construct_complex_(interior_simplex_map, boundary_simplex_map);
  }

  const Simplex_cell_map& simplex_cell_map(std::size_t cell_d) const {
    return simplex_cell_maps_[cell_d];
  }

  const Cell_point_map& cell_point_map() const {
    return cell_point_map_;
  }

private:
  std::size_t intr_d_;
  Simplex_cell_maps simplex_cell_maps_;
  Cell_simplex_map cell_simplex_map_;
  Cell_point_map cell_point_map_;
  std::vector<Hasse_cell*> hasse_cells_;
};

} // namespace coxeter_triangulation 

} // namespace Gudhi

#endif
