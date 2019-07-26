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
  typedef Gudhi::Hasse_diagram::Hasse_diagram_cell<int, double, bool> Hasse_cell;
  typedef std::map<Simplex_handle,
		   Hasse_cell*,
		   Simplex_comparator<Simplex_handle> > Simplex_cell_map;
  typedef std::vector<Simplex_cell_map> Simplex_cell_maps;
  
  typedef std::map<Hasse_cell*, Simplex_handle> Cell_simplex_map;

  typedef std::map<Hasse_cell*, Eigen::VectorXd> Cell_point_map;

  typedef std::map<Hasse_cell*, bool> Cell_bool_map;

private:
  
  enum class Result_type {self, face, coface, inserted};
  
  struct Insert_result {
    Hasse_cell* cell;
    bool success;
  };
  
  std::pair<Hasse_cell*, Result_type> insert_cell(const Simplex_handle& simplex,
						  std::size_t cell_d,
						  bool is_boundary) {
    Simplex_cell_maps& simplex_cell_maps
      = (is_boundary? boundary_simplex_cell_maps_: interior_simplex_cell_maps_);
    // std::cout << "Insert simplex for "
    // 	      << (is_boundary? "\033[1;32mB" : "\033[1;33mI")
    // 	      << simplex << "\033[0m\n";

    Simplex_cell_map& simplex_cell_map = simplex_cell_maps[cell_d];
    auto curr_it = simplex_cell_map.lower_bound(simplex.smallest_coface());
    auto upper_bound = simplex_cell_map.upper_bound(simplex);
    for (; curr_it != upper_bound; ++curr_it) {
      // std::cout << " Checking if face: "
      // 		<< (is_boundary? "\033[1;32mB" : "\033[1;33mI")
      // 		<< curr_it->first << "\033[0m\n";
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
	  if (simplex == curr_it->first) {
	    curr_it++;
	    continue;
	  }
	  if (simplex.is_face_of(curr_it->first)) {
	    // std::cout << "  Post-deleting "
	    // 	      << (is_boundary? "\033[1;32mB" : "\033[1;33mI")
	    // 	      << curr_it->first << "\033[0m\n";
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
      // std::cout << " Checking if coface: "
      // 		<< (is_boundary? "\033[1;32mB" : "\033[1;33mI")
      // 		<< curr_it->first << "\033[0m\n";
      if (curr_it->first.is_face_of(simplex)) {
	// std::cout << "  Coface!\n";
	return std::make_pair(curr_it->second, Result_type::coface);
      }
    }
    hasse_cells_.push_back(new Hasse_cell(is_boundary, cell_d));
    Hasse_cell* new_cell = hasse_cells_.back();
    cell_simplex_map_.emplace(std::make_pair(new_cell, simplex));
    simplex_cell_map.emplace(std::make_pair(simplex, new_cell));
    // std::cout << " OK for insertion!\n";
    return std::make_pair(new_cell, Result_type::inserted);
  }

  void join_collapse_level(std::size_t cell_d,
			   bool is_boundary) {
    Simplex_cell_maps& simplex_cell_maps
      = (is_boundary? boundary_simplex_cell_maps_: interior_simplex_cell_maps_);
    Simplex_cell_map& simplex_cell_map = simplex_cell_maps[cell_d];
    auto curr_it = simplex_cell_map.begin();
    while (curr_it != simplex_cell_map.end()) {
      std::vector<Simplex_handle> faces;
      const Simplex_handle& simplex = curr_it->first;
      Hasse_cell* cell = curr_it->second;
      if (cell->get_boundary().size() == 1) {
	simplex_cell_map.erase(curr_it++);
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
	Hasse_cell* b_cell = bi_pair.first;
	const Simplex_handle& face = cell_simplex_map_.at(b_cell);
	bool b_is_boundary = b_cell->get_additional_information();
	if (face == join && is_boundary == b_is_boundary) {
	  join_is_face = true;
	  break;
	}
      }
      if (join_is_face)
	simplex_cell_map.erase(curr_it++);
      else if (simplex != join) {
	simplex_cell_map.erase(curr_it++);
	insert_cell(join, cell_d, is_boundary);
      }
      else
	curr_it++;
    }
  }

  void output_layer_before(const Simplex_cell_maps& simplex_cell_maps,
			   std::size_t cell_d,
			   std::size_t i,
			   bool is_before) const {
    std::cout << "\nThe layer " << cell_d << " "
	      << (is_before? "before": "after") << " the cleanup at simplex dimension " << i << "\n";
    for (auto& sc_pair: simplex_cell_maps[cell_d]) {
      Hasse_cell* cell = sc_pair.second;
      bool is_boundary = cell->get_additional_information();
      std::cout << "  " << (is_boundary? "\033[1;32mB" : "\033[1;33mI")
		<< sc_pair.first << "\033[0m: [ ";
      if (sc_pair.second->get_boundary().empty()) {
	std::cout << "]";
	continue;
      }
      auto b_it = cell->get_boundary().begin();
      Hasse_cell* b_cell = (b_it++)->first;
      const Simplex_handle& b_simplex = cell_simplex_map_.at(b_cell);
      bool b_is_boundary = b_cell->get_additional_information();
      std::cout << (b_is_boundary? "\033[1;32mB" : "\033[1;33mI")
		<< b_simplex << "\033[0m";
      for (; b_it != cell->get_boundary().end(); ++b_it) {
	Hasse_cell* b_cell = b_it->first;
	const Simplex_handle& b_simplex = cell_simplex_map_.at(b_cell);
	bool b_is_boundary = b_cell->get_additional_information();
	std::cout << ", " << (b_is_boundary? "\033[1;32mB" : "\033[1;33mI")
		  << b_simplex << "\033[0m";
      }
      std::cout << " ]";
      std::vector<Simplex_handle> faces;
      faces.reserve(cell->get_boundary().size());
      for (auto& bi_pair: cell->get_boundary()) {
	const Simplex_handle& face = cell_simplex_map_.at(bi_pair.first);
	faces.push_back(face);
      }
      Simplex_handle join = Gudhi::coxeter_triangulation::join(faces);
      std::cout << " join = " << (is_boundary? "\033[1;32mB" : "\033[1;33mI")
		<< join << "\033[0m\n";
    }
  }
  
  void expand_level(std::size_t cell_d) {
    std::size_t min_d = cell_d + intr_d_;
    std::size_t amb_d = interior_simplex_cell_maps_[0].begin()->first.vertex().size();
    for (auto& sc_pair: interior_simplex_cell_maps_[cell_d - 1])
      min_d = std::min(min_d, sc_pair.first.dimension());
    for (auto& sc_pair: boundary_simplex_cell_maps_[cell_d - 1])
      min_d = std::min(min_d, sc_pair.first.dimension());

    for (std::size_t i = min_d + 1; i <= amb_d; ++i) {
      if (cell_d != intr_d_) {
	for (auto& sc_pair: boundary_simplex_cell_maps_[cell_d - 1]) {
	  const Simplex_handle& simplex = sc_pair.first;
	  Hasse_cell* cell = sc_pair.second;
	  if (simplex.dimension() >= i)
	    continue;
	  for (Simplex_handle coface: simplex.coface_range(i)) {
	    Hasse_cell* new_cell;
	    Result_type success;
	    std::tie(new_cell, success) = insert_cell(coface, cell_d, true);
	    if (success == Result_type::self || success == Result_type::inserted)
	      new_cell->get_boundary().emplace_back(std::make_pair(cell, 1));
	  }
	}
	// output_layer_before(boundary_simplex_cell_maps_, cell_d, i, true);
	join_collapse_level(cell_d, true);
	// output_layer_before(boundary_simplex_cell_maps_, cell_d, i, false);
      }

      for (auto& sc_pair: interior_simplex_cell_maps_[cell_d - 1]) {
	const Simplex_handle& simplex = sc_pair.first;
        Hasse_cell* cell = sc_pair.second;
  	if (simplex.dimension() >= i)
  	  continue;
  	for (Simplex_handle coface: simplex.coface_range(i)) {
  	  Hasse_cell* new_cell;
  	  Result_type success;
  	  std::tie(new_cell, success) = insert_cell(coface, cell_d, false);
	  if (success == Result_type::inserted) {
	    Simplex_cell_map b_layer = boundary_simplex_cell_maps_[cell_d-1];
	    auto curr_it = b_layer.lower_bound(coface);
	    auto upper_bound = b_layer.upper_bound(coface.greatest_face());
	    for (; curr_it != upper_bound; ++curr_it) {
	      const Simplex_handle& b_simplex = curr_it->first;
	      Hasse_cell* b_cell = curr_it->second;
	      if (b_simplex.is_face_of(coface))
		new_cell->get_boundary().emplace_back(std::make_pair(b_cell, 1));
	    }
	  }
  	  if (success == Result_type::self || success == Result_type::inserted)
  	    new_cell->get_boundary().emplace_back(std::make_pair(cell, 1));
	}
      }
      // output_layer_before(interior_simplex_cell_maps_, cell_d, i, true);
      join_collapse_level(cell_d, false);
      // output_layer_before(interior_simplex_cell_maps_, cell_d, i, false);
    }
  }
  
  void construct_complex_(const Out_simplex_map_& out_simplex_map) {
    for (auto& os_pair: out_simplex_map) {
      const Simplex_handle& simplex = os_pair.first;
      const Eigen::VectorXd& point = os_pair.second;
      Hasse_cell* new_cell;
      Result_type success;	  
      std::tie(new_cell, success) = insert_cell(simplex, 0, false);
      if (success == Result_type::inserted)
	cell_point_map_.emplace(std::make_pair(new_cell, point));
      else if (success == Result_type::face)
	cell_point_map_.at(new_cell) = os_pair.second;	
    }    
    for (std::size_t cell_d = 1;
	 cell_d < interior_simplex_cell_maps_.size() &&
	   !interior_simplex_cell_maps_[cell_d - 1].empty();
	 ++cell_d) {
      expand_level(cell_d);
    }
  }
  
  void construct_complex_(const Out_simplex_map_& interior_simplex_map,
			  const Out_simplex_map_& boundary_simplex_map) {
    for (auto& os_pair: boundary_simplex_map) {
      const Simplex_handle& simplex = os_pair.first;
      const Eigen::VectorXd& point = os_pair.second;
      Hasse_cell* new_cell;
      Result_type success;
      std::tie(new_cell, success) = insert_cell(simplex, 0, true);
      if (success == Result_type::inserted)
	cell_point_map_.emplace(std::make_pair(new_cell, point));
      else if (success == Result_type::face)
	cell_point_map_.at(new_cell) = os_pair.second;
    }
    for (auto& os_pair: interior_simplex_map) {
      const Simplex_handle& simplex = os_pair.first;
      const Eigen::VectorXd& point = os_pair.second;
      auto curr_it = boundary_simplex_cell_maps_[0].lower_bound(simplex);
      auto upper_bound = boundary_simplex_cell_maps_[0].upper_bound(simplex.greatest_face());
      bool is_coface = false;
      for (; curr_it != upper_bound; ++curr_it)
	if (curr_it->first.is_face_of(simplex)) {
	  is_coface = true;
	  break;
	}
      if (is_coface)
	continue;
      Hasse_cell* new_cell;
      Result_type success;
      std::tie(new_cell, success) = insert_cell(simplex, 0, false);
      if (success == Result_type::inserted)
	cell_point_map_.emplace(std::make_pair(new_cell, point));
      else if (success == Result_type::face)
	cell_point_map_.at(new_cell) = os_pair.second;
    }
    std::cout << "Finished building the layer 0. Simplices:\n";
    // for (auto& sc_pair: interior_simplex_cell_maps_[0])
    //   std::cout << "\033[1;33mI" << sc_pair.first << "\033[0m\n";
    // for (auto& sc_pair: boundary_simplex_cell_maps_[0])
    //   std::cout << "\033[1;32mB" << sc_pair.first << "\033[0m\n";
    std::size_t i_size = interior_simplex_cell_maps_[0].size();
    std::size_t b_size = boundary_simplex_cell_maps_[0].size();
    std::cout << "I: " << i_size << "\n"
	      << "B: " << b_size << "\n"
	      << "T: " << i_size + b_size << "\n";

    for (std::size_t cell_d = 1;
	 cell_d < interior_simplex_cell_maps_.size() &&
	   !interior_simplex_cell_maps_[cell_d - 1].empty();
	 ++cell_d) {
      expand_level(cell_d);
      
      std::cout << "\nFinished building the layer " << cell_d << ". Simplices:\n";
      // for (auto& sc_pair: interior_simplex_cell_maps_[cell_d])
      // 	std::cout << "\033[1;33mI" << sc_pair.first << "\033[0m\n";
      // if (cell_d < intr_d_)
      // 	for (auto& sc_pair: boundary_simplex_cell_maps_[cell_d])
      // 	  std::cout << "\033[1;32mB" << sc_pair.first << "\033[0m\n";
      // std::cout << "\n";
      if (cell_d < intr_d_) {
	std::size_t i_size = interior_simplex_cell_maps_[cell_d].size();
	std::size_t b_size = boundary_simplex_cell_maps_[cell_d].size();
	std::cout << "I: " << i_size << "\n"
		  << "B: " << b_size << "\n"
		  << "T: " << i_size + b_size << "\n";
      }
      else {
	std::size_t i_size = interior_simplex_cell_maps_[cell_d].size();
	std::cout << "I: " << i_size << "\n"
		  << "T: " << i_size << "\n";
      }
    }
  }
  
public:  
  
  void construct_complex(const Out_simplex_map_& out_simplex_map) {
    interior_simplex_cell_maps_.resize(intr_d_ + 1);
    construct_complex_(out_simplex_map);
  }
  
  void construct_complex(const Out_simplex_map_& out_simplex_map,
			 std::size_t limit_dimension) {
    interior_simplex_cell_maps_.resize(limit_dimension + 1);
    construct_complex_(out_simplex_map);
  }

  void construct_complex(const Out_simplex_map_& interior_simplex_map,
			 const Out_simplex_map_& boundary_simplex_map) {
    interior_simplex_cell_maps_.resize(intr_d_ + 1);
    boundary_simplex_cell_maps_.resize(intr_d_);
    construct_complex_(interior_simplex_map, boundary_simplex_map);
  }

  void construct_complex(const Out_simplex_map_& interior_simplex_map,
			 const Out_simplex_map_& boundary_simplex_map,
			 std::size_t limit_dimension) {
    interior_simplex_cell_maps_.resize(limit_dimension + 1);
    boundary_simplex_cell_maps_.resize(limit_dimension);
    construct_complex_(interior_simplex_map, boundary_simplex_map);
  }

  std::size_t intrinsic_dimension() const {
    return intr_d_;
  }

  const Simplex_cell_maps& interior_simplex_cell_maps() const {
    return interior_simplex_cell_maps_;
  }

  const Simplex_cell_maps& boundary_simplex_cell_maps() const {
    return boundary_simplex_cell_maps_;
  }
  
  const Simplex_cell_map& simplex_cell_map(std::size_t cell_d) const {
    return interior_simplex_cell_maps_[cell_d];
  }

  const Simplex_cell_map& boundary_simplex_cell_map(std::size_t cell_d) const {
    return boundary_simplex_cell_maps_[cell_d];
  }

  const Cell_point_map& cell_point_map() const {
    return cell_point_map_;
  }

  Cell_complex(std::size_t intrinsic_dimension)
    : intr_d_(intrinsic_dimension) {}
  
private:
  std::size_t intr_d_;
  Simplex_cell_maps interior_simplex_cell_maps_, boundary_simplex_cell_maps_;
  Cell_simplex_map cell_simplex_map_;
  Cell_point_map cell_point_map_;
  // Cell_bool_map cell_boundary_map_;
  std::vector<Hasse_cell*> hasse_cells_;
};

} // namespace coxeter_triangulation 

} // namespace Gudhi

#endif
