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

namespace Gudhi {

namespace coxeter_triangulation {

/**
 *  \ingroup coxeter_triangulation
 */

/** \class Cell_complex 
 *  \brief A class that constructs the cell complex from the output provided by the class 
 *   Gudhi::Manifold_tracing<Triangulation_>.
 *
 *  \tparam Out_simplex_map_ The type of a map from a simplex type that is a
 *   model of SimplexInCoxeterTriangulation to Eigen::VectorXd.
 *
 *  \ingroup coxeter_triangulation
 */
template <class Out_simplex_map_>
class Cell_complex {
public:
  
  /** \brief Type of a simplex in the ambient triangulation.
   *  Is a model of the concept SimplexInCoxeterTriangulation.
   */
  typedef typename Out_simplex_map_::key_type Simplex_handle;
  /** \brief Type of a cell in the cell complex. 
   *  Always is Gudhi::Hasse_cell from the Hasse diagram module.
   *  The additional information is the boolean that is true if and only if the cell lies
   *  on the boundary.
   */
  typedef Gudhi::Hasse_diagram::Hasse_diagram_cell<int, double, bool> Hasse_cell;
  /** \brief Type of a map from permutahedral representations of simplices in the
   *  ambient triangulation to the corresponding cells in the cell complex of some
   *  specific dimension.
   */
  typedef std::map<Simplex_handle,
		   Hasse_cell*,
		   Simplex_comparator<Simplex_handle> > Simplex_cell_map;
  /** \brief Type of a vector of maps from permutahedral representations of simplices in the
   *  ambient triangulation to the corresponding cells in the cell complex of various dimensions.
   */
  typedef std::vector<Simplex_cell_map> Simplex_cell_maps;
  
  /** \brief Type of a map from cells in the cell complex to the permutahedral representations
   *  of the corresponding simplices in the ambient triangulation.
   */
  typedef std::map<Hasse_cell*, Simplex_handle> Cell_simplex_map;

  /** \brief Type of a map from vertex cells in the cell complex to the permutahedral representations
   *  of their Cartesian coordinates.
   */
  typedef std::map<Hasse_cell*, Eigen::VectorXd> Cell_point_map;

private:
    
  Hasse_cell* insert_cell(const Simplex_handle& simplex,
			  std::size_t cell_d,
			  bool is_boundary) {
    Simplex_cell_maps& simplex_cell_maps
      = (is_boundary? boundary_simplex_cell_maps_ : interior_simplex_cell_maps_); 
#ifdef GUDHI_COX_OUTPUT_TO_HTML
    CC_detail_list& cc_detail_list =
      (is_boundary? cc_boundary_detail_lists[cell_d] : cc_interior_detail_lists[cell_d]);
    cc_detail_list.emplace_back(CC_detail_info(simplex));
#endif
    Simplex_cell_map& simplex_cell_map = simplex_cell_maps[cell_d];
    auto map_it = simplex_cell_map.find(simplex);
    if (map_it == simplex_cell_map.end()) {
      hasse_cells_.push_back(new Hasse_cell(is_boundary, cell_d));
      Hasse_cell* new_cell = hasse_cells_.back();    
      simplex_cell_map.emplace(std::make_pair(simplex, new_cell));
      cell_simplex_map_.emplace(std::make_pair(new_cell, simplex));
#ifdef GUDHI_COX_OUTPUT_TO_HTML
      cc_detail_list.back().status_ = CC_detail_info::Result_type::inserted;
#endif
      return new_cell;
    }
#ifdef GUDHI_COX_OUTPUT_TO_HTML
    CC_detail_info& cc_info = cc_detail_list.back();
    cc_info.trigger_ = to_string(map_it->first);
    cc_info.status_ = CC_detail_info::Result_type::self;
#endif
    return map_it->second;
  }

  void expand_level(std::size_t cell_d) {
    bool is_manifold_with_boundary = boundary_simplex_cell_maps_.size() > 0;
    for (auto& sc_pair: interior_simplex_cell_maps_[cell_d - 1]) {
      const Simplex_handle& simplex = sc_pair.first;
      Hasse_cell* cell = sc_pair.second;
      for (Simplex_handle coface: simplex.coface_range(cod_d_ + cell_d)) {
	Hasse_cell* new_cell = insert_cell(coface, cell_d, false);
	new_cell->get_boundary().emplace_back(std::make_pair(cell, 1));
      }
    }
    
    if (is_manifold_with_boundary) {
      for (auto& sc_pair: boundary_simplex_cell_maps_[cell_d - 1]) {
	const Simplex_handle& simplex = sc_pair.first;
	Hasse_cell* cell = sc_pair.second;
	if (cell_d != intr_d_)
	  for (Simplex_handle coface: simplex.coface_range(cod_d_ + cell_d + 1)) {
	    Hasse_cell* new_cell = insert_cell(coface, cell_d, true);
	    new_cell->get_boundary().emplace_back(std::make_pair(cell, 1));
	  }
	auto map_it = interior_simplex_cell_maps_[cell_d].find(simplex);
	if (map_it == interior_simplex_cell_maps_[cell_d].end())
	  std::cerr << "Cell_complex::expand_level error: A boundary cell does not have an interior counterpart.\n";
	else {
	  Hasse_cell* i_cell = map_it->second;
	  i_cell->get_boundary().emplace_back(std::make_pair(cell, 1));
	}
      }
    }
  }
  
  void construct_complex_(const Out_simplex_map_& out_simplex_map) {
#ifdef GUDHI_COX_OUTPUT_TO_HTML
    cc_interior_summary_lists.resize(interior_simplex_cell_maps_.size());
    cc_interior_prejoin_lists.resize(interior_simplex_cell_maps_.size());
    cc_interior_detail_lists.resize(interior_simplex_cell_maps_.size());
#endif    
    for (auto& os_pair: out_simplex_map) {
      const Simplex_handle& simplex = os_pair.first;
      const Eigen::VectorXd& point = os_pair.second;
      Hasse_cell* new_cell = insert_cell(simplex, 0, false);
      cell_point_map_.emplace(std::make_pair(new_cell, point));
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
#ifdef GUDHI_COX_OUTPUT_TO_HTML
    cc_interior_summary_lists.resize(interior_simplex_cell_maps_.size());
    cc_interior_prejoin_lists.resize(interior_simplex_cell_maps_.size());
    cc_interior_detail_lists.resize(interior_simplex_cell_maps_.size());
    cc_boundary_summary_lists.resize(boundary_simplex_cell_maps_.size());
    cc_boundary_prejoin_lists.resize(boundary_simplex_cell_maps_.size());
    cc_boundary_detail_lists.resize(boundary_simplex_cell_maps_.size());
#endif    
    for (auto& os_pair: boundary_simplex_map) {
      const Simplex_handle& simplex = os_pair.first;
      const Eigen::VectorXd& point = os_pair.second;
      Hasse_cell* new_cell = insert_cell(simplex, 0, true);
      cell_point_map_.emplace(std::make_pair(new_cell, point));
    }
    for (auto& os_pair: interior_simplex_map) {
      const Simplex_handle& simplex = os_pair.first;
      const Eigen::VectorXd& point = os_pair.second;
      Hasse_cell* new_cell = insert_cell(simplex, 0, false);
      cell_point_map_.emplace(std::make_pair(new_cell, point));
    }
#ifdef GUDHI_COX_OUTPUT_TO_HTML
    for (const auto& sc_pair: interior_simplex_cell_maps_[0])
      cc_interior_summary_lists[0].push_back(CC_summary_info(sc_pair));
    for (const auto& sc_pair: boundary_simplex_cell_maps_[0])
      cc_boundary_summary_lists[0].push_back(CC_summary_info(sc_pair));
#endif        

    for (std::size_t cell_d = 1;
	 cell_d < interior_simplex_cell_maps_.size() &&
	   !interior_simplex_cell_maps_[cell_d - 1].empty();
	 ++cell_d) {
      expand_level(cell_d);
      
#ifdef GUDHI_COX_OUTPUT_TO_HTML
      for (const auto& sc_pair: interior_simplex_cell_maps_[cell_d])
	cc_interior_summary_lists[cell_d].push_back(CC_summary_info(sc_pair));
      if (cell_d < boundary_simplex_cell_maps_.size())
	for (const auto& sc_pair: boundary_simplex_cell_maps_[cell_d])
	  cc_boundary_summary_lists[cell_d].push_back(CC_summary_info(sc_pair));
#endif
    }
  }
  
public:  
  
  /**
   * \brief Constructs the the cell complex that approximates an \f$m\f$-dimensional manifold
   *  without boundary embedded in the \f$ d \f$-dimensional Euclidean space
   *  from the output of the class Gudhi::Manifold_tracing.
   *
   * \param[in] out_simplex_map A map from simplices of dimension \f$(d-m)\f$ 
   * in the ambient triangulation that intersect the relative interior of the manifold 
   * to the intersection points.
   */
  void construct_complex(const Out_simplex_map_& out_simplex_map) {
    interior_simplex_cell_maps_.resize(intr_d_ + 1);
    if (!out_simplex_map.empty())
      cod_d_ = out_simplex_map.begin()->first.dimension();
    construct_complex_(out_simplex_map);
  }
  
  /**
   * \brief Constructs the skeleton of the cell complex that approximates
   *  an \f$m\f$-dimensional manifold without boundary embedded 
   *  in the \f$d\f$-dimensional Euclidean space
   *  up to a limit dimension from the output of the class Gudhi::Manifold_tracing.
   *
   * \param[in] out_simplex_map A map from simplices of dimension \f$(d-m)\f$ 
   * in the ambient triangulation that intersect the relative interior of the manifold 
   * to the intersection points.
   * \param[in] limit_dimension The dimension of the constructed skeleton.
   */
  void construct_complex(const Out_simplex_map_& out_simplex_map,
			 std::size_t limit_dimension) {
    interior_simplex_cell_maps_.resize(limit_dimension + 1);
    if (!out_simplex_map.empty())
      cod_d_ = out_simplex_map.begin()->first.dimension();
    construct_complex_(out_simplex_map);
  }

  /**
   * \brief Constructs the the cell complex that approximates an \f$m\f$-dimensional manifold
   *  with boundary embedded in the \f$ d \f$-dimensional Euclidean space
   *  from the output of the class Gudhi::Manifold_tracing.
   *
   * \param[in] interior_simplex_map A map from simplices of dimension \f$(d-m)\f$ 
   * in the ambient triangulation that intersect the relative interior of the manifold 
   * to the intersection points.
   * \param[in] boundary_simplex_map A map from simplices of dimension \f$(d-m+1)\f$ 
   * in the ambient triangulation that intersect the boundary of the manifold 
   * to the intersection points.
   */
  void construct_complex(const Out_simplex_map_& interior_simplex_map,
			 const Out_simplex_map_& boundary_simplex_map) {
    interior_simplex_cell_maps_.resize(intr_d_ + 1);
    boundary_simplex_cell_maps_.resize(intr_d_);
    if (!interior_simplex_map.empty())
      cod_d_ = interior_simplex_map.begin()->first.dimension();
    construct_complex_(interior_simplex_map, boundary_simplex_map);
  }

  /**
   * \brief Constructs the skeleton of the cell complex that approximates
   *  an \f$m\f$-dimensional manifold with boundary embedded 
   *  in the \f$d\f$-dimensional Euclidean space
   *  up to a limit dimension from the output of the class Gudhi::Manifold_tracing.
   *
   * \param[in] interior_simplex_map A map from simplices of dimension \f$(d-m)\f$ 
   * in the ambient triangulation that intersect the relative interior of the manifold 
   * to the intersection points.
   * \param[in] boundary_simplex_map A map from simplices of dimension \f$(d-m+1)\f$ 
   * in the ambient triangulation that intersect the boundary of the manifold 
   * to the intersection points.
   * \param[in] limit_dimension The dimension of the constructed skeleton.
   */
  void construct_complex(const Out_simplex_map_& interior_simplex_map,
			 const Out_simplex_map_& boundary_simplex_map,
			 std::size_t limit_dimension) {
    interior_simplex_cell_maps_.resize(limit_dimension + 1);
    boundary_simplex_cell_maps_.resize(limit_dimension);
    if (!interior_simplex_map.empty())
      cod_d_ = interior_simplex_map.begin()->first.dimension();
    construct_complex_(interior_simplex_map, boundary_simplex_map);
  }

  /**
   * \brief Returns the dimension of the cell complex.
   */
  std::size_t intrinsic_dimension() const {
    return intr_d_;
  }

  /**
   * \brief Returns a vector of maps from the cells of various dimensions in the interior
   *  of the cell complex of type Gudhi::Hasse_cell to the permutahedral representations 
   *  of the corresponding simplices in the ambient triangulation.
   */
  const Simplex_cell_maps& interior_simplex_cell_maps() const {
    return interior_simplex_cell_maps_;
  }

  /**
   * \brief Returns a vector of maps from the cells of various dimensions on the boundary
   *  of the cell complex of type Gudhi::Hasse_cell to the permutahedral representations 
   *  of the corresponding simplices in the ambient triangulation.
   */
  const Simplex_cell_maps& boundary_simplex_cell_maps() const {
    return boundary_simplex_cell_maps_;
  }
  
  /**
   * \brief Returns a map from the cells of a given dimension in the interior 
   *  of the cell complex of type Gudhi::Hasse_cell to the permutahedral representations 
   *  of the corresponding simplices in the ambient triangulation.
   *
   * \param[in] cell_d The dimension of the cells.
   */
  const Simplex_cell_map& interior_simplex_cell_map(std::size_t cell_d) const {
    return interior_simplex_cell_maps_[cell_d];
  }

  /**
   * \brief Returns a map from the cells of a given dimension on the boundary 
   *  of the cell complex of type Gudhi::Hasse_cell to the permutahedral representations 
   *  of the corresponding simplices in the ambient triangulation.
   *
   * \param[in] cell_d The dimension of the cells.
   */
  const Simplex_cell_map& boundary_simplex_cell_map(std::size_t cell_d) const {
    return boundary_simplex_cell_maps_[cell_d];
  }

  /**
   * \brief Returns a map from the cells in the cell complex of type Gudhi::Hasse_cell
   *  to the permutahedral representations of the corresponding simplices in the 
   *  ambient triangulation.
   */
  const Cell_simplex_map& cell_simplex_map() const {
    return cell_simplex_map_;
  }

  /**
   * \brief Returns a map from the vertex cells in the cell complex of type Gudhi::Hasse_cell
   *  to their Cartesian coordinates.
   */
  const Cell_point_map& cell_point_map() const {
    return cell_point_map_;
  }

  /**
   * \brief Conxtructor for the class Cell_complex.
   *
   * \param[in] intrinsic_dimension The dimension of the cell complex.
   */
  Cell_complex(std::size_t intrinsic_dimension)
    : intr_d_(intrinsic_dimension) {}

  ~Cell_complex() {
    for (Hasse_cell* hs_ptr: hasse_cells_)
      delete hs_ptr;
  }
  
private:
  std::size_t intr_d_, cod_d_;
  Simplex_cell_maps interior_simplex_cell_maps_, boundary_simplex_cell_maps_;
  Cell_simplex_map cell_simplex_map_;
  Cell_point_map cell_point_map_;
  std::vector<Hasse_cell*> hasse_cells_;
};

} // namespace coxeter_triangulation 

} // namespace Gudhi

#endif
