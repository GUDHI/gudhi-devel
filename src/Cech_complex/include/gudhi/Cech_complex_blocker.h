/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2018 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CECH_COMPLEX_BLOCKER_H_
#define CECH_COMPLEX_BLOCKER_H_

// TODO to remove
#include <gudhi/distance_functions.h>  // for Gudhi::Minimal_enclosing_ball_radius

#include <CGAL/number_utils.h> // for CGAL::to_double
#include <CGAL/NT_converter.h> //
// #include <CGAL/Delaunay_triangulation.h>
#include <CGAL/Epeck_d.h>  // For EXACT or SAFE version
#include <CGAL/Spatial_sort_traits_adapter_d.h>

#include <iostream>
#include <vector>
#include <cmath>  // for std::sqrt

namespace Gudhi {

namespace cech_complex {

/** \internal
 * \class Cech_blocker
 * \brief Čech complex blocker.
 *
 * \ingroup cech_complex
 *
 * \details
 * Čech blocker is an oracle constructed from a Cech_complex and a simplicial complex.
 *
 * \tparam SimplicialComplexForProximityGraph furnishes `Simplex_handle` and `Filtration_value` type definition,
 * `simplex_vertex_range(Simplex_handle sh)`and `assign_filtration(Simplex_handle sh, Filtration_value filt)` methods.
 *
 * \tparam Chech_complex is required by the blocker.
 */
template <typename SimplicialComplexForCech, typename Cech_complex, typename Kernel>
class Cech_blocker {
 private:
//   using Point_cloud = typename Cech_complex::Point_cloud;

  using Simplex_handle = typename SimplicialComplexForCech::Simplex_handle;
  using Filtration_value = typename SimplicialComplexForCech::Filtration_value;

 public:

  using Point_d = typename Kernel::Point_d;
  // Numeric type of coordinates in the kernel
  using FT = typename Kernel::FT;
  // Sphere is a pair of point and squared radius.
  using Sphere = typename std::pair<Point_d, FT>;
  
  
  // Add an int in TDS to save point index in the structure
//   using TDS = CGAL::Triangulation_data_structure<typename Kernel::Dimension,
//                                                  CGAL::Triangulation_vertex<Kernel, std::ptrdiff_t>,
//                                                  CGAL::Triangulation_full_cell<Kernel> >;
// 
//   /** \brief A (Weighted or not) Delaunay triangulation of a set of points in \f$ \mathbb{R}^D\f$.*/
//   using Triangulation = CGAL::Delaunay_triangulation<Kernel, TDS>;
//   // Vertex_iterator type from CGAL.
//   using CGAL_vertex_iterator = typename Triangulation::Vertex_iterator;

  // Structure to switch from simplex tree vertex handle to CGAL vertex iterator.
  //using Vector_vertex_iterator = std::vector< CGAL_vertex_iterator >;

  
  
    /** \brief get_point_ returns the point corresponding to the vertex given as parameter.
   * Only for internal use for faster access.
   *
   * @param[in] vertex Vertex handle of the point to retrieve.
   * @return The point found.
   */
/*  const Point_d& get_point_(std::size_t vertex) const {
    return vertex_handle_to_iterator_[vertex]->point();
  } */

   /** \internal \brief TODO
   *  \param[in] 
   *  \return  */
  template<class PointIterator>
  FT get_squared_radius(PointIterator begin, PointIterator end) const {
    return kernel_.compute_squared_radius_d_object()(begin, end);
  }
  
  /** \internal \brief TODO
   *  \param[in] 
   *  \return  */
  template<class PointIterator>
  Sphere get_sphere(PointIterator begin, PointIterator end) const {
    Point_d c = kernel_.construct_circumcenter_d_object()(begin, end);
    FT r = kernel_.squared_distance_d_object()(c, *begin);
    return std::make_pair(std::move(c), std::move(r));
  }
  
  
  /** \internal \brief Čech complex blocker operator() - the oracle - assigns the filtration value from the simplex
   * radius and returns if the simplex expansion must be blocked.
   *  \param[in] sh The Simplex_handle.
   *  \return true if the simplex radius is greater than the Cech_complex max_radius*/
  bool operator()(Simplex_handle sh) {
    using Point_cloud =  std::vector<Point_d>;
    Point_cloud points;
    
    // for each face of simplex sh, test outsider point is indeed inside enclosing ball, if yes, take it and exit loop, otherwise, new sphere is circumsphere of all vertices
    for (auto face : sc_ptr_->simplex_vertex_range(sh)) {
        /////////////////////////////////////////////////////////////////////
        
        
        
        
        ///////////////////////////////////
    }
    
    for (auto vertex : sc_ptr_->simplex_vertex_range(sh)) {
        points.push_back(cc_ptr_->get_point(vertex));
//        points.push_back(get_point_(vertex));
#ifdef DEBUG_TRACES
      std::clog << "#(" << vertex << ")#";
#endif  // DEBUG_TRACES
    }
    // TODO to remove
    //Filtration_value radius = Gudhi::Minimal_enclosing_ball_radius()(points);
    // Hind: Here change the algo of the enclosing Minimal_enclosing_ball_radius
    auto point_to_be_inserted = points.back();
    Sphere sph = get_sphere(points.cbegin(), points.cend()-1);

//     Sphere sph = get_sphere(points.cbegin(), points.cend()-1);
    CGAL::NT_converter<FT, double> cast_to_double;
//     CGAL::NT_converter<typename Cech_complex::Point, Point_d> cast_point_d_to_double;
    
    std::clog << "circumcenter: " << sph.first << ", radius: " <<  std::sqrt(cast_to_double(sph.second))<< std::endl;
    // TODO to remove
    // Filtration_value test = std::sqrt(CGAL::to_double(sph.second));

    
    // Check that the point to be inserted is already included in the sphere of the simplex containing the preceding points
    // TODO instead of Euclidean_distance ; use kernel_.squared_distance_d_object()(c, *begin);
    // Add a loop on the three faces to check sphere before computing the circumsphere
    // Add the computed sphere as cache; a vector of spheres depending on the number of faces ?
    //
//     if (Gudhi::Euclidean_distance()(cast_point_d_to_double(sph.first), point_to_be_inserted) > std::sqrt(cast_to_double(sph.second)))
//     FT r = kernel_.squared_distance_d_object()(sph.first, sph.first); //*(points.cend()-1));
    if (kernel_.squared_distance_d_object()(sph.first, point_to_be_inserted) > sph.second)
        sph = get_sphere(points.cbegin(), points.cend());

    Filtration_value radius = std::sqrt(cast_to_double(sph.second));


#ifdef DEBUG_TRACES
    if (radius > cc_ptr_->max_radius()) std::clog << "radius > max_radius => expansion is blocked\n";
#endif  // DEBUG_TRACES
    sc_ptr_->assign_filtration(sh, radius);
    return (radius > cc_ptr_->max_radius());
  }

  /** \internal \brief Čech complex blocker constructor. */
  Cech_blocker(SimplicialComplexForCech* sc_ptr, Cech_complex* cc_ptr) : sc_ptr_(sc_ptr), cc_ptr_(cc_ptr) {}

 private:
  SimplicialComplexForCech* sc_ptr_;
  Cech_complex* cc_ptr_;
  Kernel kernel_;
  //Vector_vertex_iterator vertex_handle_to_iterator_;

};

}  // namespace cech_complex

}  // namespace Gudhi

#endif  // CECH_COMPLEX_BLOCKER_H_
