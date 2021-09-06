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
#include <set>
#include <map>
#include <algorithm>
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
  
  /** \internal \brief TODO
   *  \param[in] 
   *  \return  */
  template<typename Sphere>
  class CompareSpheresRadii
  {
    public:
      CGAL::NT_converter<FT, double> cast_to_double;
      bool operator()(const Sphere& firstSphere, const Sphere& secondSphere)
      {
        return cast_to_double(firstSphere.second) < cast_to_double(secondSphere.second);
      }
  };
  
  /** \internal \brief Čech complex blocker operator() - the oracle - assigns the filtration value from the simplex
   * radius and returns if the simplex expansion must be blocked.
   *  \param[in] sh The Simplex_handle.
   *  \return true if the simplex radius is greater than the Cech_complex max_radius*/
  bool operator()(Simplex_handle sh) {
    using Point_cloud =  std::vector<Point_d>;
    CGAL::NT_converter<FT, double> cast_to_double;
    Filtration_value radius = 0.;
//     std::string key_to_permute;
    std::vector<std::string> faces_keys;
    
    // for each face of simplex sh, test outsider point is indeed inside enclosing ball, if yes, take it and exit loop, otherwise, new sphere is circumsphere of all vertices
    // std::set <Filtration_value> enclosing_ball_radii;
    std::set <Sphere, CompareSpheresRadii<Sphere>> enclosing_ball_spheres;
    for (auto face : sc_ptr_->boundary_simplex_range(sh)) {
        
        // Find which vertex of sh is missing in face. We rely on the fact that simplex_vertex_range is sorted.
        auto longlist = sc_ptr_->simplex_vertex_range(sh);
        auto shortlist = sc_ptr_->simplex_vertex_range(face);

//         std::clog << "Hind debug: within FACE loop "<< std::endl;
        // TODO to remove
//         for (auto i = std::begin(longlist); i != std::end(longlist);++i)
//             std::clog << "Hind debug: longlist: " << cc_ptr_->get_point(*i) << std::endl;
//         for (auto i = std::begin(shortlist); i != std::end(shortlist);++i)
//             std::clog << "Hind debug: shortlist: " << cc_ptr_->get_point(*i) << std::endl;

        auto longiter = std::begin(longlist);
        auto shortiter = std::begin(shortlist);
        auto enditer = std::end(shortlist);
        while(shortiter != enditer && *longiter == *shortiter) { ++longiter; ++shortiter; }
        auto extra = *longiter; // Vertex_handle

//         std::clog << "Hind debug: extra vertex: " << cc_ptr_->get_point(extra) << std::endl;
        
        Point_cloud face_points;
        std::string key, key_extra;
        for (auto vertex : sc_ptr_->simplex_vertex_range(face)) {
            face_points.push_back(cc_ptr_->get_point(vertex));
            key.append(std::to_string(vertex));
    #ifdef DEBUG_TRACES
        std::clog << "#(" << vertex << ")#";
    #endif  // DEBUG_TRACES
        }
        key_extra = key;
        key_extra.append(std::to_string(extra));
        faces_keys.push_back(key_extra);
//         key_to_permute = key_extra;
//         std::clog << "END OF VERTICES " << std::endl;
//         std::clog << "KEY is: " << key << std::endl;
//         std::clog << "KEY extra is: " << key_extra << std::endl;
        Sphere sph;
        auto it = cache_.find(key);
        if(it != cache_.end())
            sph = it->second;
        else {
            sph = get_sphere(face_points.cbegin(), face_points.cend());
        }
        if (kernel_.squared_distance_d_object()(sph.first, cc_ptr_->get_point(extra)) <= sph.second) {
            radius = std::sqrt(cast_to_double(sph.second));
            #ifdef DEBUG_TRACES
                std::clog << "circumcenter: " << sph.first << ", radius: " <<  radius << std::endl;
            #endif  // DEBUG_TRACES
//             std::clog << "distance FYI: " << kernel_.squared_distance_d_object()(sph.first, cc_ptr_->get_point(extra)) << " < " << cast_to_double(sph.second) << std::endl;
            // enclosing_ball_radii.insert(radius);
            enclosing_ball_spheres.insert(sph);
            cache_[key_extra] = sph;
        }
//         else {// TODO to remove
//             std::clog << "vertex not included BECAUSE DISTANCE: "<< kernel_.squared_distance_d_object()(sph.first, cc_ptr_->get_point(extra)) << " AND RAD SPHERE: " << sph.second << std::endl;
//         }
    }
    // Get the minimal radius of all faces enclosing balls if exists
    if (!enclosing_ball_spheres.empty()) {
        // radius = *enclosing_ball_radii.begin();
        Sphere sph_min = *enclosing_ball_spheres.begin();
        radius = std::sqrt(cast_to_double(sph_min.second));
        // std::clog << "CHECK that radius of min sphere is min radius: " << std::sqrt(cast_to_double(sph_min.second)) << "; and RADIUS min: " << radius << std::endl;
        // Set all key_to_permute permutations to min sphere in cache
//         do
//         {   
//             if (cache_.find(key_to_permute) != cache_.end()) {
//                 if (cast_to_double(cache_[key_to_permute].second) >  cast_to_double(sph_min.second))
//                     cache_[key_to_permute] = sph_min;
//             }
//             else {
//                 cache_[key_to_permute] = sph_min;
//             }
//         } while(std::next_permutation(key_to_permute.begin(), key_to_permute.end()));
        for (auto k : faces_keys) {
            cache_[k] = sph_min;
        }
    }
//     std::clog << "END OF FACES ; radius = " << radius << std::endl;
    
    if (radius == 0.) { // Spheres of each face don't contain the whole simplex
        Point_cloud points;
        for (auto vertex : sc_ptr_->simplex_vertex_range(sh)) {
            points.push_back(cc_ptr_->get_point(vertex));
        }
        Sphere sph = get_sphere(points.cbegin(), points.cend());
        radius = std::sqrt(cast_to_double(sph.second));
//         std::clog << "GLOBAL SPHERE radius = " << radius << std::endl;
        // Set all key_to_permute permutations to sphere in cache
//         do
//         {
// //             if (cache_.find(key_to_permute) != cache_.end()) {
// //                 if (cast_to_double(cache_[key_to_permute].second) >  cast_to_double(sph.second))
// //                     cache_[key_to_permute] = sph;
// //             }
// //             else {
// //                 cache_[key_to_permute] = sph;
// //             }
//         } while(std::next_permutation(key_to_permute.begin(), key_to_permute.end()));
//         for (auto k : faces_keys) {
//             cache_[k] = sph;
//         }
    }


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
  std::map<std::string, Sphere> cache_;
};

}  // namespace cech_complex

}  // namespace Gudhi

#endif  // CECH_COMPLEX_BLOCKER_H_
