/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2015  INRIA Saclay (France)
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SRC_ALPHA_SHAPES_INCLUDE_GUDHI_ALPHA_SHAPES_H_
#define SRC_ALPHA_SHAPES_INCLUDE_GUDHI_ALPHA_SHAPES_H_

// to construct a Delaunay_triangulation from a OFF file
#include <gudhi/Alpha_shapes/Delaunay_triangulation_off_io.h>

// to construct a simplex_tree from Delaunay_triangulation
#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/Simplex_tree.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>  // isnan, fmax

#include <boost/bimap.hpp>

#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/Epick_d.h>
#include <CGAL/algorithm.h>
#include <CGAL/assertions.h>
#include <CGAL/enum.h>

#include <iostream>
#include <iterator>
#include <vector>
#include <string>
#include <limits>
#include <map>

namespace Gudhi {

namespace alphashapes {

#define Kinit(f) =k.f()

/** \defgroup alpha_shapes Alpha shapes in dimension N
 *
 <DT>Implementations:</DT>
 Alpha shapes in dimension N are a subset of Delaunay Triangulation in dimension N.


 * \author    Vincent Rouvreau
 * \version   1.0
 * \date      2015
 * \copyright GNU General Public License v3.
 * @{
 */

/**
 * \brief Alpha shapes data structure.
 *
 * \details Every simplex \f$[v_0, \cdots ,v_d]\f$ admits a canonical orientation
 * induced by the order relation on vertices \f$ v_0 < \cdots < v_d \f$.
 *
 * Details may be found in \cite boissonnatmariasimplextreealgorithmica.
 *
 * \implements FilteredComplex
 *
 */
class Alpha_shapes {
 private:
  // From Simplex_tree
  /** \brief Type required to insert into a simplex_tree (with or without subfaces).*/
  typedef std::vector<Vertex_handle> typeVectorVertex;

  typedef typename Gudhi::Simplex_tree<>::Simplex_handle Simplex_handle;
  typedef typename std::pair<Simplex_handle, bool> Simplex_result;

  // From CGAL
  /** \brief Kernel for the Delaunay_triangulation.
   * Dimension can be set dynamically.
   */
  typedef CGAL::Epick_d< CGAL::Dynamic_dimension_tag > Kernel;
  /** \brief Delaunay_triangulation type required to create an alpha-shape.
   */
  typedef CGAL::Delaunay_triangulation<Kernel> Delaunay_triangulation;

  typedef typename Kernel::Compute_squared_radius_d Squared_Radius;
  typedef typename Kernel::Side_of_bounded_sphere_d Is_Gabriel;

  /** \brief Type required to insert into a simplex_tree (with or without subfaces).*/
  typedef std::vector<Kernel::Point_d> typeVectorPoint;

 private:
  /** \brief Upper bound on the simplex tree of the simplicial complex.*/
  Gudhi::Simplex_tree<> st_;

 public:

  Alpha_shapes(std::string& off_file_name) {
    // Construct a default Delaunay_triangulation (dim=0) - dim will be set in visitor reader init function
    Delaunay_triangulation dt(2);
    Gudhi::alphashapes::Delaunay_triangulation_off_reader<Delaunay_triangulation> off_reader(off_file_name, dt);
    if (!off_reader.is_valid()) {
      std::cerr << "Unable to read file " << off_file_name << std::endl;
      exit(-1); // ----- >>
    }
#ifdef DEBUG_TRACES
    std::cout << "number of vertices=" << dt.number_of_vertices() << std::endl;
    std::cout << "number of full cells=" << dt.number_of_full_cells() << std::endl;
    std::cout << "number of finite full cells=" << dt.number_of_finite_full_cells() << std::endl;
#endif  // DEBUG_TRACES
    init<Delaunay_triangulation>(dt);
  }

  template<typename T>
  Alpha_shapes(T& triangulation) {
    init<T>(triangulation);
  }

  ~Alpha_shapes() { }

 private:

  template<typename T>
  void init(T& triangulation) {
    st_.set_dimension(triangulation.maximal_dimension());

    // --------------------------------------------------------------------------------------------
    // bimap to retrieve vertex handles from points and vice versa
    typedef boost::bimap< Kernel::Point_d, Vertex_handle > bimap_points_vh;
    bimap_points_vh points_to_vh;
    // Start to insert at handle = 0 - default integer value
    Vertex_handle vertex_handle = Vertex_handle();
    // Loop on triangulation vertices list
    for (auto vit = triangulation.vertices_begin(); vit != triangulation.vertices_end(); ++vit) {
      points_to_vh.insert(bimap_points_vh::value_type(vit->point(), vertex_handle));
      vertex_handle++;
    }
    // --------------------------------------------------------------------------------------------

    // --------------------------------------------------------------------------------------------
    // Simplex_tree construction from loop on triangulation finite full cells list
    for (auto cit = triangulation.finite_full_cells_begin(); cit != triangulation.finite_full_cells_end(); ++cit) {
      typeVectorVertex vertexVector;
#ifdef DEBUG_TRACES
      std::cout << "Simplex_tree insertion ";
#endif  // DEBUG_TRACES
      for (auto vit = cit->vertices_begin(); vit != cit->vertices_end(); ++vit) {
#ifdef DEBUG_TRACES
        std::cout << " " << points_to_vh.left.at((*vit)->point());
#endif  // DEBUG_TRACES
        // Vector of vertex construction for simplex_tree structure
        vertexVector.push_back(points_to_vh.left.at((*vit)->point()));
      }
#ifdef DEBUG_TRACES
      std::cout << std::endl;
#endif  // DEBUG_TRACES
      // Insert each simplex and its subfaces in the simplex tree - filtration is NaN
      Simplex_result insert_result = st_.insert_simplex_and_subfaces(vertexVector,
                                                                     std::numeric_limits<double>::quiet_NaN());
      if (!insert_result.second) {
        std::cerr << "Alpha_shapes::init insert_simplex_and_subfaces failed" << std::endl;
      }
    }
    // --------------------------------------------------------------------------------------------

    Filtration_value filtration_max = 0.0;

    Kernel k;
    Squared_Radius squared_radius Kinit(compute_squared_radius_d_object);
    Is_Gabriel is_gabriel Kinit(side_of_bounded_sphere_d_object);
    // --------------------------------------------------------------------------------------------
    // ### For i : d -> 0
    for (int decr_dim = st_.dimension(); decr_dim >= 0; decr_dim--) {
      // ### Foreach Sigma of dim i
      for (auto f_simplex : st_.skeleton_simplex_range(decr_dim)) {
        int f_simplex_dim = st_.dimension(f_simplex);
        if (decr_dim == f_simplex_dim) {
          typeVectorPoint pointVector;
#ifdef DEBUG_TRACES
          std::cout << "Sigma of dim " << decr_dim << " is";
#endif  // DEBUG_TRACES
          for (auto vertex : st_.simplex_vertex_range(f_simplex)) {
            pointVector.push_back(points_to_vh.right.at(vertex));
#ifdef DEBUG_TRACES
            std::cout << " " << vertex;
#endif  // DEBUG_TRACES
          }
#ifdef DEBUG_TRACES
          std::cout << std::endl;
#endif  // DEBUG_TRACES
          // ### If filt(Sigma) is NaN : filt(Sigma) = alpha(Sigma)
          if (isnan(st_.filtration(f_simplex))) {
            Filtration_value alpha_shapes_filtration = 0.0;
            // No need to compute squared_radius on a single point - alpha is 0.0
            if (f_simplex_dim > 0) {
              alpha_shapes_filtration = squared_radius(pointVector.begin(), pointVector.end());
            }
            st_.assign_filtration(f_simplex, alpha_shapes_filtration);
            filtration_max = fmax(filtration_max, alpha_shapes_filtration);
#ifdef DEBUG_TRACES
            std::cout << "filt(Sigma) is NaN : filt(Sigma) =" << st_.filtration(f_simplex) << std::endl;
#endif  // DEBUG_TRACES
          }

          // ### Foreach Tau face of Sigma
          for (auto f_boundary : st_.boundary_simplex_range(f_simplex)) {
#ifdef DEBUG_TRACES
            std::cout << " | --------------------------------------------------" << std::endl;
            std::cout << " | Tau ";
            for (auto vertex : st_.simplex_vertex_range(f_boundary)) {
              std::cout << vertex << " ";
            }
            std::cout << "is a face of Sigma" << std::endl;
#endif  // DEBUG_TRACES
            // insert the Tau points in a vector for is_gabriel function
            typeVectorPoint pointVector;
            Vertex_handle vertexForGabriel = Vertex_handle();
            for (auto vertex : st_.simplex_vertex_range(f_boundary)) {
              pointVector.push_back(points_to_vh.right.at(vertex));
            }
            // Retrieve the Sigma point that is not part of Tau - parameter for is_gabriel function
            for (auto vertex : st_.simplex_vertex_range(f_simplex)) {
              if (std::find(pointVector.begin(), pointVector.end(), points_to_vh.right.at(vertex)) == pointVector.end()) {
                // vertex is not found in Tau
                vertexForGabriel = vertex;
                // No need to continue loop
                break;
              }
            }
#ifdef DEBUG_TRACES
            std::cout << " | isnan(filtration(Tau)=" << isnan(st_.filtration(f_boundary)) << std::endl;
            bool is_gab = is_gabriel(pointVector.begin(), pointVector.end(), points_to_vh.right.at(vertexForGabriel))
                != CGAL::ON_BOUNDED_SIDE;
            std::cout << " | Tau is_gabriel(Sigma)=" << is_gab << " - vertexForGabriel=" << vertexForGabriel << std::endl;
#endif  // DEBUG_TRACES
            // ### If filt(Tau) is not NaN
            // ### or Tau is not Gabriel of Sigma
            if (!isnan(st_.filtration(f_boundary)) ||
                (is_gabriel(pointVector.begin(), pointVector.end(), points_to_vh.right.at(vertexForGabriel)) == CGAL::ON_BOUNDED_SIDE)
                ) {
              // ### filt(Tau) = fmin(filt(Tau), filt(Sigma))
              Filtration_value alpha_shapes_filtration = fmin(st_.filtration(f_boundary), st_.filtration(f_simplex));
              st_.assign_filtration(f_boundary, alpha_shapes_filtration);
              filtration_max = fmax(filtration_max, alpha_shapes_filtration);
#ifdef DEBUG_TRACES
              std::cout << " | filt(Tau) = fmin(filt(Tau), filt(Sigma)) = " << st_.filtration(f_boundary) << std::endl;
#endif  // DEBUG_TRACES
            }
          }
        }
      }
    }
    // --------------------------------------------------------------------------------------------

#ifdef DEBUG_TRACES
    std::cout << "filtration_max=" << filtration_max << std::endl;
#endif  // DEBUG_TRACES
    st_.set_filtration(filtration_max);
  }

 public:

  /** \brief Returns the number of vertices in the complex. */
  size_t num_vertices() {
    return st_.num_vertices();
  }

  /** \brief Returns the number of simplices in the complex.
   *
   * Does not count the empty simplex. */
  const unsigned int& num_simplices() const {
    return st_.num_simplices();
  }

  /** \brief Returns an upper bound on the dimension of the simplicial complex. */
  int dimension() {
    return st_.dimension();
  }

  /** \brief Returns an upper bound of the filtration values of the simplices. */
  Filtration_value filtration() {
    return st_.filtration();
  }

  friend std::ostream& operator<<(std::ostream& os, const Alpha_shapes & alpha_shape) {
    // TODO: Program terminated with signal SIGABRT, Aborted - Maybe because of copy constructor
    Gudhi::Simplex_tree<> st = alpha_shape.st_;
    os << st << std::endl;
    return os;
  }
};

} // namespace alphashapes

} // namespace Gudhi

#endif  // SRC_ALPHA_SHAPES_INCLUDE_GUDHI_ALPHA_SHAPES_H_
