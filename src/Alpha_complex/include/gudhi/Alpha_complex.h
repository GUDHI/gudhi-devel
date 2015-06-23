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
#include <limits> // NaN

namespace Gudhi {

namespace alphacomplex {

#define Kinit(f) =k.f()

/** \tableofcontents
 *  \defgroup alpha_complex Alpha complex
 * 
 * \author    Vincent Rouvreau
 * 
 * @{
 * 
 * \section definition Definition
 * 
 * Alpha_complex is a Simplex_tree constructed from each finite cell of a Delaunay Triangulation.
 * 
 * The filtration value of each simplex is computed from the alpha value of the simplex if it is Gabriel or
 * from the alpha value of the simplex coface that makes the simplex not Gabriel.
 * 
 * Please refer to \cite AlphaShapesDefinition for a more complete alpha complex definition.
 * 
 * Alpha complex are interesting because it looks like an \ref alpha-shape "Alpha shape" as described in
 * \cite AlphaShapesIntroduction (an alpha complex concept vulgarization).
 * 
 * \section example Example
 * 
 * This example loads points from an OFF file, builds the Delaunay triangulation from the points, and finally
 * initialize the alpha complex with it.
 * 
 * Then, it is asked to display information about the alpha complex.
 * 
 * \include Alpha_complex_from_off.cpp
 * 
 * When launching:
 * 
 * \code $> ./alphaoffreader ../../data/points/alphacomplexdoc.off
 * \endcode
 *
 * the program output is:
 * 
 * \include alphaoffreader_for_doc.txt
 * 
 * \section algorithm Algorithm
 * 
 * <b>Data structure</b>
 * 
 * In order to build the alpha complex, first, a Simplex tree is build from the cells of a Delaunay Triangulation.
 * (The filtration value is set to NaN, which stands for unknown value):
 * \image html "alpha_complex_doc.png" "Simplex tree structure construction example"
 *
 * <b>Filtration value computation algorithm</b>
 *
 * \f{algorithm}{ 
 * \caption{Filtration value computation algorithm}\label{alpha}
 * \begin{algorithmic}
 * \For{i : dimension $\rightarrow$ 1}
 *   \ForAll{$\sigma$ of dimension i}
 *     \If {filtration($\sigma$) is NaN}
 *       \State filtration($\sigma$) = $\alpha(\sigma)$
 *     \EndIf
 *     \ForAll{$\tau$ face of $\sigma$} \Comment{propagate alpha filtration value}
 *       \If {filtration($\tau$) is not NaN} 
 *         \State filtration($\tau$) = min (filtration($\tau$), filtration($\sigma$))
 *       \Else
 *         \If {$\tau$ is not Gabriel for $\sigma$} 
 *           \State filtration($\tau$) = filtration($\sigma$)
 *         \EndIf
 *       \EndIf
 *     \EndFor
 *   \EndFor
 * \EndFor
 * \end{algorithmic}
 * \f}
 * 
 * From the example above, it means the algorithm will look into each triangulation ([1,2,3], [2,3,4], [1,3,5], ...),
 * will compute the filtration value of the triangulation, and then will propagate the filtration value as described
 * here :
 * \image html "alpha_complex_doc_135.png" "Filtration value propagation example"
 * Then, the algorithm will look into each edge ([1,2], [2,3], [1,3], ...),
 * will compute the filtration value of the edge (in this case, propagation will have no effect).
 * 
 * Finally, the algorithm will look into each vertex ([1], [2], [3], [4], [5], [6] and [7]),
 * will set the filtration value (0 in case of a vertex - propagation will have no effect).
 * 
 * \section alpha-shape Alpha shape
 * 
 * In the example above, the alpha shape of \f$\alpha_{74} < \alpha < \alpha_{73}\f$ is the alpha complex where the 
 * \f$\alpha_{74} <\f$ filtration value \f$< \alpha_{73}\f$ as described in \cite AlphaShapesIntroduction
 * 
 * \image html "alpha_complex_doc_alpha_shape.png" "Alpha shape example"
 * \copyright GNU General Public License v3.                         
 * \verbatim  Contact: gudhi-users@lists.gforge.inria.fr \endverbatim
 */
/** @} */  // end defgroup alpha_complex

/**
 * \brief Alpha complex data structure.
 *
 * \details
 * The data structure can be constructed from a CGAL Delaunay triangulation (for more informations on CGAL Delaunay 
 * triangulation, please refer to the corresponding chapter in page http://doc.cgal.org/latest/Triangulation/) or from
 * an OFF file (cf. Delaunay_triangulation_off_reader).
 * 
 * Please refer to \ref alpha_complex for examples.
 *
 */
template<typename IndexingTag = linear_indexing_tag,
typename FiltrationValue = double,
typename SimplexKey = int,  // must be a signed integer type
typename VertexHandle = int  // must be a signed integer type, int convertible to it
>
class Alpha_complex : public Simplex_tree<> {
 private:
  // From Simplex_tree
  // Type required to insert into a simplex_tree (with or without subfaces).
  typedef std::vector<Vertex_handle> Vector_vertex;

  // Simplex_result is the type returned from simplex_tree insert function.
  typedef typename std::pair<Simplex_handle, bool> Simplex_result;

  // From CGAL
  // Kernel for the Delaunay_triangulation. Dimension can be set dynamically.
  typedef CGAL::Epick_d< CGAL::Dynamic_dimension_tag > Kernel;

  // Delaunay_triangulation type required to create an alpha-complex.
  typedef CGAL::Delaunay_triangulation<Kernel> Delaunay_triangulation;

  typedef typename Kernel::Compute_squared_radius_d Squared_Radius;
  typedef typename Kernel::Side_of_bounded_sphere_d Is_Gabriel;

  // Type required to compute squared radius, or side of bounded sphere on a vector of points.
  typedef std::vector<Kernel::Point_d> Vector_of_CGAL_points;

  // Vertex_iterator type from CGAL.
  typedef Delaunay_triangulation::Vertex_iterator CGAL_vertex_iterator;

  // Boost bimap type to switch from CGAL vertex iterator to simplex tree vertex handle and vice versa.
  typedef boost::bimap< CGAL_vertex_iterator, Vertex_handle > Bimap_vertex;

 private:
  /** \brief Boost bimap to switch from CGAL vertex iterator to simplex tree vertex handle and vice versa.*/
  Bimap_vertex cgal_simplextree;
  /** \brief Pointer on the CGAL Delaunay triangulation.*/
  Delaunay_triangulation* triangulation;

 public:
  /** \brief Alpha_complex constructor from an OFF file name.
   * Uses the Delaunay_triangulation_off_reader to construct the Delaunay triangulation required to initialize 
   * the Alpha_complex.
   *
   * @param[in] off_file_name OFF file [path and] name.
   */
  Alpha_complex(std::string& off_file_name)
      : triangulation(nullptr) {
    Gudhi::Delaunay_triangulation_off_reader<Delaunay_triangulation> off_reader(off_file_name);
    if (!off_reader.is_valid()) {
      std::cerr << "Alpha_complex - Unable to read file " << off_file_name << std::endl;
      exit(-1); // ----- >>
    }
    triangulation = off_reader.get_complex();
    init();
  }

  /** \brief Alpha_complex constructor from a Delaunay triangulation.
   *
   * @param[in] triangulation_ptr Pointer on a Delaunay triangulation.
   */
  Alpha_complex(Delaunay_triangulation* triangulation_ptr)
      : triangulation(triangulation_ptr) {
    init();
  }

  /** \brief Alpha_complex destructor from a Delaunay triangulation.
   *
   * @warning Deletes the Delaunay triangulation.
   */
  ~Alpha_complex() {
    delete triangulation;
  }

 private:
  /** \brief Initialize the Alpha_complex from the Delaunay triangulation.
   *
   * @warning Delaunay triangulation must be already constructed with at least one vertex and dimension must be more 
   * than 0.
   * 
   * Initialization can be launched once.
   */
  void init() {
    if (triangulation == nullptr) {
      std::cerr << "Alpha_complex init - Cannot init from a NULL triangulation" << std::endl;
      return; // ----- >>
    }
    if (triangulation->number_of_vertices() < 1) {
      std::cerr << "Alpha_complex init - Cannot init from a triangulation without vertices" << std::endl;
      return; // ----- >>
    }
    if (triangulation->maximal_dimension() < 1) {
      std::cerr << "Alpha_complex init - Cannot init from a zero-dimension triangulation" << std::endl;
      return; // ----- >>
    }
    if (num_vertices() > 0) {
      std::cerr << "Alpha_complex init - Cannot init twice" << std::endl;
      return; // ----- >>
    }
    
    set_dimension(triangulation->maximal_dimension());

    // --------------------------------------------------------------------------------------------
    // bimap to retrieve simplex tree vertex handles from CGAL vertex iterator and vice versa
    // Start to insert at handle = 0 - default integer value
    Vertex_handle vertex_handle = Vertex_handle();
    // Loop on triangulation vertices list
    for (CGAL_vertex_iterator vit = triangulation->vertices_begin(); vit != triangulation->vertices_end(); ++vit) {
      cgal_simplextree.insert(Bimap_vertex::value_type(vit, vertex_handle));
      vertex_handle++;
    }
    // --------------------------------------------------------------------------------------------

    // --------------------------------------------------------------------------------------------
    // Simplex_tree construction from loop on triangulation finite full cells list
    for (auto cit = triangulation->finite_full_cells_begin(); cit != triangulation->finite_full_cells_end(); ++cit) {
      Vector_vertex vertexVector;
#ifdef DEBUG_TRACES
      std::cout << "Simplex_tree insertion ";
#endif  // DEBUG_TRACES
      for (auto vit = cit->vertices_begin(); vit != cit->vertices_end(); ++vit) {
#ifdef DEBUG_TRACES
        std::cout << " " << cgal_simplextree.left.at(*vit);
#endif  // DEBUG_TRACES
        // Vector of vertex construction for simplex_tree structure
        vertexVector.push_back(cgal_simplextree.left.at(*vit));
      }
#ifdef DEBUG_TRACES
      std::cout << std::endl;
#endif  // DEBUG_TRACES
      // Insert each simplex and its subfaces in the simplex tree - filtration is NaN
      Simplex_result insert_result = insert_simplex_and_subfaces(vertexVector,
                                                                     std::numeric_limits<double>::quiet_NaN());
      if (!insert_result.second) {
        std::cerr << "Alpha_complex::init insert_simplex_and_subfaces failed" << std::endl;
      }
    }
    // --------------------------------------------------------------------------------------------

    Filtration_value filtration_max = 0.0;
    // --------------------------------------------------------------------------------------------
    // ### For i : d -> 0
    for (int decr_dim = dimension(); decr_dim >= 0; decr_dim--) {
      // ### Foreach Sigma of dim i
      for (auto f_simplex : skeleton_simplex_range(decr_dim)) {
        int f_simplex_dim = dimension(f_simplex);
        if (decr_dim == f_simplex_dim) {
          Vector_of_CGAL_points pointVector;
#ifdef DEBUG_TRACES
          std::cout << "Sigma of dim " << decr_dim << " is";
#endif  // DEBUG_TRACES
          for (auto vertex : simplex_vertex_range(f_simplex)) {
            pointVector.push_back((cgal_simplextree.right.at(vertex))->point());
#ifdef DEBUG_TRACES
            std::cout << " " << vertex;
#endif  // DEBUG_TRACES
          }
#ifdef DEBUG_TRACES
          std::cout << std::endl;
#endif  // DEBUG_TRACES
          // ### If filt(Sigma) is NaN : filt(Sigma) = alpha(Sigma)
          if (isnan(filtration(f_simplex))) {
            Filtration_value alpha_complex_filtration = 0.0;
            // No need to compute squared_radius on a single point - alpha is 0.0
            if (f_simplex_dim > 0) {
              // squared_radius function initialization
              Kernel k;
              Squared_Radius squared_radius Kinit(compute_squared_radius_d_object);

              alpha_complex_filtration = squared_radius(pointVector.begin(), pointVector.end());
            }
            assign_filtration(f_simplex, alpha_complex_filtration);
            filtration_max = fmax(filtration_max, alpha_complex_filtration);
#ifdef DEBUG_TRACES
            std::cout << "filt(Sigma) is NaN : filt(Sigma) =" << filtration(f_simplex) << std::endl;
#endif  // DEBUG_TRACES
          }
          propagate_alpha_filtration(f_simplex, decr_dim);
        }
      }
    }
    // --------------------------------------------------------------------------------------------

#ifdef DEBUG_TRACES
    std::cout << "filtration_max=" << filtration_max << std::endl;
#endif  // DEBUG_TRACES
    set_filtration(filtration_max);
  }

  template<typename Simplex_handle>
  void propagate_alpha_filtration(Simplex_handle f_simplex, int decr_dim) {
    // ### Foreach Tau face of Sigma
    for (auto f_boundary : boundary_simplex_range(f_simplex)) {
#ifdef DEBUG_TRACES
      std::cout << " | --------------------------------------------------" << std::endl;
      std::cout << " | Tau ";
      for (auto vertex : simplex_vertex_range(f_boundary)) {
        std::cout << vertex << " ";
      }
      std::cout << "is a face of Sigma" << std::endl;
      std::cout << " | isnan(filtration(Tau)=" << isnan(filtration(f_boundary)) << std::endl;
#endif  // DEBUG_TRACES
      // ### If filt(Tau) is not NaN
      if (!isnan(filtration(f_boundary))) {
        // ### filt(Tau) = fmin(filt(Tau), filt(Sigma))
        Filtration_value alpha_complex_filtration = fmin(filtration(f_boundary), filtration(f_simplex));
        assign_filtration(f_boundary, alpha_complex_filtration);
        // No need to check for filtration_max, alpha_complex_filtration is a min of an existing filtration value
#ifdef DEBUG_TRACES
        std::cout << " | filt(Tau) = fmin(filt(Tau), filt(Sigma)) = " << filtration(f_boundary) << std::endl;
#endif  // DEBUG_TRACES
        // ### Else
      } else {
        // No need to compute is_gabriel for dimension <= 2
        // i.e. : Sigma = (3,1) => Tau = 1
        if (decr_dim > 1) {
          // insert the Tau points in a vector for is_gabriel function
          Vector_of_CGAL_points pointVector;
          Vertex_handle vertexForGabriel = Vertex_handle();
          for (auto vertex : simplex_vertex_range(f_boundary)) {
            pointVector.push_back((cgal_simplextree.right.at(vertex))->point());
          }
          // Retrieve the Sigma point that is not part of Tau - parameter for is_gabriel function
          for (auto vertex : simplex_vertex_range(f_simplex)) {
            if (std::find(pointVector.begin(), pointVector.end(), (cgal_simplextree.right.at(vertex))->point())
                == pointVector.end()) {
              // vertex is not found in Tau
              vertexForGabriel = vertex;
              // No need to continue loop
              break;
            }
          }
          // is_gabriel function initialization
          Kernel k;
          Is_Gabriel is_gabriel Kinit(side_of_bounded_sphere_d_object);
#ifdef DEBUG_TRACES
          bool is_gab = is_gabriel(pointVector.begin(), pointVector.end(), (cgal_simplextree.right.at(vertexForGabriel))->point())
              != CGAL::ON_BOUNDED_SIDE;
          std::cout << " | Tau is_gabriel(Sigma)=" << is_gab << " - vertexForGabriel=" << vertexForGabriel << std::endl;
#endif  // DEBUG_TRACES
          // ### If Tau is not Gabriel of Sigma
          if ((is_gabriel(pointVector.begin(), pointVector.end(), (cgal_simplextree.right.at(vertexForGabriel))->point())
               == CGAL::ON_BOUNDED_SIDE)) {
            // ### filt(Tau) = filt(Sigma)
            Filtration_value alpha_complex_filtration = filtration(f_simplex);
            assign_filtration(f_boundary, alpha_complex_filtration);
            // No need to check for filtration_max, alpha_complex_filtration is an existing filtration value
#ifdef DEBUG_TRACES
            std::cout << " | filt(Tau) = filt(Sigma) = " << filtration(f_boundary) << std::endl;
#endif  // DEBUG_TRACES
          }
        }
      }
    }
  }
};

} // namespace alphacomplex

} // namespace Gudhi

#endif  // SRC_ALPHA_COMPLEX_INCLUDE_GUDHI_ALPHA_COMPLEX_H_
