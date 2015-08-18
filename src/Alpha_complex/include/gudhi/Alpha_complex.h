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

#ifndef ALPHA_COMPLEX_H_
#define ALPHA_COMPLEX_H_

// to construct a simplex_tree from Delaunay_triangulation
#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/Simplex_tree.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>  // isnan, fmax

#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/Epick_d.h>
#include <CGAL/algorithm.h>
#include <CGAL/assertions.h>
#include <CGAL/enum.h>

#include <iostream>
#include <vector>
#include <string>
#include <limits>  // NaN
#include <map>
#include <utility>  // std::pair

namespace Gudhi {

namespace alphacomplex {

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
template<class Kernel>
class Alpha_complex : public Simplex_tree<> {
 private:
  // From Simplex_tree
  // Type required to insert into a simplex_tree (with or without subfaces).
  typedef std::vector<Vertex_handle> Vector_vertex;

  // Simplex_result is the type returned from simplex_tree insert function.
  typedef typename std::pair<Simplex_handle, bool> Simplex_result;

  // Delaunay_triangulation type required to create an alpha-complex.
  typedef typename CGAL::Delaunay_triangulation<Kernel> Delaunay_triangulation;

  typedef typename Kernel::Compute_squared_radius_d Squared_Radius;
  typedef typename Kernel::Side_of_bounded_sphere_d Is_Gabriel;

  typedef typename Kernel::Point_d Point_d;

  // Type required to compute squared radius, or side of bounded sphere on a vector of points.
  typedef typename std::vector<Point_d> Vector_of_CGAL_points;

  // Vertex_iterator type from CGAL.
  typedef typename Delaunay_triangulation::Vertex_iterator CGAL_vertex_iterator;

  // size_type type from CGAL.
  typedef typename Delaunay_triangulation::size_type size_type;

  // Double map type to switch from CGAL vertex iterator to simplex tree vertex handle and vice versa.
  typedef typename std::map< CGAL_vertex_iterator, Vertex_handle > Map_vertex_iterator_to_handle;
  typedef typename std::map< Vertex_handle, CGAL_vertex_iterator > Map_vertex_handle_to_iterator;

 private:
  /** \brief Map to switch from CGAL vertex iterator to simplex tree vertex handle.*/
  Map_vertex_iterator_to_handle vertex_iterator_to_handle_;
  /** \brief Map to switch from simplex tree vertex handle to CGAL vertex iterator.*/
  Map_vertex_handle_to_iterator vertex_handle_to_iterator_;
  /** \brief Pointer on the CGAL Delaunay triangulation.*/
  Delaunay_triangulation* triangulation_;
  /** \brief Kernel for triangulation_ functions access.*/
  Kernel kernel_;
  /** \brief Maximum value for alpha square.*/
  Filtration_value max_alpha_square_;

 public:
  /** \brief Alpha_complex constructor from an OFF file name.
   * Uses the Delaunay_triangulation_off_reader to construct the Delaunay triangulation required to initialize 
   * the Alpha_complex.
   *
   * @param[in] off_file_name OFF file [path and] name.
   */
  Alpha_complex(const std::string& off_file_name, Filtration_value max_alpha_square)
      : triangulation_(nullptr),
      max_alpha_square_(max_alpha_square) {
    Gudhi::Delaunay_triangulation_off_reader<Delaunay_triangulation> off_reader(off_file_name);
    if (!off_reader.is_valid()) {
      std::cerr << "Alpha_complex - Unable to read file " << off_file_name << std::endl;
      exit(-1);  // ----- >>
    }
    triangulation_ = off_reader.get_complex();
    init();
  }

  /** \brief Alpha_complex constructor from a Delaunay triangulation.
   *
   * @param[in] triangulation_ptr Pointer on a Delaunay triangulation.
   */
  Alpha_complex(Delaunay_triangulation* triangulation_ptr, Filtration_value max_alpha_square)
      : triangulation_(triangulation_ptr),
      max_alpha_square_(max_alpha_square) {
    init();
  }

  /** \brief Alpha_complex constructor from a list of points.
   * Uses the Delaunay_triangulation_off_reader to construct the Delaunay triangulation required to initialize 
   * the Alpha_complex.
   *
   * @param[in] dimension Dimension of points to be inserted.
   * @param[in] size Number of points to be inserted.
   * @param[in] firstPoint Iterator on the first point to be inserted.
   * @param[in] last Point Iterator on the last point to be inserted.
   */
  template<typename ForwardIterator >
  Alpha_complex(int dimension, size_type size, ForwardIterator firstPoint, ForwardIterator lastPoint,
                Filtration_value max_alpha_square)
      : triangulation_(nullptr),
      max_alpha_square_(max_alpha_square) {
    triangulation_ = new Delaunay_triangulation(dimension);
    size_type inserted = triangulation_->insert(firstPoint, lastPoint);
    if (inserted != size) {
      std::cerr << "Alpha_complex - insertion failed " << inserted << " != " << size << std::endl;
      exit(-1);  // ----- >>
    }
    init();
  }

  /** \brief Alpha_complex destructor from a Delaunay triangulation.
   *
   * @warning Deletes the Delaunay triangulation.
   */
  ~Alpha_complex() {
    delete triangulation_;
  }

  /** \brief get_point returns the point corresponding to the vertex given as parameter.
   *
   * @param[in] vertex Vertex handle of the point to retrieve.
   * @return The founded point.
   */
  Point_d get_point(Vertex_handle vertex) {
    Point_d point(dimension());
    try {
      if (vertex_handle_to_iterator_[vertex] != nullptr) {
        point = vertex_handle_to_iterator_[vertex]->point();
      }
    } catch (...) {
      std::cerr << "Alpha_complex - getPoint not found on vertex " << vertex << std::endl;
    }
    return point;
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
    if (triangulation_ == nullptr) {
      std::cerr << "Alpha_complex init - Cannot init from a NULL triangulation" << std::endl;
      return;  // ----- >>
    }
    if (triangulation_->number_of_vertices() < 1) {
      std::cerr << "Alpha_complex init - Cannot init from a triangulation without vertices" << std::endl;
      return;  // ----- >>
    }
    if (triangulation_->maximal_dimension() < 1) {
      std::cerr << "Alpha_complex init - Cannot init from a zero-dimension triangulation" << std::endl;
      return;  // ----- >>
    }
    if (num_vertices() > 0) {
      std::cerr << "Alpha_complex init - Cannot init twice" << std::endl;
      return;  // ----- >>
    }

    set_dimension(triangulation_->maximal_dimension());

    // --------------------------------------------------------------------------------------------
    // double map to retrieve simplex tree vertex handles from CGAL vertex iterator and vice versa
    // Start to insert at handle = 0 - default integer value
    Vertex_handle vertex_handle = Vertex_handle();
    // Loop on triangulation vertices list
    for (CGAL_vertex_iterator vit = triangulation_->vertices_begin(); vit != triangulation_->vertices_end(); ++vit) {
#ifdef DEBUG_TRACES
      std::cout << "Vertex insertion - " << vertex_handle << " -> " << vit->point() << std::endl;
#endif  // DEBUG_TRACES
      vertex_iterator_to_handle_[vit] = vertex_handle;
      vertex_handle_to_iterator_[vertex_handle] = vit;
      vertex_handle++;
    }
    // --------------------------------------------------------------------------------------------

    Filtration_value filtration_max = 0.0;
    // --------------------------------------------------------------------------------------------
    // Simplex_tree construction from loop on triangulation finite full cells list
    for (auto cit = triangulation_->finite_full_cells_begin(); cit != triangulation_->finite_full_cells_end(); ++cit) {
      Vector_vertex vertex_full_cell;
#ifdef DEBUG_TRACES
      std::cout << "Simplex_tree insertion ";
#endif  // DEBUG_TRACES
      for (auto vit = cit->vertices_begin(); vit != cit->vertices_end(); ++vit) {
        if (*vit != nullptr) {
#ifdef DEBUG_TRACES
          std::cout << " " << vertex_iterator_to_handle_[*vit];
#endif  // DEBUG_TRACES
          // Vector of vertex construction for simplex_tree structure
          vertex_full_cell.push_back(vertex_iterator_to_handle_[*vit]);
        }
      }
#ifdef DEBUG_TRACES
      std::cout << std::endl;
#endif  // DEBUG_TRACES

      Simplex_tree<> full_cell;
      full_cell.set_dimension(triangulation_->maximal_dimension());
      // Create a simplex tree containing only one of the full cells
      Simplex_result insert_result = full_cell.insert_simplex_and_subfaces(vertex_full_cell);
      if (!insert_result.second) {
        std::cerr << "Alpha_complex::init insert_simplex_and_subfaces failed" << std::endl;
        exit(-1);  // ----->>
      }
      // +++ For i : d -> 0
      for (int fc_decr_dim = full_cell.dimension(); (fc_decr_dim >= 0); fc_decr_dim--) {
        // +++ Foreach Sigma of dim i
        // No need to skip this loop in case alpha²(Sigma) > max_alpha_square_ because of
        // if (fc_decr_dim == f_simplex_dim) which means "only for a full cell"
        for (auto fc_simplex : full_cell.skeleton_simplex_range(fc_decr_dim)) {
          int f_simplex_dim = full_cell.dimension(fc_simplex);
          if (fc_decr_dim == f_simplex_dim) {
            Vector_of_CGAL_points pointVector;
            Vector_vertex current_vertex;
#ifdef DEBUG_TRACES
            std::cout << "Sigma of dim " << fc_decr_dim << " is";
#endif  // DEBUG_TRACES
            for (auto vertex : full_cell.simplex_vertex_range(fc_simplex)) {
              pointVector.push_back(get_point(vertex));
              current_vertex.push_back(vertex);
#ifdef DEBUG_TRACES
              std::cout << " " << vertex;
#endif  // DEBUG_TRACES
            }
#ifdef DEBUG_TRACES
            std::cout << std::endl;
#endif  // DEBUG_TRACES
            Simplex_handle sigma_handle = find(current_vertex);
            bool skip_propagation = false;
            // +++ If filt(Sigma) is NaN : filt(Sigma) = alpha²(Sigma)
            if ((sigma_handle == null_simplex()) || isnan(filtration(sigma_handle))) {
              Filtration_value alpha_complex_filtration = compute_alpha_square(pointVector.begin(), pointVector.end(),
                                                                               f_simplex_dim);
              if (alpha_complex_filtration <= max_alpha_square_) {
                // Only insert Sigma in Simplex tree if alpha²(Sigma) <= max_alpha_square_
                if (sigma_handle == null_simplex()) {
#ifdef DEBUG_TRACES
                  std::cout << "Alpha_complex::init Sigma not found" << std::endl;
#endif  // DEBUG_TRACES
                  insert_result = insert_simplex(current_vertex, std::numeric_limits<double>::quiet_NaN());
                  if (!insert_result.second) {
                    std::cerr << "Alpha_complex::init insert_simplex failed" << std::endl;
                    exit(-1);  // ----->>
                  }
                  // Sigma is the new inserted simplex handle
                  sigma_handle = insert_result.first;
                }
#ifdef DEBUG_TRACES
                std::cout << "Alpha_complex::init filtration = " << alpha_complex_filtration << std::endl;
#endif  // DEBUG_TRACES
                assign_filtration(sigma_handle, alpha_complex_filtration);
                filtration_max = fmax(filtration_max, alpha_complex_filtration);
              } else {
                // if alpha²(Sigma) > max_alpha_square_ skip propagation
                skip_propagation = true;
#ifdef DEBUG_TRACES
                std::cout << "Alpha_complex::init skip propagation on this full cell" << std::endl;
#endif  // DEBUG_TRACES
              }
            }  // --- If filt(Sigma) is NaN : filt(Sigma) = alpha(Sigma)
            if ((filtration(sigma_handle) <= max_alpha_square_) && !skip_propagation) {
              // Propagate alpha filtration value in Simplex tree if alpha²(Sigma) <= max_alpha_square_
              // in case Sigma is not found AND not inserted (alpha_complex_filtration > max_alpha_square_),
              // filtration(null_simplex()) returns INFINITY => no propagation
              propagate_alpha_filtration(full_cell, fc_simplex, fc_decr_dim, sigma_handle);
            }
          }
        }  // --- Foreach Sigma of dim i
      }  // --- For i : d -> 0
    }
    // --------------------------------------------------------------------------------------------

#ifdef DEBUG_TRACES
    std::cout << "filtration_max=" << filtration_max << std::endl;
#endif  // DEBUG_TRACES
    set_filtration(filtration_max);
  }

  template<typename ForwardIterator >
  Filtration_value compute_alpha_square(ForwardIterator firstPoint, ForwardIterator lastPoint, int f_simplex_dim) {
    Filtration_value alpha_square_value = 0.0;
    // No need to compute squared_radius on a single point - alpha is 0.0
    if (f_simplex_dim > 0) {
      // squared_radius function initialization
      Squared_Radius squared_radius = kernel_.compute_squared_radius_d_object();

      alpha_square_value = squared_radius(firstPoint, lastPoint);
    }
    return alpha_square_value;
  }

  void propagate_alpha_filtration(Simplex_tree& full_cell, Simplex_handle fc_simplex, int fc_decr_dim,
                                  Simplex_handle sigma_handle) {
    // ### Foreach Tau face of Sigma
    for (auto f_boundary : full_cell.boundary_simplex_range(fc_simplex)) {
#ifdef DEBUG_TRACES
      std::cout << " | --------------------------------------------------" << std::endl;
      std::cout << " | Tau ";
#endif  // DEBUG_TRACES
      Vector_vertex tau_vertex;
      for (auto vertex : full_cell.simplex_vertex_range(f_boundary)) {
        tau_vertex.push_back(vertex);
#ifdef DEBUG_TRACES
        std::cout << vertex << " ";
#endif  // DEBUG_TRACES
      }
#ifdef DEBUG_TRACES
      std::cout << "is a face of Sigma" << std::endl;
#endif  // DEBUG_TRACES
      Simplex_handle tau_handle = find(tau_vertex);
      // ### If filt(Tau) is not NaN

      if ((tau_handle != null_simplex()) && (!isnan(filtration(tau_handle)))) {
        // ### filt(Tau) = fmin(filt(Tau), filt(Sigma))
        Filtration_value alpha_complex_filtration = fmin(filtration(tau_handle), filtration(sigma_handle));
        assign_filtration(tau_handle, alpha_complex_filtration);
        // No need to check for filtration_max, alpha_complex_filtration is a min of an existing filtration value
#ifdef DEBUG_TRACES
        std::cout << " | filt(Tau) = fmin(filt(Tau), filt(Sigma)) = " << alpha_complex_filtration << std::endl;
#endif  // DEBUG_TRACES
      } else {
        // No need to compute is_gabriel for dimension <= 2
        // i.e. : Sigma = (3,1) => Tau = 1
        if (fc_decr_dim > 1) {
          // insert the Tau points in a vector for is_gabriel function
          Vector_of_CGAL_points pointVector;
          Vertex_handle vertexForGabriel = Vertex_handle();
          for (auto vertex : full_cell.simplex_vertex_range(f_boundary)) {
            pointVector.push_back(get_point(vertex));
          }
          // Retrieve the Sigma point that is not part of Tau - parameter for is_gabriel function
          for (auto vertex : simplex_vertex_range(sigma_handle)) {
            if (std::find(pointVector.begin(), pointVector.end(), get_point(vertex)) == pointVector.end()) {
              // vertex is not found in Tau
              vertexForGabriel = vertex;
              // No need to continue loop
              break;
            }
          }
          // is_gabriel function initialization
          Is_Gabriel is_gabriel = kernel_.side_of_bounded_sphere_d_object();
          bool is_gab = is_gabriel(pointVector.begin(), pointVector.end(), get_point(vertexForGabriel))
              != CGAL::ON_BOUNDED_SIDE;
#ifdef DEBUG_TRACES
          std::cout << " | Tau is_gabriel(Sigma)=" << is_gab << " - vertexForGabriel=" << vertexForGabriel << std::endl;
#endif  // DEBUG_TRACES
          // ### If Tau is not Gabriel of Sigma
          if (false == is_gab) {
            if (tau_handle == null_simplex()) {
#ifdef DEBUG_TRACES
              std::cout << " | Tau not found" << std::endl;
#endif  // DEBUG_TRACES
              // in case Tau is not yet created
              Simplex_result insert_result = insert_simplex(tau_vertex, std::numeric_limits<double>::quiet_NaN());
              if (!insert_result.second) {
                std::cerr << "Alpha_complex::propagate_alpha_filtration insert_simplex failed" << std::endl;
                exit(-1);  // ----->>
              }
              // Sigma is the new inserted simplex handle
              tau_handle = insert_result.first;
            }
            // ### filt(Tau) = filt(Sigma)
            assign_filtration(tau_handle, filtration(sigma_handle));
            // No need to check for filtration_max, alpha_complex_filtration is an existing filtration value
#ifdef DEBUG_TRACES
            std::cout << " | filt(Tau) = filt(Sigma) = " << filtration(sigma_handle) << std::endl;
#endif  // DEBUG_TRACES
          }
        }
      }
    }
  }
};

}  // namespace alphacomplex

}  // namespace Gudhi

#endif  // ALPHA_COMPLEX_H_
