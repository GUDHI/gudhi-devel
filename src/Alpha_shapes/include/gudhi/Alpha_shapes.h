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

#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/Epick_d.h>
#include <CGAL/algorithm.h>
#include <CGAL/assertions.h>

#include <iostream>
#include <iterator>
#include <vector>
#include <string>

namespace Gudhi {

namespace alphashapes {

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

  // From CGAL
  /** \brief Kernel for the Delaunay_triangulation.
   * Dimension can be set dynamically.
   */
  typedef CGAL::Epick_d< CGAL::Dynamic_dimension_tag > Kernel;
  /** \brief Delaunay_triangulation type required to create an alpha-shape.
   */
  typedef CGAL::Delaunay_triangulation<Kernel> Delaunay_triangulation;

 private:
  /** \brief Upper bound on the simplex tree of the simplicial complex.*/
  Gudhi::Simplex_tree<> _st;

 public:

  Alpha_shapes(std::string off_file_name, int dimension) {
    Delaunay_triangulation dt(dimension);
    Gudhi::alphashapes::Delaunay_triangulation_off_reader<Delaunay_triangulation>
        off_reader(off_file_name, dt, true, true);
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
  Alpha_shapes(T triangulation) {
    init<T>(triangulation);
  }

  ~Alpha_shapes() { }

 private:

  template<typename T>
  void init(T triangulation) {
    _st.set_dimension(triangulation.maximal_dimension());
    _st.set_filtration(0.0);
    // triangulation points list
    for (auto vit = triangulation.finite_vertices_begin();
         vit != triangulation.finite_vertices_end(); ++vit) {
      typeVectorVertex vertexVector;
      Vertex_handle vertexHdl = std::distance(triangulation.finite_vertices_begin(), vit);
      vertexVector.push_back(vertexHdl);

      // Insert each point in the simplex tree
      _st.insert_simplex(vertexVector, 0.0);

#ifdef DEBUG_TRACES
      std::cout << "P" << vertexHdl << ":";
      for (auto Coord = vit->point().cartesian_begin(); Coord != vit->point().cartesian_end(); ++Coord) {
        std::cout << *Coord << " ";
      }
      std::cout << std::endl;
#endif  // DEBUG_TRACES
    }
    // triangulation finite full cells list
    for (auto cit = triangulation.finite_full_cells_begin();
         cit != triangulation.finite_full_cells_end(); ++cit) {
      typeVectorVertex vertexVector;
      for (auto vit = cit->vertices_begin(); vit != cit->vertices_end(); ++vit) {
        // Vertex handle is distance - 1
        Vertex_handle vertexHdl = std::distance(triangulation.vertices_begin(), *vit) - 1;
        vertexVector.push_back(vertexHdl);
      }
      // Insert each point in the simplex tree
      _st.insert_simplex_and_subfaces(vertexVector, 0.0);

#ifdef DEBUG_TRACES
      std::cout << "C" << std::distance(triangulation.finite_full_cells_begin(), cit) << ":";
      for (auto value : vertexVector) {
        std::cout << value << ' ';
      }
      std::cout << std::endl;
#endif  // DEBUG_TRACES
    }
  }

 public:

  /** \brief Returns the number of vertices in the complex. */
  size_t num_vertices() {
    return _st.num_vertices();
  }

  /** \brief Returns the number of simplices in the complex.
   *
   * Does not count the empty simplex. */
  unsigned num_simplices() const {
    return _st.num_simplices();
  }

  /** \brief Returns an upper bound on the dimension of the simplicial complex. */
  int dimension() {
    return _st.dimension();
  }

  /** \brief Returns an upper bound of the filtration values of the simplices. */
  Filtration_value filtration() {
    return _st.filtration();
  }

  friend std::ostream& operator<<(std::ostream& os, const Alpha_shapes& alpha_shape) {
    // TODO: Program terminated with signal SIGABRT, Aborted - Maybe because of copy constructor
    Gudhi::Simplex_tree<> st = alpha_shape._st;
    os << st << std::endl;
    return os;
  }
};

} // namespace alphashapes

} // namespace Gudhi

#endif  // SRC_ALPHA_SHAPES_INCLUDE_GUDHI_ALPHA_SHAPES_H_
