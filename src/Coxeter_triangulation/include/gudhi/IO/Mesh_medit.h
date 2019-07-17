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

#ifndef IO_MESH_MEDIT_H_
#define IO_MESH_MEDIT_H_

/* \class Mesh_medit
 * \brief Structure to store a mesh that can be output in Medit .mesh file format
 *  using the output_meshes_to_medit method.
 *
 * \ingroup coxeter_triangulation
 */
struct Mesh_medit {
  /** \brief Type of a range of vertices. */
  typedef std::vector<Eigen::VectorXd> Vertex_points;
  /** \brief Type of a mesh element.
   *  A pair consisting of a vector of vertex indices of type std::size_t
   *  and of an integer that represents the common reference number for
   *  the mesh elements of this type. */
  typedef std::pair<std::vector<std::size_t>, std::size_t> Mesh_element;
  /** \brief Type of a range of mesh elements. */
  typedef std::vector<Mesh_element> Mesh_elements;
  /** \brief Type of a range of scalar field . */
  typedef std::vector<double> Scalar_field_range;

  /** \brief Range of vertices of type Eigen::VectorXd to output. */
  Vertex_points vertex_points;
  /** \brief Range of edges. */
  Mesh_elements edges;
  /** \brief Range of triangles. */
  Mesh_elements triangles;
  /** \brief Range of tetrahedra. */
  Mesh_elements tetrahedra;
  /** \brief Range of scalar values over triangles. */
  Scalar_field_range triangles_scalar_range;
  /** \brief Range of scalar values over tetrahedra. */
  Scalar_field_range tetrahedra_scalar_range;
};

#endif
