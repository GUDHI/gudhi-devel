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
#ifndef DELAUNAY_TRIANGULATION_OFF_IO_H_
#define DELAUNAY_TRIANGULATION_OFF_IO_H_

#include <string>
#include <vector>
#include <fstream>
#include <map>

#include "gudhi/Off_reader.h"

namespace Gudhi {

/** 
 * \class Delaunay_triangulation_off_visitor_reader Delaunay_triangulation_off_io.h gudhi/Delaunay_triangulation_off_io.h
 * \brief OFF file visitor implementation according to Off_reader in order to construct a CGAL Delaunay triangulation.
 * 
 * For more informations on CGAL Delaunay triangulation, please refer to the corresponding chapter in page
 *  http://doc.cgal.org/latest/Triangulation/
 */
template<typename Complex>
class Delaunay_triangulation_off_visitor_reader {
 private:
  Complex* complex_;
  typedef typename Complex::Point Point;
  std::vector<Point> point_cloud;

 public:
  // TODO(VR) : Pass a Complex as a parameter is required, even if not used. Otherwise, compilation is KO.

  /** \brief Delaunay_triangulation_off_visitor_reader constructor
   *
   * @param[in] complex_ptr_ pointer on a Delaunay triangulation.
   */
  Delaunay_triangulation_off_visitor_reader(Complex* complex_ptr_)
      : complex_(nullptr) { }

  /** \brief Off_reader visitor init implementation. 
   * 
   * The init parameters are set from OFF file header.
   * Dimension value is required in order to construct Delaunay triangulation.
   *
   * @param[in] dim space dimension of vertices.
   * @param[in] num_vertices number of vertices in the OFF file (not used).
   * @param[in] num_faces number of faces in the OFF file (not used).
   * @param[in] num_edges number of edges in the OFF file (not used).
   */
  void init(int dim, int num_vertices, int num_faces, int num_edges) {
#ifdef DEBUG_TRACES
    std::cout << "Delaunay_triangulation_off_visitor_reader::init - dim=" << dim << " - num_vertices=" <<
        num_vertices << " - num_faces=" << num_faces << " - num_edges=" << num_edges << std::endl;
#endif  // DEBUG_TRACES
    if (num_faces > 0) {
      std::cerr << "Delaunay_triangulation_off_visitor_reader::init faces are not taken into account from OFF " <<
          "file for Delaunay triangulation - faces are computed.\n";
    }
    if (num_edges > 0) {
      std::cerr << "Delaunay_triangulation_off_visitor_reader::init edges are not taken into account from OFF " <<
          "file for Delaunay triangulation - edges are computed.\n";
    }
    // Complex construction with dimension from file
    complex_ = new Complex(dim);
  }

  /** \brief Off_reader visitor point implementation. 
   * 
   * The point function is called on each vertex line from OFF file.
   * This function inserts the vertex in the Delaunay triangulation.
   *
   * @param[in] point vector of vertex coordinates.
   */
  void point(const std::vector<double>& point) {
#ifdef DEBUG_TRACES
    std::cout << "Delaunay_triangulation_off_visitor_reader::point ";
    for (auto coordinate : point) {
      std::cout << coordinate << " | ";
    }
    std::cout << std::endl;
#endif  // DEBUG_TRACES
    // Fill the point cloud
    // VR: complex_->insert(Point(point.size(), point.begin(), point.end()));
    point_cloud.push_back(Point(point.size(), point.begin(), point.end()));
  }

  // Off_reader visitor maximal_face implementation - not used
  void maximal_face(const std::vector<int>& face) {
    // For Delaunay Triangulation, only points are read
  }

  // Off_reader visitor done implementation
  void done() {
    // It is advised to insert all the points at a time in a Delaunay Triangulation because points are sorted at the
    // beginning of the insertion
    complex_->insert(point_cloud.begin(), point_cloud.end());
  }

  /** \brief Returns the constructed Delaunay triangulation.
   *
   * @return A pointer on the Delaunay triangulation.
   * 
   * @warning The returned pointer can be nullptr.
   */
  Complex* get_complex() const {
    return complex_;
  }
};

/** 
 * \class Delaunay_triangulation_off_reader Delaunay_triangulation_off_io.h gudhi/Delaunay_triangulation_off_io.h
 * \brief OFF file reader implementation in order to construct a Delaunay triangulation.
 * 
 * This class is using the Delaunay_triangulation_off_visitor_reader to visit the OFF file according to Off_reader.
 * 
 * For more informations on CGAL Delaunay triangulation, please refer to the corresponding chapter in page
 *  http://doc.cgal.org/latest/Triangulation/
 * 
 * \section Example
 *
 * This example loads points from an OFF file and builds the Delaunay triangulation.
 * Then, it is asked to display the number of vertices and finites full cells from the Delaunay triangulation.
 * 
 * \include Delaunay_triangulation_off_rw.cpp
 * 
 * When launching:
 * 
 * \code $> ./dtoffrw ../../data/points/alphacomplexdoc.off triangulated.off
 * \endcode
 *
 * the program output is:
 * 
 * \include dtoffrw_alphashapedoc_result.txt
 */
template<typename Complex>
class Delaunay_triangulation_off_reader {
 public:
  /** \brief Reads the OFF file and constructs the Delaunay triangulation from the points
   * that are in the OFF file.
   *
   * @param[in] name_file OFF file to read.
   *
   * @warning Check with is_valid() function to see if read operation was successful.
   */
  Delaunay_triangulation_off_reader(const std::string & name_file)
  : valid_(false) {
    std::ifstream stream(name_file);
    if (stream.is_open()) {
      Delaunay_triangulation_off_visitor_reader<Complex> off_visitor(complex_);
      Off_reader off_reader(stream);
      valid_ = off_reader.read(off_visitor);
      if (valid_) {
        complex_ = off_visitor.get_complex();
        if (complex_ == nullptr) {
          std::cerr << "Delaunay_triangulation_off_reader::Delaunay_triangulation_off_reader off_visitor returns " <<
              "an empty pointer\n";
          valid_ = false;
        }
      }
    } else {
      std::cerr << "Delaunay_triangulation_off_reader::Delaunay_triangulation_off_reader could not open file " <<
          name_file << "\n";
    }
  }

  /** \brief Returns if the OFF file read operation was successful or not.
   *
   * @return OFF file read status.
   */
  bool is_valid() const {
    return valid_;
  }

  /** \brief Returns the constructed Delaunay triangulation.
   *
   * @return A pointer on the Delaunay triangulation.
   * 
   * @warning The returned pointer can be nullptr.
   */
  Complex* get_complex() const {
    if (valid_)
      return complex_;
    return nullptr;
  }

 private:
  /** \brief OFF file read status.*/
  bool valid_;
  /** \brief A pointer on the Delaunay triangulation.*/
  Complex* complex_;
};

/**
 * \class Delaunay_triangulation_off_writer Delaunay_triangulation_off_io.h gudhi/Delaunay_triangulation_off_io.h
 * \brief OFF file writer from a Delaunay triangulation.
 * 
 * This class constructs the OFF file header according to http://www.geomview.org/docs/html/OFF.html
 * 
 * The header is followed by the list of points coordinates (Delaunay triangulation vertices)
 * 
 * And finally is followed by the list of faces (Delaunay triangulation finite full cells)
 * 
 * For more informations on CGAL Delaunay triangulation, please refer to the corresponding chapter in page
 * http://doc.cgal.org/latest/Triangulation/
 * 
 * \section Example
 *
 * This example loads points from an OFF file and builds the Delaunay triangulation.
 * Then, the Delaunay triangulation is saved in a new file including the triangulation as a list of faces.
 * 
 * \include Delaunay_triangulation_off_rw.cpp
 * 
 * When launching:
 * 
 * \code $> ./dtoffrw ../../data/points/alphashapedoc.off triangulated.off
 * \endcode
 *
 * The result will be an OFF file of dimension 2 with the 7 points from alphashapedoc.off followed by the 6
 * triangulations of dimension 3 (the first value on each faces):
 * \include dtoffrw_alphashapedoc_result.off
 */
template<typename Complex>
class Delaunay_triangulation_off_writer {
 public:
  typedef typename Complex::Point Point;

  /** \brief Writes the OFF file from the Delaunay triangulation
   *
   * @param[in] name_file OFF file to write.
   * @param[in] complex_ptr pointer on a Delaunay triangulation.
   *
   * @warning Check with is_valid() function to see if write operation was successful.
   */
  Delaunay_triangulation_off_writer(const std::string & name_file, Complex* complex_ptr)
      : valid_(false) {
    std::ofstream stream(name_file);
    if (stream.is_open()) {
      if (complex_ptr->current_dimension() == 3) {
        // OFF header
        stream << "OFF" << std::endl;
        // no endl on next line - don't know why...
        stream << complex_ptr->number_of_vertices() << " " << complex_ptr->number_of_finite_full_cells() << " 0";
      } else {
        // nOFF header
        stream << "nOFF" << std::endl;
        // no endl on next line - don't know why...
        stream << complex_ptr->current_dimension() << " " << complex_ptr->number_of_vertices() << " " <<
            complex_ptr->number_of_finite_full_cells() << " 0";
      }

      // bimap to retrieve vertex handles from points and vice versa
      std::map< Point, int > points_to_vh;
      // Start to insert at default handle value
      int vertex_handle = int();

      // Points list
      for (auto vit = complex_ptr->vertices_begin(); vit != complex_ptr->vertices_end(); ++vit) {
        for (auto Coord = vit->point().cartesian_begin(); Coord != vit->point().cartesian_end(); ++Coord) {
          stream << *Coord << " ";
        }
        stream << std::endl;
        points_to_vh[vit->point()] = vertex_handle;
        vertex_handle++;
      }

      for (auto cit = complex_ptr->finite_full_cells_begin(); cit != complex_ptr->finite_full_cells_end(); ++cit) {
        stream << std::distance(cit->vertices_begin(), cit->vertices_end()) << " ";
        for (auto vit = cit->vertices_begin(); vit != cit->vertices_end(); ++vit) {
          stream << points_to_vh[(*vit)->point()] - 1 << " ";
        }
        stream << std::endl;
      }
      stream.close();
      valid_ = true;
    } else {
      std::cerr << "Delaunay_triangulation_off_writer::Delaunay_triangulation_off_writer could not open file " <<
          name_file << "\n";
    }
  }

  /** \brief Returns if the OFF write operation was successful or not.
   *
   * @return OFF file write status.
   */
  bool is_valid() const {
    return valid_;
  }

 private:
  /* \brief OFF file write status. */
  bool valid_;
};

}  // namespace Gudhi

#endif  // DELAUNAY_TRIANGULATION_OFF_IO_H_
