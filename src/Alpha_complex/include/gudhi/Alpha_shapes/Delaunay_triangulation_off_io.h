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
#ifndef SRC_ALPHA_SHAPES_INCLUDE_GUDHI_ALPHA_SHAPES_DELAUNAY_TRIANGULATION_OFF_IO_H_
#define SRC_ALPHA_SHAPES_INCLUDE_GUDHI_ALPHA_SHAPES_DELAUNAY_TRIANGULATION_OFF_IO_H_

#include <string>
#include <vector>
#include <fstream>
#include <map>

#include "gudhi/Off_reader.h"

namespace Gudhi {

namespace alphacomplex {

/**
 *@brief Off reader visitor with flag that can be passed to Off_reader to read a Delaunay_triangulation_complex.
 */
template<typename Complex>
class Delaunay_triangulation_off_visitor_reader {
 private:
  Complex* _complex;
  typedef typename Complex::Point Point;

 public:

  // Pass a Complex as a parameter is required, even if not used. Otherwise, compilation is KO.
  Delaunay_triangulation_off_visitor_reader(Complex* _complex_ptr)
      : _complex(nullptr) { }

  void init(int dim, int num_vertices, int num_faces, int num_edges) {
#ifdef DEBUG_TRACES
    std::cout << "Delaunay_triangulation_off_visitor_reader::init - dim=" << dim << " - num_vertices=" <<
        num_vertices << " - num_faces=" << num_faces << " - num_edges=" << num_edges << std::endl;
#endif  // DEBUG_TRACES
    if (num_faces > 0) {
      std::cerr << "Delaunay_triangulation_off_visitor_reader::init faces are not taken into account from OFF " <<
          "file for Delaunay triangulation - faces are computed." << std::endl;
    }
    if (num_edges > 0) {
      std::cerr << "Delaunay_triangulation_off_visitor_reader::init edges are not taken into account from OFF " <<
          "file for Delaunay triangulation - edges are computed." << std::endl;
    }
    // Complex construction with dimension from file
    _complex = new Complex(dim);
  }

  void point(const std::vector<double>& point) {
#ifdef DEBUG_TRACES
    std::cout << "Delaunay_triangulation_off_visitor_reader::point ";
    for (auto coordinate : point) {
      std::cout << coordinate << " | ";
    }
    std::cout << std::endl;
#endif  // DEBUG_TRACES
    _complex->insert(Point(point.size(), point.begin(), point.end()));
  }

  void maximal_face(const std::vector<int>& face) {
    // For Delaunay Triangulation, only points are read
#ifdef DEBUG_TRACES
    std::cout << "Delaunay_triangulation_off_visitor_reader::face ";
    for (auto vertex : face) {
      std::cout << vertex << " | ";
    }
    std::cout << std::endl;
#endif  // DEBUG_TRACES
  }

  void done() {
#ifdef DEBUG_TRACES
    std::cout << "Delaunay_triangulation_off_visitor_reader::done" << std::endl;
#endif  // DEBUG_TRACES
  }

  Complex* get_complex() {
    return _complex;
  }
};

/**
 *@brief Class that allows to load a Delaunay_triangulation_complex from an off file.
 */
template<typename Complex>
class Delaunay_triangulation_off_reader {
 public:

  /**
   * name_file : file to read
   * read_complex : complex that will receive the file content
   * read_only_points : specify true if only the points must be read
   */
  Delaunay_triangulation_off_reader(const std::string & name_file) : valid_(false) {
    std::ifstream stream(name_file);
    if (stream.is_open()) {
      Delaunay_triangulation_off_visitor_reader<Complex> off_visitor(_complex);
      Off_reader off_reader(stream);
      valid_ = off_reader.read(off_visitor);
      if (valid_) {
        _complex = off_visitor.get_complex();
      }
    } else {
      std::cerr << "Delaunay_triangulation_off_reader::Delaunay_triangulation_off_reader could not open file " <<
          name_file << std::endl;
    }

  }

  /**
   * return true if reading did not meet problems.
   */
  bool is_valid() const {
    return valid_;
  }

  Complex* get_complex() {
    if (valid_)
      return _complex;
    return nullptr;

  }

 private:
  bool valid_;
  Complex* _complex;
};

template<typename Complex>
class Delaunay_triangulation_off_writer {
 public:
  typedef typename Complex::Point Point;

  /**
   * name_file : file where the off will be written
   * save_complex : complex that be outputted in the file
   * for now only save triangles.
   */
  Delaunay_triangulation_off_writer(const std::string & name_file, Complex* complex_ptr) {
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
        std::vector<int> vertexVector;
        stream << std::distance(cit->vertices_begin(), cit->vertices_end()) << " ";
        for (auto vit = cit->vertices_begin(); vit != cit->vertices_end(); ++vit) {
          stream << points_to_vh[(*vit)->point()] << " ";
        }
        stream << std::endl;
      }
      stream.close();
    } else {
      std::cerr << "Delaunay_triangulation_off_writer::Delaunay_triangulation_off_writer could not open file " <<
          name_file << std::endl;
    }
  }
};

} // namespace alphacomplex

} // namespace Gudhi

#endif  // SRC_ALPHA_SHAPES_INCLUDE_GUDHI_ALPHA_SHAPES_DELAUNAY_TRIANGULATION_OFF_IO_H_
