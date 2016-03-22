/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014  INRIA Sophia Antipolis-Mediterranee (France)
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
#ifndef SKELETON_BLOCKER_SKELETON_BLOCKER_OFF_IO_H_
#define SKELETON_BLOCKER_SKELETON_BLOCKER_OFF_IO_H_

#include <gudhi/Off_reader.h>

#include <string>
#include <vector>
#include <map>

namespace Gudhi {

namespace skbl {

/**
 *@brief Off reader visitor that can be passed to Off_reader to read a Skeleton_blocker_complex.
 */
template<typename Complex>
class Skeleton_blocker_off_flag_visitor_reader {
  Complex& complex_;
  typedef typename Complex::Vertex_handle Vertex_handle;
  typedef typename Complex::Point Point;

  const bool load_only_points_;

 public:
  explicit Skeleton_blocker_off_flag_visitor_reader(Complex& complex, bool load_only_points = false) :
      complex_(complex),
      load_only_points_(load_only_points) { }

  void init(int dim, int num_vertices, int num_faces, int num_edges) {
    // todo do an assert to check that this number are correctly read
    // todo reserve size for vector points
  }

  void point(const std::vector<double>& point) {
    complex_.add_vertex(Point(point.begin(), point.end()));
  }

  void maximal_face(const std::vector<int>& face) {
    if (!load_only_points_) {
      for (size_t i = 0; i < face.size(); ++i)
        for (size_t j = i + 1; j < face.size(); ++j) {
          complex_.add_edge_without_blockers(Vertex_handle(face[i]), Vertex_handle(face[j]));
        }
    }
  }

  void done() { }
};

/**
 *@brief Off reader visitor that can be passed to Off_reader to read a Skeleton_blocker_complex.
 */
template<typename Complex>
class Skeleton_blocker_off_visitor_reader {
  Complex& complex_;
  typedef typename Complex::Vertex_handle Vertex_handle;
  typedef typename Complex::Simplex Simplex;
  typedef typename Complex::Point Point;

  const bool load_only_points_;
  std::vector<Point> points_;
  std::vector<Simplex> maximal_faces_;

 public:
  explicit Skeleton_blocker_off_visitor_reader(Complex& complex, bool load_only_points = false) :
      complex_(complex),
      load_only_points_(load_only_points) { }

  void init(int dim, int num_vertices, int num_faces, int num_edges) {
    maximal_faces_.reserve(num_faces);
    points_.reserve(num_vertices);
  }

  void point(const std::vector<double>& point) {
    points_.emplace_back(point.begin(), point.end());
  }

  void maximal_face(const std::vector<int>& face) {
    if (!load_only_points_) {
      Simplex s;
      for (auto x : face)
        s.add_vertex(Vertex_handle(x));
      maximal_faces_.emplace_back(s);
    }
  }

  void done() {
    complex_ = make_complex_from_top_faces<Complex>(maximal_faces_.begin(), maximal_faces_.end(),
                                           points_.begin(), points_.end() );
  }
};

/**
 *@brief Class that allows to load a Skeleton_blocker_complex from an off file.
 */
template<typename Complex>
class Skeleton_blocker_off_reader {
 public:
  /**
   * name_file : file to read
   * read_complex : complex that will receive the file content
   * read_only_points : specify true if only the points must be read
   */
  Skeleton_blocker_off_reader(const std::string & name_file, Complex& read_complex,
                              bool read_only_points = false, bool is_flag = false) : valid_(false) {
    std::ifstream stream(name_file);
    if (stream.is_open()) {
      if (is_flag || read_only_points) {
        Skeleton_blocker_off_flag_visitor_reader<Complex> off_visitor(read_complex, read_only_points);
        Off_reader off_reader(stream);
        valid_ = off_reader.read(off_visitor);
      } else {
        Skeleton_blocker_off_visitor_reader<Complex> off_visitor(read_complex, read_only_points);
        Off_reader off_reader(stream);
        valid_ = off_reader.read(off_visitor);
      }
    }
  }

  /**
   * return true iff reading did not meet problems.
   */
  bool is_valid() const {
    return valid_;
  }

 private:
  bool valid_;
};

template<typename Complex>
class Skeleton_blocker_off_writer {
 public:
  /**
   * name_file : file where the off will be written
   * save_complex : complex that be outputted in the file
   * for now only save triangles.
   */
  Skeleton_blocker_off_writer(const std::string & name_file, const Complex& save_complex) {
    typedef typename Complex::Vertex_handle Vertex_handle;

    std::ofstream stream(name_file);
    if (stream.is_open()) {
      stream << "OFF\n";
      size_t num_triangles = std::distance(save_complex.triangle_range().begin(), save_complex.triangle_range().end());
      stream << save_complex.num_vertices() << " " << num_triangles << " 0 \n";

      // in case the complex has deactivated some vertices, eg only has vertices 0 2 5 7 for instance
      // we compute a map from 0 2 5 7 to 0 1 2 3
      std::map<Vertex_handle, size_t> vertex_num;
      size_t current_vertex = 0;

      for (auto v : save_complex.vertex_range()) {
        vertex_num[v] = current_vertex++;
        const auto& pt(save_complex.point(v));
        for (auto x : pt)
          stream << x << " ";
        stream << std::endl;
      }

      for (const auto & t : save_complex.triangle_range()) {
        stream << "3 ";
        for (auto x : t)
          stream << vertex_num[x] << " ";
        stream << std::endl;
      }
      stream.close();
    } else {
      std::cerr << "could not open file " << name_file << std::endl;
    }
  }
};

}  // namespace skbl

}  // namespace Gudhi

#endif  // SKELETON_BLOCKER_SKELETON_BLOCKER_OFF_IO_H_
