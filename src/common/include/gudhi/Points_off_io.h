/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2015 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */
#ifndef POINTS_OFF_IO_H_
#define POINTS_OFF_IO_H_

#include <gudhi/Off_reader.h>

#include <string>
#include <vector>
#include <fstream>
#include <map>

namespace Gudhi {

/** 
 * \brief OFF file visitor implementation according to Off_reader in order to read points from an OFF file.
 */
template<typename Point_d>
class Points_off_visitor_reader {
 private:
  std::vector<Point_d> point_cloud;

 public:
  /** \brief Off_reader visitor init implementation. 
   * 
   * The init parameters are set from OFF file header.
   * Dimension value is required in order to construct a vector of points.
   *
   * @param[in] dim space dimension of vertices.
   * @param[in] num_vertices number of vertices in the OFF file (not used).
   * @param[in] num_faces number of faces in the OFF file (not used).
   * @param[in] num_edges number of edges in the OFF file (not used).
   */
  void init(int dim, int num_vertices, int num_faces, int num_edges) {
#ifdef DEBUG_TRACES
    std::clog << "Points_off_visitor_reader::init - dim=" << dim << " - num_vertices=" <<
        num_vertices << " - num_faces=" << num_faces << " - num_edges=" << num_edges << std::endl;
#endif  // DEBUG_TRACES
    if (num_faces > 0) {
      std::cerr << "Points_off_visitor_reader::init faces are not taken into account from OFF file for Points.\n";
    }
    if (num_edges > 0) {
      std::cerr << "Points_off_visitor_reader::init edges are not taken into account from OFF file for Points.\n";
    }
  }

  /** @brief Off_reader visitor point implementation. 
   * 
   * The point function is called on each vertex line from OFF file.
   * This function inserts the vertex in the vector of points.
   *
   * @param[in] point vector of vertex coordinates.
   * 
   * @details
   * Point_d must have a constructor with the following form:
   * 
   * @code template<class InputIterator > Point_d::Point_d(InputIterator first, InputIterator last) @endcode
   * 
   */
  void point(const std::vector<double>& point) {
#ifdef DEBUG_TRACES
    std::clog << "Points_off_visitor_reader::point ";
    for (auto coordinate : point) {
      std::clog << coordinate << " | ";
    }
    std::clog << std::endl;
#endif  // DEBUG_TRACES
    // Fill the point cloud
    point_cloud.push_back(Point_d(point.begin(), point.end()));
  }

  // Off_reader visitor maximal_face implementation - Only points are read
  void maximal_face(const std::vector<int>& face) { }

  // Off_reader visitor done implementation - Only points are read
  void done() { }

  /** \brief Point cloud getter.
   *
   * @return point_cloud.
   */
  const std::vector<Point_d>& get_point_cloud() const {
    return point_cloud;
  }
};

/** 
 * \brief OFF file reader implementation in order to read points from an OFF file.
 * 
 * This class is using the Points_off_visitor_reader to visit the OFF file according to Off_reader.
 * 
 * Point_d must have a constructor with the following form:
 * 
 * \code template<class InputIterator > Point_d::Point_d(int d, InputIterator first, InputIterator last) \endcode
 * 
 * where d is the point dimension. 		
 * 
 * \section pointoffioexample Example
 *
 * This example loads points from an OFF file and builds a vector of points (vector of double).
 * Then, it is asked to display the points.
 * 
 * \include example_vector_double_points_off_reader.cpp
 * 
 * When launching:
 * 
 * \code $> ./vector_double_off_reader ../../data/points/alphacomplexdoc.off
 * \endcode
 *
 * the program outputs a file ../../data/points/alphacomplexdoc.off.txt:
 * 
 * \include vectordoubleoffreader_result.txt
 */
template<typename Point_d>
class Points_off_reader {
 public:
  /** \brief Reads the OFF file and constructs a vector of points from the points
   * that are in the OFF file.
   *
   * @param[in] name_file OFF file to read.
   * 
   * \post Check with is_valid() function to see if read operation was successful.
   */
  Points_off_reader(const std::string& name_file)
  : valid_(false) {
    std::ifstream stream(name_file);
    if (stream.is_open()) {
      Off_reader off_reader(stream);
      Points_off_visitor_reader<Point_d> off_visitor;
      valid_ = off_reader.read(off_visitor);
      if (valid_) {
        point_cloud = off_visitor.get_point_cloud();
      }
    } else {
      std::cerr << "Points_off_reader::Points_off_reader could not open file " << name_file << "\n";
    }
  }

  /** \brief Returns if the OFF file read operation was successful or not.
   *
   * @return OFF file read status.
   */
  bool is_valid() const {
    return valid_;
  }

   /** \brief Point cloud getter.
   *
   * @return point_cloud.
   */
  const std::vector<Point_d>& get_point_cloud() const {
    return point_cloud;
  }

 private:
  /** \brief point_cloud.*/
  std::vector<Point_d> point_cloud;
  /** \brief OFF file read status.*/
  bool valid_;
};

}  // namespace Gudhi

#endif  // POINTS_OFF_IO_H_
