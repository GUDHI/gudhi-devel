/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2015 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */
#ifndef POINTS_3D_OFF_IO_H_
#define POINTS_3D_OFF_IO_H_

#include <gudhi/Off_reader.h>

#include <string>
#include <vector>
#include <fstream>
#include <map>

namespace Gudhi {

/** 
 * @brief OFF file visitor implementation according to Off_reader in order to read points from an OFF file.
 */
template<typename Point_3>
class Points_3D_off_visitor_reader {
 private:
  std::vector<Point_3> point_cloud_;
  bool valid_;

 public:
  /** @brief Off_reader visitor init implementation. 
   * 
   * The init parameters are set from OFF file header.
   * Dimension value is required and the value must be 3.
   *
   * @param[in] dim space dimension of vertices.
   * @param[in] num_vertices number of vertices in the OFF file (not used).
   * @param[in] num_faces number of faces in the OFF file (not used).
   * @param[in] num_edges number of edges in the OFF file (not used).
   */
  void init(int dim, int num_vertices, int num_faces, int num_edges) {
#ifdef DEBUG_TRACES
    std::cout << "Points_3D_off_visitor_reader::init - dim=" << dim << " - num_vertices=" <<
        num_vertices << " - num_faces=" << num_faces << " - num_edges=" << num_edges << std::endl;
#endif  // DEBUG_TRACES
    if (dim == 3) {
      valid_ = true;
    } else {
      valid_ = false;
      std::cerr << "Points_3D_off_reader::Points_3D_off_reader cannot read OFF files in dimension " << dim << "\n";
    }

    if (num_faces > 0) {
      std::cerr << "Points_3D_off_visitor_reader::init faces are not taken into account from OFF file for Points.\n";
    }
    if (num_edges > 0) {
      std::cerr << "Points_3D_off_visitor_reader::init edges are not taken into account from OFF file for Points.\n";
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
   * Point_3 must have a constructor with the following form:
   * 
   * @code template<class InputIterator > Point_3::Point_3(double x, double y, double z) @endcode
   */
  void point(const std::vector<double>& point) {
    if (valid_) {
#ifdef DEBUG_TRACES
      std::cout << "Points_3D_off_visitor_reader::point ";
      for (auto coordinate : point) {
        std::cout << coordinate << " | ";
      }
      std::cout << std::endl;
#endif  // DEBUG_TRACES
      // Fill the point cloud
      point_cloud_.push_back(Point_3(point[0], point[1], point[2]));
    }
  }

  // Off_reader visitor maximal_face implementation - Only points are read

  void maximal_face(const std::vector<int>& face) { }

  // Off_reader visitor done implementation - Only points are read

  void done() { }

  /** @brief Point cloud getter.
   *
   * @return The point cloud.
   */
  const std::vector<Point_3>& get_point_cloud() const {
    return point_cloud_;
  }

  /** @brief Returns if the OFF file read operation was successful or not.
   *
   * @return OFF file read status.
   */
  bool is_valid() const {
    return valid_;
  }
};

/** 
 * \@brief OFF file reader implementation in order to read dimension 3 points from an OFF file.
 * 
 * @details
 * This class is using the Points_3D_off_visitor_reader to visit the OFF file according to Off_reader.
 * 
 * Point_3 must have a constructor with the following form:
 * 
 * @code template<class InputIterator > Point_3::Point_3(double x, double y, double z) @endcode
 * 
 * @section point3doffioexample Example
 *
 * This example loads points from an OFF file and builds a vector of CGAL points in dimension 3.
 * Then, it is asked to display the points.
 * 
 * @include common/example_CGAL_3D_points_off_reader.cpp
 * 
 * When launching:
 * 
 * @code $> ./cgal3Doffreader ../../data/points/tore3D_300.off
 * @endcode
 *
 * the program output is:
 * 
 * @include common/cgal3Doffreader_result.txt
 */
template<typename Point_3>
class Points_3D_off_reader {
 public:
  /** @brief Reads the OFF file and constructs a vector of points from the points
   * that are in the OFF file.
   *
   * @param[in] name_file OFF file to read.
   * 
   * @post Check with is_valid() function to see if read operation was successful.
   */
  Points_3D_off_reader(const std::string& name_file)
      : valid_(false) {
    std::ifstream stream(name_file);
    if (stream.is_open()) {
      Off_reader off_reader(stream);
      Points_3D_off_visitor_reader<Point_3> off_visitor;
      valid_ = off_reader.read(off_visitor);
      valid_ = valid_ && off_visitor.is_valid();
      if (valid_) {
        point_cloud = off_visitor.get_point_cloud();
      }
    } else {
      std::cerr << "Points_3D_off_reader::Points_3D_off_reader could not open file " << name_file << "\n";
    }
  }

  /** @brief Returns if the OFF file read operation was successful or not.
   *
   * @return OFF file read status.
   */
  bool is_valid() const {
    return valid_;
  }

  /** @brief Point cloud getter.
   *
   * @return point_cloud.
   */
  const std::vector<Point_3>& get_point_cloud() const {
    return point_cloud;
  }

 private:
  /** @brief point_cloud.*/
  std::vector<Point_3> point_cloud;
  /** @brief OFF file read status.*/
  bool valid_;
};

}  // namespace Gudhi

#endif  // POINTS_3D_OFF_IO_H_
