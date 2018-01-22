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
 *
 */


#ifndef OFF_POINT_ITERATOR_H_
#define OFF_POINT_ITERATOR_H_


#include <sstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include <fstream>

namespace Gudhi {

bool goto_next_uncomment_line(std::string& uncomment_line, std::istream& stream) {
  uncomment_line.clear();
  do
    std::getline(stream, uncomment_line); while (uncomment_line[0] == '%');
  return (uncomment_line.size() > 0 && uncomment_line[0] != '%');
}
  
/** \brief OFF file reader range. 
 * 
 * OFF file must be conform to format described here : 
 * http://www.geomview.org/docs/html/OFF.html
 */
template <class Point_d>
class Off_point_range {
public:

  class iterator_ : public boost::iterator_facade< iterator_,
                                                   Point_d const,
                                                   boost::forward_traversal_tag> {
    private:
    std::ifstream& stream_;
    Point_d value_;
    unsigned vertices_left_;
    
    friend class boost::iterator_core_access;

    bool update_value() {
      if (vertices_left_) {
        std::string line;
        if (!goto_next_uncomment_line(line, stream_)) return false;
        std::vector<double> point;
        std::istringstream iss(line);
        point.assign(std::istream_iterator<double>(iss), std::istream_iterator<double>());
        value_ = Point_d(point);
        vertices_left_--;
      }
      return true;
    }
    
    bool equal(iterator_ const& other) const {
      return vertices_left_ == other.vertices_left_;
    }

    Point_d const& dereference() const {
      return value_;
    }

    void increment() {
      update_value();
    }
    
  public:
    iterator_(std::ifstream& stream,
              unsigned vertices_left)
        : stream_(stream), vertices_left_(vertices_left) {
        update_value();
    }

    void update_number_of_vertices(unsigned vertices_left) {
      vertices_left_ = vertices_left;
      update_value();
    }

    int dimension() const {
      return value_.size();
    }
    
  };

  typedef iterator_ iterator;
  typedef iterator_ const_iterator;  
  
  Off_point_range(const std::string& name_file_off)
    : begin_(stream_,0), end_(stream_,0) {
    stream_ = std::ifstream(name_file_off);
    if (!stream_.is_open()) {
      std::cerr << "could not open file \n";
    }
    unsigned number_of_vertices = 0;
    bool success_read_off_preambule = read_off_preambule(number_of_vertices);
    if (!success_read_off_preambule)
      std::cerr << "could not read off preambule\n";
    begin_.update_number_of_vertices(number_of_vertices);
    dim_ = begin_.dimension();
  }

  ~Off_point_range() {
    stream_.close();
  }

  iterator begin() const {
    return begin_;
  }

  iterator end() const {
    return end_;
  }

  int dimension() {
    return dim_;
  }

private:
  bool read_off_preambule(unsigned& number_of_vertices) {
    std::string line;
    if (!goto_next_uncomment_line(line, stream_)) return false;

    bool is_off_file = (line.find("OFF") != std::string::npos);
    bool is_noff_file = (line.find("nOFF") != std::string::npos);

    if (!is_off_file && !is_noff_file) {
      std::cerr << line << std::endl;
      std::cerr << "missing off header\n";
      return false;
    }

    if (!goto_next_uncomment_line(line, stream_)) return false;
    std::istringstream iss(line);
    int num_faces, num_edges, dim;
    if ((is_off_file) && (!is_noff_file)) {
      if (!(iss >> number_of_vertices >> num_faces >> num_edges)) {
        std::cerr << "incorrect number of vertices/faces/edges\n";
        return false;
      }
    } else {
      if (!(iss >> dim >> number_of_vertices >> num_faces >> num_edges)) {
      std::cerr << "incorrect number of vertices/faces/edges\n";
      return false;
      }
    }
    return true;
  }

 private:
  std::ifstream stream_;
  iterator begin_, end_;
  int dim_;
};

}  // namespace Gudhi

#endif  // OFF_READER_H_
