/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */


#ifndef OFF_READER_H_
#define OFF_READER_H_


#include <sstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include <fstream>

namespace Gudhi {

/** \brief OFF file reader top class visitor. 
 * 
 * OFF file must be conform to \ref FileFormatsOFF
 */
class Off_reader {
 public:
  Off_reader(std::ifstream& stream) : stream_(stream) { }

  ~Off_reader() {
    stream_.close();
  }

  /** \brief
   * Read an OFF file and calls the following methods :
   * 
   * <CODE>void init(int dim,int num_vertices,int num_faces,int num_edges);  // from file header - num_edges may not be set
   * 
   * void point(const std::vector<double>& point);  // for each point read
   * 
   * void maximal_face(const std::list<int>& face);  // for each face read
   * 
   * void done();  // upon file read is finished</CODE>
   * 
   * of the visitor when reading a point or a maximal face. Edges are not taken into account.
   */
  template<typename OffVisitor>
  bool read(OffVisitor& off_visitor) {
    bool success_read_off_preamble = read_off_preamble(off_visitor);
    if (!success_read_off_preamble) {
      std::cerr << "could not read off preambule\n";
      return false;
    }

    bool success_read_off_points = read_off_points(off_visitor);
    if (!success_read_off_points) {
      std::cerr << "could not read off points\n";
      return false;
    }

    bool success_read_off_faces = read_off_faces(off_visitor);
    if (!success_read_off_faces) {
      std::cerr << "could not read off faces\n";
      return false;
    }

    off_visitor.done();
    return success_read_off_preamble && success_read_off_points && success_read_off_faces;
  }

 private:
  std::ifstream& stream_;

  struct Off_info {
    int dim;
    int num_vertices;
    int num_edges;
    int num_faces;
  };

  Off_info off_info_;

  template<typename OffVisitor>
  bool read_off_preamble(OffVisitor& off_visitor) {
    std::string line;
    if (!goto_next_uncomment_line(line)) return false;

    bool is_off_file = (line.find("OFF") != std::string::npos);
    bool is_noff_file = (line.find("nOFF") != std::string::npos);



    if (!is_off_file && !is_noff_file) {
      std::cerr << line << std::endl;
      std::cerr << "missing off header\n";
      return false;
    }

    if (is_noff_file) {
      // Should be on a separate line, but we accept it on the same line as the number of vertices
      stream_ >> off_info_.dim;
    } else {
      off_info_.dim = 3;
    }

    if (!goto_next_uncomment_line(line)) return false;
    std::istringstream iss(line);
    if (!(iss >> off_info_.num_vertices >> off_info_.num_faces >> off_info_.num_edges)) {
      std::cerr << "incorrect number of vertices/faces/edges\n";
      return false;
    }
    off_visitor.init(off_info_.dim, off_info_.num_vertices, off_info_.num_faces, off_info_.num_edges);

    return true;
  }

  bool goto_next_uncomment_line(std::string& uncomment_line) {
    do {
      // skip whitespace, including empty lines
      if (!std::ifstream::sentry(stream_)) return false;
      std::getline(stream_, uncomment_line);
    } while (uncomment_line[0] == '#');
    return static_cast<bool>(stream_);
  }

  template<typename OffVisitor>
  bool read_off_points(OffVisitor& visitor) {
    int num_vertices_to_read = off_info_.num_vertices;
    while (num_vertices_to_read--) {
      std::string line;
      if (!goto_next_uncomment_line(line)) return false;
      std::vector<double> point;
      std::istringstream iss(line);
      point.assign(std::istream_iterator<double>(iss), std::istream_iterator<double>());
      // if(point.size() != off_info_.dim) return false;
      visitor.point(point);
    }
    return true;
  }

  template<typename OffVisitor>
  bool read_off_faces(OffVisitor& visitor) {
    std::string line;
    while (goto_next_uncomment_line(line)) {
      std::istringstream iss(line);
      int num_face_vertices;
      iss >> num_face_vertices;
      std::vector<int> face;
      face.assign(std::istream_iterator<int>(iss), std::istream_iterator<int>());
      // if (face.size() != (off_info_.dim + 1)) return false;
      visitor.maximal_face(face);
    }
    return true;
  }
};

template<typename OFFVisitor>
void read_off(const std::string& name_file_off, OFFVisitor& vis) {
  std::ifstream stream(name_file_off);
  if (!stream.is_open()) {
    std::cerr << "could not open file \n";
  } else {
    Off_reader off_reader(stream);
    off_reader.read(vis);
  }
}

}  // namespace Gudhi

#endif  // OFF_READER_H_
