/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <iostream>
#include <ostream>
#include <fstream>
#include <string>

#include <gudhi/filtered_zigzag_persistence.h>

using ZP = Gudhi::zigzag_persistence::Filtered_zigzag_persistence<>;
using id_handle = ZP::face_key;
using filtration_value = ZP::filtration_value;
using dimension_type = ZP::dimension_type;

enum lineType : int { INCLUSION, REMOVAL, COMMENT };

lineType read_operation(std::string& line, std::vector<id_handle>& faces, double& timestamp) {
  lineType type;
  faces.clear();
  id_handle num;

  size_t current = line.find_first_not_of(' ', 0);
  if (current == std::string::npos) return COMMENT;

  if (line[current] == 'i')
    type = INCLUSION;
  else if (line[current] == 'r')
    type = REMOVAL;
  else if (line[current] == '#')
    return COMMENT;
  else {
    std::clog << "(1) Syntaxe error in file." << std::endl;
    exit(0);
  }

  current = line.find_first_not_of(' ', current + 1);
  if (current == std::string::npos) {
    std::clog << "(2) Syntaxe error in file." << std::endl;
    exit(0);
  }
  size_t next = line.find_first_of(' ', current);
  timestamp = std::stod(line.substr(current, next - current));

  current = line.find_first_not_of(' ', next);
  while (current != std::string::npos) {
    next = line.find_first_of(' ', current);
    num = std::stoi(line.substr(current, next - current));
    faces.push_back(num);
    current = line.find_first_not_of(' ', next);
  }

  return type;
}

//example of input file: example/zigzag_filtration_example.txt
int main(int argc, char* const argv[]) {
  if (argc != 2) {
    if (argc < 2)
      std::clog << "Missing argument: input file name is needed." << std::endl;
    else
      std::clog << "Too many arguments: only input file name is needed." << std::endl;
    return 0;
  }

  std::string line;
  std::ifstream file(argv[1]);

  //std::cout could be replaced by any other output stream
  ZP zp([](dimension_type dim, filtration_value birth, filtration_value death) {
    std::cout << "[" << dim << "] ";
    std::cout << birth << " - " << death;
    std::cout << std::endl;
  });

  if (file.is_open()) {
    std::vector<id_handle> data;
    unsigned int id = 0;
    double timestamp;
    lineType type;

    while (getline(file, line, '\n') && read_operation(line, data, timestamp) == COMMENT);
    double lastTimestamp = timestamp;
    // first operation has to be an insertion.
    zp.insert_face(id, data, 0, timestamp);

    while (getline(file, line, '\n')) {
      type = read_operation(line, data, timestamp);
      if (type != COMMENT && lastTimestamp != timestamp) {
        lastTimestamp = timestamp;
      }

      if (type == INCLUSION) {
        ++id;
        int dim = data.size() == 0 ? 0 : data.size() - 1;
        zp.insert_face(id, data, dim, timestamp);
      } else if (type == REMOVAL) {
        ++id;
        zp.remove_face(data[0], data[1], timestamp);
      }
    }

    file.close();
  } else {
    std::clog << "Unable to open input file." << std::endl;
    file.setstate(std::ios::failbit);
  }

  //retrieve infinit bars remaining at the end
  //again std::cout could be replaced by any other output stream
  zp.get_current_infinite_intervals([](dimension_type dim, filtration_value birth) {
    std::cout << "[" << dim << "] ";
    std::cout << birth << " - inf";
    std::cout << std::endl;
  });

  return 0;
}