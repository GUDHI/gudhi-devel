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
#include <string>

#include <gudhi/Zigzag_persistence.h>
#include <gudhi/Simplex_tree.h>

using ST = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_wide_indexation>;
using ZP = Gudhi::zigzag_persistence::Zigzag_persistence<ST>;
using Vertex_handle = ST::Vertex_handle;
using Filtration_value = ST::Filtration_value;
using interval_filtration = ZP::filtration_value_interval;

enum lineType : int { INCLUSION, REMOVAL, COMMENT };

void print_barcode(ZP& zp) {
  std::clog << std::endl << "Current barcode:" << std::endl;
  for (auto& bar : zp.get_persistence_diagram(0, true)) {
    std::clog << std::floor(bar.birth()) << " - ";
    if (bar.death() == std::numeric_limits<Filtration_value>::infinity()) {
      std::clog << "inf";
    } else {
      std::clog << std::floor(bar.death());
    }
    std::clog << " (" << bar.dim() << ")" << std::endl;
  }
  std::clog << std::endl;
}

lineType read_operation(std::string& line, std::vector<Vertex_handle>& vertices, double& timestamp) {
  lineType type;
  vertices.clear();
  Vertex_handle num;

  size_t current = line.find_first_not_of(' ', 0);
  if (current == std::string::npos) return COMMENT;

  if (line[current] == 'i')
    type = INCLUSION;
  else if (line[current] == 'r')
    type = REMOVAL;
  else if (line[current] == '#')
    return COMMENT;
  else {
    std::clog << "Syntaxe error in file." << std::endl;
    exit(0);
  }

  current = line.find_first_not_of(' ', current + 1);
  if (current == std::string::npos) {
    std::clog << "Syntaxe error in file." << std::endl;
    exit(0);
  }
  size_t next = line.find_first_of(' ', current);
  timestamp = std::stod(line.substr(current, next - current));

  current = line.find_first_not_of(' ', next);
  if (current == std::string::npos) {
    std::clog << "Syntaxe error in file." << std::endl;
    exit(0);
  }

  do {
    next = line.find_first_of(' ', current);
    num = std::stoi(line.substr(current, next - current));
    vertices.push_back(num);
    current = line.find_first_not_of(' ', next);
  } while (current != std::string::npos);

  return type;
}

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
  ZP zp;

  if (file.is_open()) {
    std::vector<Vertex_handle> vertices;
    double timestamp;
    lineType type;

    while (getline(file, line, '\n') && read_operation(line, vertices, timestamp) == COMMENT);
    double lastTimestamp = timestamp;
    // first operation has to be an insertion.
    zp.insert_simplex(vertices, timestamp);
    std::cout << line << std::endl;

    while (getline(file, line, '\n')) {
      type = read_operation(line, vertices, timestamp);
      if (type != COMMENT && lastTimestamp != timestamp) {
        print_barcode(zp);
        lastTimestamp = timestamp;
      }
      if (type != COMMENT) std::cout << line << std::endl;

      if (type == INCLUSION) {
        zp.insert_simplex(vertices, timestamp);
      } else if (type == REMOVAL) {
        zp.remove_simplex(vertices, timestamp);
      }
    }
    print_barcode(zp);

    file.close();
  } else {
    std::clog << "Unable to open input file." << std::endl;
    file.setstate(std::ios::failbit);
  }

  return 0;
}