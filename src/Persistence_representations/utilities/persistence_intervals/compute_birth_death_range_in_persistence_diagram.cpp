/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Persistence_intervals.h>

#include <iostream>
#include <vector>
#include <limits>
#include <utility>

using Persistence_intervals = Gudhi::Persistence_representations::Persistence_intervals;

int main(int argc, char** argv) {
  std::cout << "This program computes the range of birth and death times of persistence pairs in diagrams provided as "
            << "an input.\n"
            << "The first parameter is the dimension of persistence to be used to create persistence intervals. "
            << "If your file contains the information about dimension of persistence pairs, please provide here the "
            << "dimension of persistence pairs you want to use. "
            << "If your input files consist only of birth-death pairs, please set this first parameter to -1.\n"
            << "The remaining parameters of the program are the names of files with persistence diagrams.\n";

  if (argc < 3) {
    std::cout << "Wrong parameter list, the program will now terminate \n";
    return 1;
  }

  unsigned dimension = std::numeric_limits<unsigned>::max();
  int dim = atoi(argv[1]);
  if (dim >= 0) {
    dimension = (unsigned)dim;
  }
  std::vector<const char*> filenames;
  for (int i = 2; i < argc; ++i) {
    filenames.push_back(argv[i]);
  }

  double min_ = std::numeric_limits<double>::max();
  double max_ = -std::numeric_limits<double>::max();

  for (size_t file_no = 0; file_no != filenames.size(); ++file_no) {
    std::cout << "Creating diagram based on a file : " << filenames[file_no] << std::endl;
    Persistence_intervals p(filenames[file_no], dimension);
    std::pair<double, double> min_max_ = p.get_x_range();
    if (min_max_.first < min_) min_ = min_max_.first;
    if (min_max_.second > max_) max_ = min_max_.second;
  }
  std::cout << "Birth-death range : min: " << min_ << ", max: " << max_ << std::endl;
  return 0;
}
