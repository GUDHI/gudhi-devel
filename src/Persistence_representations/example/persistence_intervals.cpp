/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Persistence_intervals.h>

#include <iostream>
#include <utility>
#include <vector>

using Persistence_intervals = Gudhi::Persistence_representations::Persistence_intervals;

int main(int argc, char** argv) {
  if (argc != 2) {
    std::clog << "To run this program, please provide the name of a file with persistence diagram \n";
    return 1;
  }

  Persistence_intervals p(argv[1]);
  std::pair<double, double> min_max_ = p.get_x_range();
  std::clog << "Birth-death range : " << min_max_.first << " " << min_max_.second << std::endl;

  std::vector<double> dominant_ten_intervals_length = p.length_of_dominant_intervals(10);
  std::clog << "Length of ten dominant intervals : " << std::endl;
  for (size_t i = 0; i != dominant_ten_intervals_length.size(); ++i) {
    std::clog << dominant_ten_intervals_length[i] << std::endl;
  }

  std::vector<std::pair<double, double> > ten_dominant_intervals = p.dominant_intervals(10);
  std::clog << "Here are the dominant intervals : " << std::endl;
  for (size_t i = 0; i != ten_dominant_intervals.size(); ++i) {
    std::clog << "( " << ten_dominant_intervals[i].first << "," << ten_dominant_intervals[i].second << std::endl;
  }

  std::vector<size_t> histogram = p.histogram_of_lengths(10);
  std::clog << "Here is the histogram of barcode's length : " << std::endl;
  for (size_t i = 0; i != histogram.size(); ++i) {
    std::clog << histogram[i] << " ";
  }
  std::clog << std::endl;

  std::vector<size_t> cumulative_histogram = p.cumulative_histogram_of_lengths(10);
  std::clog << "Cumulative histogram : " << std::endl;
  for (size_t i = 0; i != cumulative_histogram.size(); ++i) {
    std::clog << cumulative_histogram[i] << " ";
  }
  std::clog << std::endl;

  std::vector<double> char_funct_diag = p.characteristic_function_of_diagram(min_max_.first, min_max_.second);
  std::clog << "Characteristic function of diagram : " << std::endl;
  for (size_t i = 0; i != char_funct_diag.size(); ++i) {
    std::clog << char_funct_diag[i] << " ";
  }
  std::clog << std::endl;

  std::vector<double> cumul_char_funct_diag =
      p.cumulative_characteristic_function_of_diagram(min_max_.first, min_max_.second);
  std::clog << "Cumulative characteristic function of diagram : " << std::endl;
  for (size_t i = 0; i != cumul_char_funct_diag.size(); ++i) {
    std::clog << cumul_char_funct_diag[i] << " ";
  }
  std::clog << std::endl;

  std::clog << "Persistence Betti numbers \n";
  std::vector<std::pair<double, size_t> > pbns = p.compute_persistent_betti_numbers();
  for (size_t i = 0; i != pbns.size(); ++i) {
    std::clog << pbns[i].first << " " << pbns[i].second << std::endl;
  }

  return 0;
}
