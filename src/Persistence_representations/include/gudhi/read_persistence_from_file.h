/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - 2025/06 Hannah Schreiber: Various small bug fixes (missing `inline`s, `DEBUG_TRACES`s etc.)
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef READ_PERSISTENCE_FROM_FILE_H_
#define READ_PERSISTENCE_FROM_FILE_H_

#ifdef DEBUG_TRACES
#include <iostream> // std::clog
#endif
#include <utility>  // std::pair
#include <limits>   // for std::numeric_limits
#include <vector>
#include <string>

#include <gudhi/reader_utils.h>
#include <gudhi/Debug_utils.h>

namespace Gudhi {
namespace Persistence_representations {

/**
 * Universal procedure to read files with persistence. It ignores the lines starting from # (treat them as comments).
 * It reads the fist line which is not a comment and assume that there are some numerical entries over there. The
 * program assume that each other line in the file, which is not a comment, have the same number of numerical entries
 * (2, 3 or 4).
 * If there are two numerical entries per line, then the function assume that they are birth/death coordinates.
 * If there are three numerical entries per line, then the function assume that they are: dimension and birth/death
 * coordinates.
 * If there are four numerical entries per line, then the function assume that they are: the characteristic of a filed
 * over which persistence was computed, dimension and birth/death coordinates.
 * The 'inf' string can appear only as a last element of a line.
 * The procedure returns vector of persistence pairs.
 *
 * @ingroup Persistence_representations
 **/
inline std::vector<std::pair<double, double> > read_persistence_intervals_in_one_dimension_from_file(
    std::string const& filename,
    int dimension = -1,
    double what_to_substitute_for_infinite_bar = -1)
{
  using Bar = std::pair<double, double>;
  using Barcode = std::vector<Bar>;

  // TODO: `when what_to_substitute_for_infinite_bar` != -1, make changes directly in `barcode_initial` to avoid
  // allocating space? 
  // `final_barcode` makes sense when there are a lot of infinite bars which gets ignored (other
  // solution would be to swap those bars to the end and end with a resize(). The now unused space won't be freed,
  // but this do not differ from `final_barcode.reserve(barcode_initial.size())` used below).
  Barcode barcode_initial = read_persistence_intervals_in_dimension(filename, dimension);
  Barcode final_barcode;
  final_barcode.reserve(barcode_initial.size());

#ifdef DEBUG_TRACES
  std::clog << "Here are the intervals that we read from the file : \n";
  for (std::size_t i = 0; i != barcode_initial.size(); ++i) {
    std::clog << barcode_initial[i].first << " " << barcode_initial[i].second << std::endl;
  }
#endif

  for (std::size_t i = 0; i != barcode_initial.size(); ++i) {
#ifdef DEBUG_TRACES
    std::clog << "Considering interval : " << barcode_initial[i].first << " " << barcode_initial[i].second << std::endl;
#endif

    if (barcode_initial[i].first > barcode_initial[i].second) {
      // note that in this case barcode_initial[i].second != std::numeric_limits<double>::infinity()
#ifdef DEBUG_TRACES
      std::clog << "Swap and enter \n";
#endif
      // swap them to make sure that birth < death
      final_barcode.push_back(Bar(barcode_initial[i].second, barcode_initial[i].first));
    } else {
      if (barcode_initial[i].second != std::numeric_limits<double>::infinity()) {
#ifdef DEBUG_TRACES
        std::clog << "Simply enters\n";
#endif
        // in this case, due to the previous conditions we know that
        // barcode_initial[i].first < barcode_initial[i].second, so we put them as they are
        final_barcode.push_back(Bar(barcode_initial[i].first, barcode_initial[i].second));
      } else if (what_to_substitute_for_infinite_bar != -1){
        if (barcode_initial[i].first < what_to_substitute_for_infinite_bar) {
          // if only birth < death.
          final_barcode.push_back(Bar(barcode_initial[i].first, what_to_substitute_for_infinite_bar));
        }
      }
      // if the variable what_to_substitute_for_infinite_bar is not set, then we ignore all the infinite bars.
    }
  }

#ifdef DEBUG_TRACES
  std::clog << "Here are the final bars that we are sending further : \n";
  for (std::size_t i = 0; i != final_barcode.size(); ++i) {
    std::clog << final_barcode[i].first << " " << final_barcode[i].second << std::endl;
  }
  std::clog << "final_barcode.size() : " << final_barcode.size() << std::endl;
#endif

  return final_barcode;
}  // read_persistence_intervals_in_one_dimension_from_file

}  // namespace Persistence_representations
}  // namespace Gudhi

#endif  // READ_PERSISTENCE_FROM_FILE_H_
