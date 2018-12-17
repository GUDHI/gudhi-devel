/*    This file is part of the Gudhi hiLibrary. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016 Inria
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

#ifndef PERSISTENCE_INTERVALS_H_
#define PERSISTENCE_INTERVALS_H_

// gudhi include
#include <gudhi/read_persistence_from_file.h>

// standard include
#include <limits>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <functional>
#include <utility>
#include <string>

namespace Gudhi {
namespace Persistence_representations {

/**
 * This class implements the following concepts: Vectorized_topological_data, Topological_data_with_distances,
 *Real_valued_topological_data
**/
class Persistence_intervals {
 public:
  /**
   * This is a constructor of a class Persistence_intervals from a text file. Each line of the input file is supposed to
   *contain two numbers of a type double (or convertible to double)
   * representing the birth and the death of the persistence interval. If the pairs are not sorted so that birth <=
   *death, then the constructor will sort then that way.
   * * The second parameter of a constructor is a dimension of intervals to be read from a file. If your file contains
   *only birth-death pairs, use the default value.
  **/
  Persistence_intervals(const char* filename, unsigned dimension = std::numeric_limits<unsigned>::max());

  /**
   * This is a constructor of a class Persistence_intervals from a vector of pairs. Each pair is assumed to represent a
   *persistence interval. We assume that the first elements of pairs
   * are smaller or equal the second elements of pairs.
  **/
  Persistence_intervals(const std::vector<std::pair<double, double> >& intervals);

  /**
       * This procedure returns x-range of a given persistence diagram.
      **/
  std::pair<double, double> get_x_range() const {
    double min_ = std::numeric_limits<int>::max();
    double max_ = -std::numeric_limits<int>::max();
    for (size_t i = 0; i != this->intervals.size(); ++i) {
      if (this->intervals[i].first < min_) min_ = this->intervals[i].first;
      if (this->intervals[i].second > max_) max_ = this->intervals[i].second;
    }
    return std::make_pair(min_, max_);
  }

  /**
   * This procedure returns y-range of a given persistence diagram.
  **/
  std::pair<double, double> get_y_range() const {
    double min_ = std::numeric_limits<int>::max();
    double max_ = -std::numeric_limits<int>::max();
    for (size_t i = 0; i != this->intervals.size(); ++i) {
      if (this->intervals[i].second < min_) min_ = this->intervals[i].second;
      if (this->intervals[i].second > max_) max_ = this->intervals[i].second;
    }
    return std::make_pair(min_, max_);
  }

  /**
   * Procedure that compute the vector of lengths of the dominant (i.e. the longest) persistence intervals. The list is
   *truncated at the parameter of the call where_to_cut (set by default to 100).
  **/
  std::vector<double> length_of_dominant_intervals(size_t where_to_cut = 100) const;

  /**
       * Procedure that compute the vector of the dominant (i.e. the longest) persistence intervals. The parameter of
    *the procedure (set by default to 100) is the number of dominant intervals returned by the procedure.
      **/
  std::vector<std::pair<double, double> > dominant_intervals(size_t where_to_cut = 100) const;

  /**
       * Procedure to compute a histogram of interval's length. A histogram is a block plot. The number of blocks is
    *determined by the first parameter of the function (set by default to 10).
       * For the sake of argument let us assume that the length of the longest interval is 1 and the number of bins is
    *10. In this case the i-th block correspond to a range between i-1/10 and i10.
       * The vale of a block supported at the interval is the number of persistence intervals of a length between x_0
    *and x_1.
      **/
  std::vector<size_t> histogram_of_lengths(size_t number_of_bins = 10) const;

  /**
   * Based on a histogram of intervals lengths computed by the function histogram_of_lengths H the procedure below
   *computes the cumulative histogram. The i-th position of the resulting histogram
   * is the sum of values of H for the positions from 0 to i.
  **/
  std::vector<size_t> cumulative_histogram_of_lengths(size_t number_of_bins = 10) const;

  /**
   * In this procedure we assume that each barcode is a characteristic function of a hight equal to its length. The
   *persistence diagram is a sum of such a functions. The procedure below construct a function being a
   * sum of the characteristic functions of persistence intervals. The first two parameters are the range in which the
   *function is to be computed and the last parameter is the number of bins in
   * the discretization of the interval [_min,_max].
  **/
  std::vector<double> characteristic_function_of_diagram(double x_min, double x_max, size_t number_of_bins = 10) const;

  /**
   * Cumulative version of the function characteristic_function_of_diagram
  **/
  std::vector<double> cumulative_characteristic_function_of_diagram(double x_min, double x_max,
                                                                    size_t number_of_bins = 10) const;

  /**
   * Compute the function of persistence Betti numbers. The returned value is a vector of pair. First element of each
   *pair is a place where persistence Betti numbers change.
   * Second element of each pair is the value of Persistence Betti numbers at that point.
  **/
  std::vector<std::pair<double, size_t> > compute_persistent_betti_numbers() const;

  /**
   *This is a non optimal procedure that compute vector of distances from each point of diagram to its k-th nearest
   *neighbor (k is a parameter of the program). The resulting vector is by default truncated to 10
   *elements (this value can be changed by using the second parameter of the program). The points are returned in order
   *from the ones which are farthest away from their k-th nearest neighbors.
  **/
  std::vector<double> k_n_n(size_t k, size_t where_to_cut = 10) const;

  /**
* Operator that send the diagram to a stream.
**/
  friend std::ostream& operator<<(std::ostream& out, const Persistence_intervals& intervals) {
    for (size_t i = 0; i != intervals.intervals.size(); ++i) {
      out << intervals.intervals[i].first << " " << intervals.intervals[i].second << std::endl;
    }
    return out;
  }

  /**
   * Generating gnuplot script to plot the interval.
  **/
  void plot(const char* filename, double min_x = std::numeric_limits<double>::max(),
            double max_x = std::numeric_limits<double>::max(), double min_y = std::numeric_limits<double>::max(),
            double max_y = std::numeric_limits<double>::max()) const {
    // this program create a gnuplot script file that allows to plot persistence diagram.
    std::ofstream out;

    std::stringstream gnuplot_script;
    gnuplot_script << filename << "_GnuplotScript";

    out.open(gnuplot_script.str().c_str());

    std::pair<double, double> min_max_values = this->get_x_range();
    if (min_x == max_x) {
      out << "set xrange [" << min_max_values.first - 0.1 * (min_max_values.second - min_max_values.first) << " : "
          << min_max_values.second + 0.1 * (min_max_values.second - min_max_values.first) << " ]" << std::endl;
      out << "set yrange [" << min_max_values.first - 0.1 * (min_max_values.second - min_max_values.first) << " : "
          << min_max_values.second + 0.1 * (min_max_values.second - min_max_values.first) << " ]" << std::endl;
    } else {
      out << "set xrange [" << min_x << " : " << max_x << " ]" << std::endl;
      out << "set yrange [" << min_y << " : " << max_y << " ]" << std::endl;
    }
    out << "plot '-' using 1:2 notitle \"" << filename << "\", \\" << std::endl;
    out << "     '-' using 1:2 notitle with lp" << std::endl;
    for (size_t i = 0; i != this->intervals.size(); ++i) {
      out << this->intervals[i].first << " " << this->intervals[i].second << std::endl;
    }
    out << "EOF" << std::endl;
    out << min_max_values.first - 0.1 * (min_max_values.second - min_max_values.first) << " "
        << min_max_values.first - 0.1 * (min_max_values.second - min_max_values.first) << std::endl;
    out << min_max_values.second + 0.1 * (min_max_values.second - min_max_values.first) << " "
        << min_max_values.second + 0.1 * (min_max_values.second - min_max_values.first) << std::endl;

    out.close();

    std::cout << "To visualize, install gnuplot and type the command: gnuplot -persist -e \"load \'"
              << gnuplot_script.str().c_str() << "\'\"" << std::endl;
  }

  /**
* Return number of points in the diagram.
**/
  size_t size() const { return this->intervals.size(); }

  /**
   * Return the persistence interval at the given position. Note that intervals are not sorted with respect to their
   *lengths.
  **/
  inline std::pair<double, double> operator[](size_t i) const {
    if (i >= this->intervals.size()) throw("Index out of range! Operator [], one_d_gaussians class\n");
    return this->intervals[i];
  }

  // Implementations of functions for various concepts.
  /**
   * This is a simple function projecting the persistence intervals to a real number. The function we use here is a sum
   *of squared lengths of intervals. It can be naturally interpreted as
   * sum of step function, where the step hight it equal to the length of the interval.
   * At the moment this function is not tested, since it is quite likely to be changed in the future. Given this, when
   *using it, keep in mind that it
   * will be most likely changed in the next versions.
   **/
  double project_to_R(int number_of_function) const;
  /**
   * The function gives the number of possible projections to R. This function is required by the
   *Real_valued_topological_data concept.
  **/
  size_t number_of_projections_to_R() const { return this->number_of_functions_for_projections_to_reals; }

  /**
   * Return a family of vectors obtained from the persistence diagram. The i-th vector consist of the length of i
   *dominant persistence intervals.
  **/
  std::vector<double> vectorize(int number_of_function) const {
    return this->length_of_dominant_intervals(number_of_function);
  }
  /**
      * This function return the number of functions that allows vectorization of a persistence diagram. It is required
    *in a concept Vectorized_topological_data.
      **/
  size_t number_of_vectorize_functions() const { return this->number_of_functions_for_vectorization; }

  // end of implementation of functions needed for concepts.

  // For visualization use output from vectorize and build histograms.
  std::vector<std::pair<double, double> > output_for_visualization() { return this->intervals; }

 protected:
  void set_up_numbers_of_functions_for_vectorization_and_projections_to_reals() {
    // warning, this function can be only called after filling in the intervals vector.
    this->number_of_functions_for_vectorization = this->intervals.size();
    this->number_of_functions_for_projections_to_reals = 1;
  }

  std::vector<std::pair<double, double> > intervals;
  size_t number_of_functions_for_vectorization;
  size_t number_of_functions_for_projections_to_reals;
};

Persistence_intervals::Persistence_intervals(const char* filename, unsigned dimension) {
  if (dimension == std::numeric_limits<unsigned>::max()) {
    this->intervals = read_persistence_intervals_in_one_dimension_from_file(filename);
  } else {
    this->intervals = read_persistence_intervals_in_one_dimension_from_file(filename, dimension);
  }
  this->set_up_numbers_of_functions_for_vectorization_and_projections_to_reals();
}  // Persistence_intervals

Persistence_intervals::Persistence_intervals(const std::vector<std::pair<double, double> >& intervals_)
    : intervals(intervals_) {
  this->set_up_numbers_of_functions_for_vectorization_and_projections_to_reals();
}

std::vector<double> Persistence_intervals::length_of_dominant_intervals(size_t where_to_cut) const {
  std::vector<double> result(this->intervals.size());
  for (size_t i = 0; i != this->intervals.size(); ++i) {
    result[i] = this->intervals[i].second - this->intervals[i].first;
  }
  std::sort(result.begin(), result.end(), std::greater<double>());

  result.resize(std::min(where_to_cut, result.size()));
  return result;
}  // length_of_dominant_intervals

bool compare(const std::pair<size_t, double>& first, const std::pair<size_t, double>& second) {
  return first.second > second.second;
}

std::vector<std::pair<double, double> > Persistence_intervals::dominant_intervals(size_t where_to_cut) const {
  bool dbg = false;
  std::vector<std::pair<size_t, double> > position_length_vector(this->intervals.size());
  for (size_t i = 0; i != this->intervals.size(); ++i) {
    position_length_vector[i] = std::make_pair(i, this->intervals[i].second - this->intervals[i].first);
  }

  std::sort(position_length_vector.begin(), position_length_vector.end(), compare);

  std::vector<std::pair<double, double> > result;
  result.reserve(std::min(where_to_cut, position_length_vector.size()));

  for (size_t i = 0; i != std::min(where_to_cut, position_length_vector.size()); ++i) {
    result.push_back(this->intervals[position_length_vector[i].first]);
    if (dbg)
      std::cerr << "Position : " << position_length_vector[i].first << " length : " << position_length_vector[i].second
                << std::endl;
  }

  return result;
}  // dominant_intervals

std::vector<size_t> Persistence_intervals::histogram_of_lengths(size_t number_of_bins) const {
  bool dbg = false;

  if (dbg) std::cerr << "this->intervals.size() : " << this->intervals.size() << std::endl;
  // first find the length of the longest interval:
  double lengthOfLongest = 0;
  for (size_t i = 0; i != this->intervals.size(); ++i) {
    if ((this->intervals[i].second - this->intervals[i].first) > lengthOfLongest) {
      lengthOfLongest = this->intervals[i].second - this->intervals[i].first;
    }
  }

  if (dbg) {
    std::cerr << "lengthOfLongest : " << lengthOfLongest << std::endl;
  }

  // this is a container we will use to store the resulting histogram
  std::vector<size_t> result(number_of_bins + 1, 0);

  // for every persistence interval in our collection.
  for (size_t i = 0; i != this->intervals.size(); ++i) {
    // compute its length relative to the length of the dominant interval:
    double relative_length_of_this_interval = (this->intervals[i].second - this->intervals[i].first) / lengthOfLongest;

    // given the relative length (between 0 and 1) compute to which bin should it contribute.
    size_t position = (size_t)(relative_length_of_this_interval * number_of_bins);

    ++result[position];

    if (dbg) {
      std::cerr << "i : " << i << std::endl;
      std::cerr << "Interval : [" << this->intervals[i].first << " , " << this->intervals[i].second << " ] \n";
      std::cerr << "relative_length_of_this_interval : " << relative_length_of_this_interval << std::endl;
      std::cerr << "position : " << position << std::endl;
      getchar();
    }
  }

  if (dbg) {
    for (size_t i = 0; i != result.size(); ++i) std::cerr << result[i] << std::endl;
  }
  return result;
}

std::vector<size_t> Persistence_intervals::cumulative_histogram_of_lengths(size_t number_of_bins) const {
  std::vector<size_t> histogram = this->histogram_of_lengths(number_of_bins);
  std::vector<size_t> result(histogram.size());

  size_t sum = 0;
  for (size_t i = 0; i != histogram.size(); ++i) {
    sum += histogram[i];
    result[i] = sum;
  }
  return result;
}

std::vector<double> Persistence_intervals::characteristic_function_of_diagram(double x_min, double x_max,
                                                                              size_t number_of_bins) const {
  bool dbg = false;

  std::vector<double> result(number_of_bins);
  std::fill(result.begin(), result.end(), 0);

  for (size_t i = 0; i != this->intervals.size(); ++i) {
    if (dbg) {
      std::cerr << "Interval : " << this->intervals[i].first << " , " << this->intervals[i].second << std::endl;
    }

    size_t beginIt = 0;
    if (this->intervals[i].first < x_min) beginIt = 0;
    if (this->intervals[i].first >= x_max) beginIt = result.size();
    if ((this->intervals[i].first > x_min) && (this->intervals[i].first < x_max)) {
      beginIt = number_of_bins * (this->intervals[i].first - x_min) / (x_max - x_min);
    }

    size_t endIt = 0;
    if (this->intervals[i].second < x_min) endIt = 0;
    if (this->intervals[i].second >= x_max) endIt = result.size();
    if ((this->intervals[i].second > x_min) && (this->intervals[i].second < x_max)) {
      endIt = number_of_bins * (this->intervals[i].second - x_min) / (x_max - x_min);
    }

    if (beginIt > endIt) {
      beginIt = endIt;
    }

    if (dbg) {
      std::cerr << "beginIt : " << beginIt << std::endl;
      std::cerr << "endIt : " << endIt << std::endl;
    }

    for (size_t pos = beginIt; pos != endIt; ++pos) {
      result[pos] += ((x_max - x_min) / static_cast<double>(number_of_bins)) *
                     (this->intervals[i].second - this->intervals[i].first);
    }
    if (dbg) {
      std::cerr << "Result at this stage \n";
      for (size_t aa = 0; aa != result.size(); ++aa) {
        std::cerr << result[aa] << " ";
      }
      std::cerr << std::endl;
    }
  }
  return result;
}  // characteristic_function_of_diagram

std::vector<double> Persistence_intervals::cumulative_characteristic_function_of_diagram(double x_min, double x_max,
                                                                                         size_t number_of_bins) const {
  std::vector<double> intsOfBars = this->characteristic_function_of_diagram(x_min, x_max, number_of_bins);
  std::vector<double> result(intsOfBars.size());
  double sum = 0;
  for (size_t i = 0; i != intsOfBars.size(); ++i) {
    sum += intsOfBars[i];
    result[i] = sum;
  }
  return result;
}  // cumulative_characteristic_function_of_diagram

template <typename T>
bool compare_first_element_of_pair(const std::pair<T, bool>& f, const std::pair<T, bool>& s) {
  return (f.first < s.first);
}

std::vector<std::pair<double, size_t> > Persistence_intervals::compute_persistent_betti_numbers() const {
  std::vector<std::pair<double, bool> > places_where_pbs_change(2 * this->intervals.size());

  for (size_t i = 0; i != this->intervals.size(); ++i) {
    places_where_pbs_change[2 * i] = std::make_pair(this->intervals[i].first, true);
    places_where_pbs_change[2 * i + 1] = std::make_pair(this->intervals[i].second, false);
  }

  std::sort(places_where_pbs_change.begin(), places_where_pbs_change.end(), compare_first_element_of_pair<double>);
  size_t pbn = 0;
  std::vector<std::pair<double, size_t> > pbns(places_where_pbs_change.size());
  for (size_t i = 0; i != places_where_pbs_change.size(); ++i) {
    if (places_where_pbs_change[i].second == true) {
      ++pbn;
    } else {
      --pbn;
    }
    pbns[i] = std::make_pair(places_where_pbs_change[i].first, pbn);
  }
  return pbns;
}

inline double compute_euclidean_distance(const std::pair<double, double>& f, const std::pair<double, double>& s) {
  return sqrt((f.first - s.first) * (f.first - s.first) + (f.second - s.second) * (f.second - s.second));
}

std::vector<double> Persistence_intervals::k_n_n(size_t k, size_t where_to_cut) const {
  bool dbg = false;
  if (dbg) {
    std::cerr << "Here are the intervals : \n";
    for (size_t i = 0; i != this->intervals.size(); ++i) {
      std::cerr << "[ " << this->intervals[i].first << " , " << this->intervals[i].second << "] \n";
    }
    getchar();
  }

  std::vector<double> result;
  // compute all to all distance between point in the diagram. Also, consider points in the diagonal with the infinite
  // multiplicity.
  std::vector<std::vector<double> > distances(this->intervals.size());
  for (size_t i = 0; i != this->intervals.size(); ++i) {
    std::vector<double> aa(this->intervals.size());
    std::fill(aa.begin(), aa.end(), 0);
    distances[i] = aa;
  }
  std::vector<double> distances_from_diagonal(this->intervals.size());
  std::fill(distances_from_diagonal.begin(), distances_from_diagonal.end(), 0);

  for (size_t i = 0; i != this->intervals.size(); ++i) {
    std::vector<double> distancesFromI;
    for (size_t j = i + 1; j != this->intervals.size(); ++j) {
      distancesFromI.push_back(compute_euclidean_distance(this->intervals[i], this->intervals[j]));
    }
    // also add a distance from this guy to diagonal:
    double distanceToDiagonal = compute_euclidean_distance(
        this->intervals[i], std::make_pair(0.5 * (this->intervals[i].first + this->intervals[i].second),
                                           0.5 * (this->intervals[i].first + this->intervals[i].second)));
    distances_from_diagonal[i] = distanceToDiagonal;

    if (dbg) {
      std::cerr << "Here are the distances form the point : [" << this->intervals[i].first << " , "
                << this->intervals[i].second << "] in the diagram \n";
      for (size_t aa = 0; aa != distancesFromI.size(); ++aa) {
        std::cerr << "To : " << i + aa << " : " << distancesFromI[aa] << " ";
      }
      std::cerr << std::endl;
      getchar();
    }

    // filling in the distances matrix:
    for (size_t j = i + 1; j != this->intervals.size(); ++j) {
      distances[i][j] = distancesFromI[j - i - 1];
      distances[j][i] = distancesFromI[j - i - 1];
    }
  }
  if (dbg) {
    std::cerr << "Here is the distance matrix : \n";
    for (size_t i = 0; i != distances.size(); ++i) {
      for (size_t j = 0; j != distances.size(); ++j) {
        std::cerr << distances[i][j] << " ";
      }
      std::cerr << std::endl;
    }
    std::cerr << std::endl << std::endl << "And here are the distances to the diagonal : " << std::endl;
    for (size_t i = 0; i != distances_from_diagonal.size(); ++i) {
      std::cerr << distances_from_diagonal[i] << " ";
    }
    std::cerr << std::endl << std::endl;
    getchar();
  }

  for (size_t i = 0; i != this->intervals.size(); ++i) {
    std::vector<double> distancesFromI = distances[i];
    distancesFromI.push_back(distances_from_diagonal[i]);

    // sort it:
    std::sort(distancesFromI.begin(), distancesFromI.end(), std::greater<double>());

    if (k > distancesFromI.size()) {
      if (dbg) {
        std::cerr << "There are not enough neighbors in your set. We set the result to plus infty \n";
      }
      result.push_back(std::numeric_limits<double>::max());
    } else {
      if (distances_from_diagonal[i] > distancesFromI[k]) {
        if (dbg) {
          std::cerr << "The k-th n.n. is on a diagonal. Therefore we set up a distance to diagonal \n";
        }
        result.push_back(distances_from_diagonal[i]);
      } else {
        result.push_back(distancesFromI[k]);
      }
    }
  }
  std::sort(result.begin(), result.end(), std::greater<double>());
  result.resize(std::min(result.size(), where_to_cut));

  return result;
}

double Persistence_intervals::project_to_R(int number_of_function) const {
  double result = 0;

  for (size_t i = 0; i != this->intervals.size(); ++i) {
    result +=
        (this->intervals[i].second - this->intervals[i].first) * (this->intervals[i].second - this->intervals[i].first);
  }

  return result;
}

}  // namespace Persistence_representations
}  // namespace Gudhi

#endif  // PERSISTENCE_INTERVALS_H_
