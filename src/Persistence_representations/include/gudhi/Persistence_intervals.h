/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - 2019/12 Vincent Rouvreau: Fix #118 - Make histogram_of_lengths and cumulative_histogram_of_lengths
 *          return the exact number_of_bins (was failing on x86)
 *      - 2025/06 Hannah Schreiber: Various small bug fixes (missing `inline`s, `DEBUG_TRACES`s etc.)
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PERSISTENCE_INTERVALS_H_
#define PERSISTENCE_INTERVALS_H_

// standard include
#ifdef DEBUG_TRACES
#include <iostream>   // std::clog
#endif
#include <cstddef>    // std::size_t
#include <ostream>    // std::ostream
#include <fstream>    // std::ofstream
#include <sstream>    // std::stringstream
#include <limits>     // std::numeric_limits
#include <algorithm>  // std::sort
#include <cmath>      // std::sqrt
#include <utility>    // std::pair
#include <vector>

// gudhi include
#include <gudhi/read_persistence_from_file.h>
#include <gudhi/Debug_utils.h>

namespace Gudhi {
namespace Persistence_representations {

// TODO: it would have been better to have this file in a subfolder "Persistence_representations"
// to avoid including it with "<gudhi/Persistence_intervals.h>" which makes it sound universal within gudhi
// even though it is only used in this format within this module.
// How critical would it be for retro-compatibility to change that?

/**
 * @class Persistence_intervals Persistence_intervals.h gudhi/Persistence_intervals.h
 * @brief This class implements the following concepts: Vectorized_topological_data, Topological_data_with_distances,
 * Real_valued_topological_data
 *
 * @ingroup Persistence_representations
 **/
class Persistence_intervals
{
 public:
  /**
   * @brief Constructor from a text file.
   *
   * @param filename Each line of the input file is supposed to contain two numbers of a type double (or convertible
   * to double) representing the birth and the death of the persistence interval. If the pairs are not sorted so that
   * birth <= death, then the constructor will sort then that way. It is possible to have a third value at the begining
   * of a line which represent the dimension of the interval.
   * @param dimension If anything other than `std::numeric_limits<unsigned int>::max()`, only the intervals in
   * this given dimension are token into account if the dimension is given in the file lines. Ignored otherwise.
   * Default value: `std::numeric_limits<double>::max()`.
   */
  Persistence_intervals(const char* filename, unsigned int dimension = std::numeric_limits<unsigned int>::max())
      : intervals_(read_persistence_intervals_in_one_dimension_from_file(
            filename,
            dimension == std::numeric_limits<unsigned int>::max() ? -1 : static_cast<int>(dimension))),
        number_of_functions_for_vectorization_(intervals_.size()),
        number_of_functions_for_projections_to_reals_(1)
  {}

  /**
   * @brief Constructor from a vector of pairs.
   *
   * @param intervals Each pair in the vector is assumed to represent a persistence interval \f$ (b, d) \f$,
   * such that \f$ b \leq d \f$.
   */
  Persistence_intervals(const std::vector<std::pair<double, double> >& intervals)
      : intervals_(intervals),
        number_of_functions_for_vectorization_(intervals.size()),
        number_of_functions_for_projections_to_reals_(1)
  {}

  /**
   * This procedure returns x-range of a given persistence diagram.
   **/
  std::pair<double, double> get_x_range() const
  {
    double min = std::numeric_limits<int>::max();
    double max = -std::numeric_limits<int>::max();
    for (std::size_t i = 0; i != this->intervals_.size(); ++i) {
      if (this->intervals_[i].first < min) min = this->intervals_[i].first;
      if (this->intervals_[i].second > max) max = this->intervals_[i].second;
    }
    return std::make_pair(min, max);
  }

  /**
   * This procedure returns y-range of a given persistence diagram.
   **/
  std::pair<double, double> get_y_range() const
  {
    double min = std::numeric_limits<int>::max();
    double max = -std::numeric_limits<int>::max();
    for (std::size_t i = 0; i != this->intervals_.size(); ++i) {
      if (this->intervals_[i].second < min) min = this->intervals_[i].second;
      if (this->intervals_[i].second > max) max = this->intervals_[i].second;
    }
    return std::make_pair(min, max);
  }

  /**
   * Procedure that compute the vector of lengths of the dominant (i.e. the longest) persistence intervals. The list is
   * truncated at the parameter of the call where_to_cut (set by default to 100).
   **/
  std::vector<double> length_of_dominant_intervals(std::size_t where_to_cut = 100) const;

  /**
   * Procedure that compute the vector of the dominant (i.e. the longest) persistence intervals. The parameter of
   * the procedure (set by default to 100) is the number of dominant intervals returned by the procedure.
   **/
  std::vector<std::pair<double, double> > dominant_intervals(std::size_t where_to_cut = 100) const;

  /**
   * Procedure to compute a histogram of interval's length. A histogram is a block plot. The number of blocks is
   * determined by the first parameter of the function (set by default to 10).
   * For the sake of argument let us assume that the length of the longest interval is 1 and the number of bins is
   * 10. In this case the i-th block correspond to a range between i-1/10 and i10.
   * The vale of a block supported at the interval is the number of persistence intervals of a length between x_0
   * and x_1.
   **/
  std::vector<std::size_t> histogram_of_lengths(std::size_t number_of_bins = 10) const;

  /**
   * Based on a histogram of intervals lengths computed by the function histogram_of_lengths H the procedure below
   * computes the cumulative histogram. The i-th position of the resulting histogram is the sum of values of H for
   * the positions from 0 to i.
   **/
  std::vector<std::size_t> cumulative_histogram_of_lengths(std::size_t number_of_bins = 10) const;

  /**
   * In this procedure we assume that each barcode is a characteristic function of a height equal to its length. The
   * persistence diagram is a sum of such a functions. The procedure below construct a function being a sum of the
   * characteristic functions of persistence intervals. The first two parameters are the range in which the function is
   * to be computed and the last parameter is the number of bins in the discretization of the interval [_min,_max].
   **/
  std::vector<double> characteristic_function_of_diagram(double x_min,
                                                         double x_max,
                                                         std::size_t number_of_bins = 10) const;

  /**
   * Cumulative version of the function characteristic_function_of_diagram
   **/
  std::vector<double> cumulative_characteristic_function_of_diagram(double x_min,
                                                                    double x_max,
                                                                    std::size_t number_of_bins = 10) const;

  /**
   * Compute the function of persistence Betti numbers. The returned value is a vector of pair. First element of each
   * pair is a place where persistence Betti numbers change.
   * Second element of each pair is the value of Persistence Betti numbers at that point.
   **/
  std::vector<std::pair<double, std::size_t> > compute_persistent_betti_numbers() const;

  /**
   * This is a non optimal procedure that compute vector of distances from each point of diagram to its k-th nearest
   * neighbor (k is a parameter of the program). The resulting vector is by default truncated to 10 elements
   * (this value can be changed by using the second parameter of the program). The points are returned in order
   * from the ones which are farthest away from their k-th nearest neighbors.
   **/
  std::vector<double> k_n_n(std::size_t k, std::size_t where_to_cut = 10) const;

  /**
   * Operator that send the diagram to a stream.
   **/
  friend std::ostream& operator<<(std::ostream& out, const Persistence_intervals& intervals)
  {
    for (std::size_t i = 0; i != intervals.intervals_.size(); ++i) {
      out << intervals.intervals_[i].first << " " << intervals.intervals_[i].second << std::endl;
    }
    return out;
  }

  /**
   * Generating gnuplot script to plot the interval.
   **/
  void plot(const char* filename,
            double min_x = std::numeric_limits<double>::max(),
            double max_x = std::numeric_limits<double>::max(),
            double min_y = std::numeric_limits<double>::max(),
            double max_y = std::numeric_limits<double>::max()) const
  {
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
    for (std::size_t i = 0; i != this->intervals_.size(); ++i) {
      out << this->intervals_[i].first << " " << this->intervals_[i].second << std::endl;
    }
    out << "EOF" << std::endl;
    out << min_max_values.first - 0.1 * (min_max_values.second - min_max_values.first) << " "
        << min_max_values.first - 0.1 * (min_max_values.second - min_max_values.first) << std::endl;
    out << min_max_values.second + 0.1 * (min_max_values.second - min_max_values.first) << " "
        << min_max_values.second + 0.1 * (min_max_values.second - min_max_values.first) << std::endl;

    out.close();

#ifdef DEBUG_TRACES
    std::clog << "To visualize, install gnuplot and type the command: gnuplot -persist -e \"load \'"
              << gnuplot_script.str().c_str() << "\'\"" << std::endl;
#endif
  }

  /**
   * Return number of points in the diagram.
   **/
  std::size_t size() const { return this->intervals_.size(); }

  /**
   * Return the persistence interval at the given position. Note that intervals are not sorted with respect to their
   * lengths.
   **/
  const std::pair<double, double>& operator[](std::size_t i) const
  {
    if (i >= this->intervals_.size()) throw("Index out of range! Operator [], one_d_gaussians class\n");
    return this->intervals_[i];
  }

  // Implementations of functions for various concepts.

  /**
   * This is a simple function projecting the persistence intervals to a real number. The function we use here is a sum
   * of squared lengths of intervals. It can be naturally interpreted as sum of step function, where the step height
   * it equal to the length of the interval.
   * At the moment this function is not tested, since it is quite likely to be changed in the future. Given this, when
   * using it, keep in mind that it will be most likely changed in the next versions.
   **/
  double project_to_R(int number_of_function) const;

  /**
   * The function gives the number of possible projections to R. This function is required by the
   * Real_valued_topological_data concept.
   **/
  std::size_t number_of_projections_to_R() const { return this->number_of_functions_for_projections_to_reals_; }

  /**
   * Return a family of vectors obtained from the persistence diagram. The i-th vector consist of the length of i
   * dominant persistence intervals.
   **/
  std::vector<double> vectorize(int number_of_function) const
  {
    return this->length_of_dominant_intervals(number_of_function);
  }

  /**
   * This function return the number of functions that allows vectorization of a persistence diagram. It is required
   * in a concept Vectorized_topological_data.
   **/
  std::size_t number_of_vectorize_functions() const { return this->number_of_functions_for_vectorization_; }

  // end of implementation of functions needed for concepts.

  // For visualization use output from vectorize and build histograms.
  const std::vector<std::pair<double, double> >& output_for_visualization() { return this->intervals_; }

 protected:
  std::vector<std::pair<double, double> > intervals_;

 private:
  std::size_t number_of_functions_for_vectorization_;
  std::size_t number_of_functions_for_projections_to_reals_;
};

inline std::vector<double> Persistence_intervals::length_of_dominant_intervals(std::size_t where_to_cut) const
{
  std::vector<double> result(this->intervals_.size());
  for (std::size_t i = 0; i != this->intervals_.size(); ++i) {
    result[i] = this->intervals_[i].second - this->intervals_[i].first;
  }
  std::sort(result.begin(), result.end(), std::greater<double>());

  result.resize(std::min(where_to_cut, result.size()));
  return result;
}

inline std::vector<std::pair<double, double> > Persistence_intervals::dominant_intervals(std::size_t where_to_cut) const
{
  std::vector<std::pair<std::size_t, double> > position_length_vector(this->intervals_.size());
  for (std::size_t i = 0; i != this->intervals_.size(); ++i) {
    position_length_vector[i] = std::make_pair(i, this->intervals_[i].second - this->intervals_[i].first);
  }

  std::sort(position_length_vector.begin(),
            position_length_vector.end(),
            [](const std::pair<std::size_t, double>& first, const std::pair<std::size_t, double>& second) -> bool {
              return first.second > second.second;
            });

  std::vector<std::pair<double, double> > result;
  result.reserve(std::min(where_to_cut, position_length_vector.size()));

  for (std::size_t i = 0; i != std::min(where_to_cut, position_length_vector.size()); ++i) {
    result.push_back(this->intervals_[position_length_vector[i].first]);
#ifdef DEBUG_TRACES
    std::clog << "Position : " << position_length_vector[i].first << " length : " << position_length_vector[i].second
              << std::endl;
#endif
  }

  return result;
}  // dominant_intervals

inline std::vector<std::size_t> Persistence_intervals::histogram_of_lengths(std::size_t number_of_bins) const
{
#ifdef DEBUG_TRACES
  std::clog << "this->intervals.size() : " << this->intervals_.size() << std::endl;
#endif

  // first find the length of the longest interval:
  double lengthOfLongest = 0;
  for (std::size_t i = 0; i != this->intervals_.size(); ++i) {
    if ((this->intervals_[i].second - this->intervals_[i].first) > lengthOfLongest) {
      lengthOfLongest = this->intervals_[i].second - this->intervals_[i].first;
    }
  }

#ifdef DEBUG_TRACES
  std::clog << "lengthOfLongest : " << lengthOfLongest << std::endl;
#endif

  // this is a container we will use to store the resulting histogram
  std::vector<std::size_t> result(number_of_bins + 1, 0);

  // for every persistence interval in our collection.
  for (std::size_t i = 0; i != this->intervals_.size(); ++i) {
    // compute its length relative to the length of the dominant interval:
    double relative_length_of_this_interval =
        (this->intervals_[i].second - this->intervals_[i].first) / lengthOfLongest;

    // given the relative length (between 0 and 1) compute to which bin should it contribute.
    std::size_t position = (std::size_t)(relative_length_of_this_interval * number_of_bins);

    ++result[position];

#ifdef DEBUG_TRACES
    std::clog << "i : " << i << std::endl;
    std::clog << "Interval : [" << this->intervals_[i].first << " , " << this->intervals_[i].second << " ] \n";
    std::clog << "relative_length_of_this_interval : " << relative_length_of_this_interval << std::endl;
    std::clog << "position : " << position << std::endl;
#endif
  }
  // we want number of bins equals to number_of_bins (some unexpected results on x86)
  result[number_of_bins - 1] += result[number_of_bins];
  result.resize(number_of_bins);

#ifdef DEBUG_TRACES
  for (std::size_t i = 0; i != result.size(); ++i) std::clog << result[i] << std::endl;
#endif
  return result;
}

inline std::vector<std::size_t> Persistence_intervals::cumulative_histogram_of_lengths(std::size_t number_of_bins) const
{
  std::vector<std::size_t> histogram = this->histogram_of_lengths(number_of_bins);
  std::vector<std::size_t> result(histogram.size());

  std::size_t sum = 0;
  for (std::size_t i = 0; i != histogram.size(); ++i) {
    sum += histogram[i];
    result[i] = sum;
  }
  return result;
}

inline std::vector<double> Persistence_intervals::characteristic_function_of_diagram(double x_min,
                                                                                     double x_max,
                                                                                     std::size_t number_of_bins) const
{
  std::vector<double> result(number_of_bins);
  std::fill(result.begin(), result.end(), 0);

  for (std::size_t i = 0; i != this->intervals_.size(); ++i) {
#ifdef DEBUG_TRACES
    std::clog << "Interval : " << this->intervals_[i].first << " , " << this->intervals_[i].second << std::endl;
#endif

    std::size_t beginIt = 0;
    if (this->intervals_[i].first < x_min) beginIt = 0;
    if (this->intervals_[i].first >= x_max) beginIt = result.size();
    if ((this->intervals_[i].first > x_min) && (this->intervals_[i].first < x_max)) {
      beginIt = number_of_bins * (this->intervals_[i].first - x_min) / (x_max - x_min);
    }

    std::size_t endIt = 0;
    if (this->intervals_[i].second < x_min) endIt = 0;
    if (this->intervals_[i].second >= x_max) endIt = result.size();
    if ((this->intervals_[i].second > x_min) && (this->intervals_[i].second < x_max)) {
      endIt = number_of_bins * (this->intervals_[i].second - x_min) / (x_max - x_min);
    }

    if (beginIt > endIt) {
      beginIt = endIt;
    }

#ifdef DEBUG_TRACES
    std::clog << "beginIt : " << beginIt << std::endl;
    std::clog << "endIt : " << endIt << std::endl;
#endif

    for (std::size_t pos = beginIt; pos != endIt; ++pos) {
      result[pos] += ((x_max - x_min) / static_cast<double>(number_of_bins)) *
                     (this->intervals_[i].second - this->intervals_[i].first);
    }
#ifdef DEBUG_TRACES
    std::clog << "Result at this stage \n";
    for (std::size_t aa = 0; aa != result.size(); ++aa) {
      std::clog << result[aa] << " ";
    }
    std::clog << std::endl;
#endif
  }
  return result;
}  // characteristic_function_of_diagram

inline std::vector<double> Persistence_intervals::cumulative_characteristic_function_of_diagram(
    double x_min,
    double x_max,
    std::size_t number_of_bins) const
{
  std::vector<double> intsOfBars = this->characteristic_function_of_diagram(x_min, x_max, number_of_bins);
  std::vector<double> result(intsOfBars.size());
  double sum = 0;
  for (std::size_t i = 0; i != intsOfBars.size(); ++i) {
    sum += intsOfBars[i];
    result[i] = sum;
  }
  return result;
}

inline std::vector<std::pair<double, std::size_t> > Persistence_intervals::compute_persistent_betti_numbers() const
{
  std::vector<std::pair<double, bool> > places_where_pbs_change(2 * this->intervals_.size());

  for (std::size_t i = 0; i != this->intervals_.size(); ++i) {
    places_where_pbs_change[2 * i] = std::make_pair(this->intervals_[i].first, true);
    places_where_pbs_change[2 * i + 1] = std::make_pair(this->intervals_[i].second, false);
  }

  std::sort(
      places_where_pbs_change.begin(),
      places_where_pbs_change.end(),
      [](const std::pair<double, bool>& f, const std::pair<double, bool>& s) -> bool { return f.first < s.first; });

  std::size_t pbn = 0;
  std::vector<std::pair<double, std::size_t> > pbns(places_where_pbs_change.size());
  for (std::size_t i = 0; i != places_where_pbs_change.size(); ++i) {
    if (places_where_pbs_change[i].second == true) {
      ++pbn;
    } else {
      --pbn;
    }
    pbns[i] = std::make_pair(places_where_pbs_change[i].first, pbn);
  }
  return pbns;
}

inline std::vector<double> Persistence_intervals::k_n_n(std::size_t k, std::size_t where_to_cut) const
{
#ifdef DEBUG_TRACES
  std::clog << "Here are the intervals : \n";
  for (std::size_t i = 0; i != this->intervals_.size(); ++i) {
    std::clog << "[ " << this->intervals_[i].first << " , " << this->intervals_[i].second << "] \n";
  }
#endif

  auto compute_euclidean_distance = [](const std::pair<double, double>& f,
                                       const std::pair<double, double>& s) -> double {
    return std::sqrt((f.first - s.first) * (f.first - s.first) + (f.second - s.second) * (f.second - s.second));
  };

  std::vector<double> result;
  // compute all to all distance between point in the diagram. Also, consider points in the diagonal with the infinite
  // multiplicity.
  std::vector<std::vector<double> > distances(this->intervals_.size());
  for (std::size_t i = 0; i != this->intervals_.size(); ++i) {
    std::vector<double> aa(this->intervals_.size());
    std::fill(aa.begin(), aa.end(), 0);
    distances[i] = aa;
  }
  std::vector<double> distances_from_diagonal(this->intervals_.size());
  std::fill(distances_from_diagonal.begin(), distances_from_diagonal.end(), 0);

  for (std::size_t i = 0; i != this->intervals_.size(); ++i) {
    std::vector<double> distancesFromI;
    for (std::size_t j = i + 1; j != this->intervals_.size(); ++j) {
      distancesFromI.push_back(compute_euclidean_distance(this->intervals_[i], this->intervals_[j]));
    }
    // also add a distance from this guy to diagonal:
    double distanceToDiagonal =
        compute_euclidean_distance(this->intervals_[i],
                                   std::make_pair(0.5 * (this->intervals_[i].first + this->intervals_[i].second),
                                                  0.5 * (this->intervals_[i].first + this->intervals_[i].second)));
    distances_from_diagonal[i] = distanceToDiagonal;

#ifdef DEBUG_TRACES
    std::clog << "Here are the distances form the point : [" << this->intervals_[i].first << " , "
              << this->intervals_[i].second << "] in the diagram \n";
    for (std::size_t aa = 0; aa != distancesFromI.size(); ++aa) {
      std::clog << "To : " << i + aa << " : " << distancesFromI[aa] << " ";
    }
    std::clog << std::endl;
#endif

    // filling in the distances matrix:
    for (std::size_t j = i + 1; j != this->intervals_.size(); ++j) {
      distances[i][j] = distancesFromI[j - i - 1];
      distances[j][i] = distancesFromI[j - i - 1];
    }
  }

#ifdef DEBUG_TRACES
  std::clog << "Here is the distance matrix : \n";
  for (std::size_t i = 0; i != distances.size(); ++i) {
    for (std::size_t j = 0; j != distances.size(); ++j) {
      std::clog << distances[i][j] << " ";
    }
    std::clog << std::endl;
  }
  std::clog << std::endl << std::endl << "And here are the distances to the diagonal : " << std::endl;
  for (std::size_t i = 0; i != distances_from_diagonal.size(); ++i) {
    std::clog << distances_from_diagonal[i] << " ";
  }
  std::clog << std::endl << std::endl;
#endif

  for (std::size_t i = 0; i != this->intervals_.size(); ++i) {
    std::vector<double> distancesFromI = distances[i];
    distancesFromI.push_back(distances_from_diagonal[i]);

    // sort it:
    std::sort(distancesFromI.begin(), distancesFromI.end(), std::greater<double>());

    if (k > distancesFromI.size()) {
#ifdef DEBUG_TRACES
      std::clog << "There are not enough neighbors in your set. We set the result to plus infty \n";
#endif
      result.push_back(std::numeric_limits<double>::max());
    } else {
      if (distances_from_diagonal[i] > distancesFromI[k]) {
#ifdef DEBUG_TRACES
        std::clog << "The k-th n.n. is on a diagonal. Therefore we set up a distance to diagonal \n";
#endif
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

inline double Persistence_intervals::project_to_R(int number_of_function) const
{
  double result = 0;

  for (std::size_t i = 0; i != this->intervals_.size(); ++i) {
    auto diff = this->intervals_[i].second - this->intervals_[i].first;
    result += diff * diff;
  }

  return result;
}

}  // namespace Persistence_representations
}  // namespace Gudhi

#endif  // PERSISTENCE_INTERVALS_H_
