/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PERSISTENCE_LANDSCAPE_H_
#define PERSISTENCE_LANDSCAPE_H_

// gudhi include
#include <gudhi/read_persistence_from_file.h>
#include <gudhi/common_persistence_representations.h>

// standard include
#include <cmath>
#include <iostream>
#include <vector>
#include <limits>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include <utility>
#include <functional>

namespace Gudhi {
namespace Persistence_representations {

// pre declaration
class Persistence_landscape;
template <typename operation>
Persistence_landscape operation_on_pair_of_landscapes(const Persistence_landscape& land1,
                                                      const Persistence_landscape& land2);

/**
 * \class Persistence_landscape Persistence_landscape.h gudhi/Persistence_landscape.h
 * \brief A class implementing persistence landscapes data structures.
 *
 * \ingroup Persistence_representations
 *
 * \details
 * For theoretical description, please consult <i>Statistical topological data analysis using persistence
 * landscapes</i>\cite bubenik_landscapes_2015 , and for details of algorithms,
 * <i>A persistence landscapes toolbox for topological statistics</i>\cite bubenik_dlotko_landscapes_2016.
 *
 * Persistence landscapes allow vectorization, computations of distances, computations of projections to Real,
 * computations of averages and scalar products. Therefore they implement suitable interfaces.
 * It implements the following concepts: Vectorized_topological_data, Topological_data_with_distances,
 * Real_valued_topological_data, Topological_data_with_averages, Topological_data_with_scalar_product
 *
 * Note that at the moment, due to rounding errors during the construction of persistence landscapes, elements which
 * are different by 0.000005 are considered the same. If the scale in your persistence diagrams is comparable to this
 * value, please rescale them before use this code.
 *
**/
class Persistence_landscape {
 public:
  /**
   * Default constructor.
  **/
  Persistence_landscape() { this->set_up_numbers_of_functions_for_vectorization_and_projections_to_reals(); }

  /**
  * Constructor that takes as an input a vector of birth-death pairs.
  **/
  Persistence_landscape(const std::vector<std::pair<double, double> >& p,
                        size_t number_of_levels = std::numeric_limits<size_t>::max());

  /**
       * Constructor that reads persistence intervals from file and creates persistence landscape. The format of the
    *input file is the following: in each line we put birth-death pair. Last line is assumed
       * to be empty. Even if the points within a line are not ordered, they will be ordered while the input is read.
      **/
  Persistence_landscape(const char* filename, size_t dimension = std::numeric_limits<unsigned>::max(),
                        size_t number_of_levels = std::numeric_limits<size_t>::max());

  /**
   * This procedure loads a landscape from file. It erase all the data that was previously stored in this landscape.
  **/
  void load_landscape_from_file(const char* filename);

  /**
   * The procedure stores a landscape to a file. The file can be later used by a procedure load_landscape_from_file.
  **/
  void print_to_file(const char* filename) const;

  /**
   * This function compute integral of the landscape (defined formally as sum of integrals on R of all landscape
   *functions)
  **/
  double compute_integral_of_landscape() const;

  /**
       * This function compute integral of the 'level'-level of a landscape.
      **/
  double compute_integral_of_a_level_of_a_landscape(size_t level) const;

  /**
       * This function compute integral of the landscape p-th power of a landscape (defined formally as sum of integrals
    *on R of p-th powers of all landscape functions)
      **/
  double compute_integral_of_landscape(double p) const;  // this function compute integral of p-th power of landscape.

  /**
   * A function that computes the value of a landscape at a given point. The parameters of the function are: unsigned
   *level and double x.
   * The procedure will compute the value of the level-landscape at the point x.
  **/
  double compute_value_at_a_given_point(unsigned level, double x) const;

  /**
   * Writing landscape into a stream. A i-th level landscape starts with a string "lambda_i". Then the discontinuity
   *points of the landscapes follows.
   * Shall those points be joined with lines, we will obtain the i-th landscape function.
  **/
  friend std::ostream& operator<<(std::ostream& out, Persistence_landscape& land);

  template <typename operation>
  friend Persistence_landscape operation_on_pair_of_landscapes(const Persistence_landscape& land1,
                                                               const Persistence_landscape& land2);

  /**
   *\private A function that compute sum of two landscapes.
  **/
  friend Persistence_landscape add_two_landscapes(const Persistence_landscape& land1,
                                                  const Persistence_landscape& land2) {
    return operation_on_pair_of_landscapes<std::plus<double> >(land1, land2);
  }

  /**
       *\private A function that compute difference of two landscapes.
      **/
  friend Persistence_landscape subtract_two_landscapes(const Persistence_landscape& land1,
                                                       const Persistence_landscape& land2) {
    return operation_on_pair_of_landscapes<std::minus<double> >(land1, land2);
  }

  /**
   * An operator +, that compute sum of two landscapes.
  **/
  friend Persistence_landscape operator+(const Persistence_landscape& first, const Persistence_landscape& second) {
    return add_two_landscapes(first, second);
  }

  /**
   * An operator -, that compute difference of two landscapes.
  **/
  friend Persistence_landscape operator-(const Persistence_landscape& first, const Persistence_landscape& second) {
    return subtract_two_landscapes(first, second);
  }

  /**
   * An operator * that allows multiplication of a landscape by a real number.
  **/
  friend Persistence_landscape operator*(const Persistence_landscape& first, double con) {
    return first.multiply_lanscape_by_real_number_not_overwrite(con);
  }

  /**
   * An operator * that allows multiplication of a landscape by a real number (order of parameters swapped).
  **/
  friend Persistence_landscape operator*(double con, const Persistence_landscape& first) {
    return first.multiply_lanscape_by_real_number_not_overwrite(con);
  }

  /**
   * Operator +=. The second parameter is persistence landscape.
  **/
  Persistence_landscape operator+=(const Persistence_landscape& rhs) {
    *this = *this + rhs;
    return *this;
  }

  /**
   * Operator -=. The second parameter is a persistence landscape.
  **/
  Persistence_landscape operator-=(const Persistence_landscape& rhs) {
    *this = *this - rhs;
    return *this;
  }

  /**
   * Operator *=. The second parameter is a real number by which the y values of all landscape functions are multiplied.
   *The x-values remain unchanged.
  **/
  Persistence_landscape operator*=(double x) {
    *this = *this * x;
    return *this;
  }

  /**
   * Operator /=. The second parameter is a real number.
  **/
  Persistence_landscape operator/=(double x) {
    if (x == 0) throw("In operator /=, division by 0. Program terminated.");
    *this = *this * (1 / x);
    return *this;
  }

  /**
   * An operator to compare two persistence landscapes.
  **/
  bool operator==(const Persistence_landscape& rhs) const;

  /**
       * An operator to compare two persistence landscapes.
      **/
  bool operator!=(const Persistence_landscape& rhs) const { return !((*this) == rhs); }

  /**
   * Computations of maximum (y) value of landscape.
  **/
  double compute_maximum() const {
    double maxValue = 0;
    if (this->land.size()) {
      maxValue = -std::numeric_limits<int>::max();
      for (size_t i = 0; i != this->land[0].size(); ++i) {
        if (this->land[0][i].second > maxValue) maxValue = this->land[0][i].second;
      }
    }
    return maxValue;
  }

  /**
       *\private Computations of minimum (y) value of landscape.
      **/
  double compute_minimum() const {
    double minValue = 0;
    if (this->land.size()) {
      minValue = std::numeric_limits<int>::max();
      for (size_t i = 0; i != this->land[0].size(); ++i) {
        if (this->land[0][i].second < minValue) minValue = this->land[0][i].second;
      }
    }
    return minValue;
  }

  /**
   *\private Computations of a \f$L^i\f$ norm of landscape, where i is the input parameter.
  **/
  double compute_norm_of_landscape(double i) {
    Persistence_landscape l;
    if (i < std::numeric_limits<double>::max()) {
      return compute_distance_of_landscapes(*this, l, i);
    } else {
      return compute_max_norm_distance_of_landscapes(*this, l);
    }
  }

  /**
   * An operator to compute the value of a landscape in the level 'level' at the argument 'x'.
  **/
  double operator()(unsigned level, double x) const { return this->compute_value_at_a_given_point(level, x); }

  /**
   *\private Computations of \f$L^{\infty}\f$ distance between two landscapes.
  **/
  friend double compute_max_norm_distance_of_landscapes(const Persistence_landscape& first,
                                                        const Persistence_landscape& second);

  /**
   *\private Computations of \f$L^{p}\f$ distance between two landscapes. p is the parameter of the procedure.
  **/
  friend double compute_distance_of_landscapes(const Persistence_landscape& first, const Persistence_landscape& second,
                                               double p);

  /**
   * Function to compute absolute value of a PL function. The representation of persistence landscapes allow to store
   *general PL-function. When computing distance between two landscapes, we compute difference between
   * them. In this case, a general PL-function with negative value can appear as a result. Then in order to compute
   *distance, we need to take its absolute value. This is the purpose of this procedure.
  **/
  Persistence_landscape abs();

  Persistence_landscape* new_abs();

  /**
   * Computes the number of landscape functions.
  **/
  size_t size() const { return this->land.size(); }

  /**
   *  Compute maximal value of lambda-level landscape.
  **/
  double find_max(unsigned lambda) const;

  /**
   *\private Function to compute inner (scalar) product of two landscapes.
  **/
  friend double compute_inner_product(const Persistence_landscape& l1, const Persistence_landscape& l2);

  // Implementations of functions for various concepts.

  /**
   * The number of projections to R is defined to the number of nonzero landscape functions. I-th projection is an
   *integral of i-th landscape function over whole R.
   * This function is required by the Real_valued_topological_data concept.
   * At the moment this function is not tested, since it is quite likely to be changed in the future. Given this, when
   *using it, keep in mind that it
   * will be most likely changed in the next versions.
  **/
  double project_to_R(int number_of_function) const {
    return this->compute_integral_of_a_level_of_a_landscape((size_t)number_of_function);
  }

  /**
   * The function gives the number of possible projections to R. This function is required by the
   *Real_valued_topological_data concept.
  **/
  size_t number_of_projections_to_R() const { return this->number_of_functions_for_projections_to_reals; }

  /**
   * This function produce a vector of doubles based on a landscape. It is required in a concept
   * Vectorized_topological_data
  */
  std::vector<double> vectorize(int number_of_function) const {
    // TODO(PD) think of something smarter over here
    std::vector<double> v;
    if ((size_t)number_of_function > this->land.size()) {
      return v;
    }
    v.reserve(this->land[number_of_function].size());
    for (size_t i = 0; i != this->land[number_of_function].size(); ++i) {
      v.push_back(this->land[number_of_function][i].second);
    }
    return v;
  }
  /**
   * This function return the number of functions that allows vectorization of persistence landscape. It is required in
   *a concept Vectorized_topological_data.
   **/
  size_t number_of_vectorize_functions() const { return this->number_of_functions_for_vectorization; }

  /**
   * A function to compute averaged persistence landscape, based on vector of persistence landscapes.
   * This function is required by Topological_data_with_averages concept.
  **/
  void compute_average(const std::vector<Persistence_landscape*>& to_average) {
    bool dbg = false;

    if (dbg) {
      std::clog << "to_average.size() : " << to_average.size() << std::endl;
    }

    std::vector<Persistence_landscape*> nextLevelMerge(to_average.size());
    for (size_t i = 0; i != to_average.size(); ++i) {
      nextLevelMerge[i] = to_average[i];
    }
    bool is_this_first_level = true;  // in the loop, we will create dynamically a number of intermediate complexes. We
                                      // have to clean that up, but we cannot erase the initial landscapes we have
    // to average. In this case, we simply check if the nextLevelMerge are the input landscapes or the ones created in
    // that loop by using this extra variable.

    while (nextLevelMerge.size() != 1) {
      if (dbg) {
        std::clog << "nextLevelMerge.size() : " << nextLevelMerge.size() << std::endl;
      }
      std::vector<Persistence_landscape*> nextNextLevelMerge;
      nextNextLevelMerge.reserve(to_average.size());
      for (size_t i = 0; i < nextLevelMerge.size(); i = i + 2) {
        if (dbg) {
          std::clog << "i : " << i << std::endl;
        }
        Persistence_landscape* l = new Persistence_landscape;
        if (i + 1 != nextLevelMerge.size()) {
          (*l) = (*nextLevelMerge[i]) + (*nextLevelMerge[i + 1]);
        } else {
          (*l) = *nextLevelMerge[i];
        }
        nextNextLevelMerge.push_back(l);
      }
      if (dbg) {
        std::clog << "After this iteration \n";
        getchar();
      }

      if (!is_this_first_level) {
        // deallocate the memory if the vector nextLevelMerge do not consist of the initial landscapes
        for (size_t i = 0; i != nextLevelMerge.size(); ++i) {
          delete nextLevelMerge[i];
        }
      }
      is_this_first_level = false;
      nextLevelMerge.swap(nextNextLevelMerge);
    }
    (*this) = (*nextLevelMerge[0]);
    (*this) *= 1 / static_cast<double>(to_average.size());
  }

  /**
  * A function to compute distance between persistence landscape.
  * The parameter of this function is a Persistence_landscape.
  * This function is required in Topological_data_with_distances concept.
  * For max norm distance, set power to std::numeric_limits<double>::max()
  **/
  double distance(const Persistence_landscape& second, double power = 1) const {
    if (power < std::numeric_limits<double>::max()) {
      return compute_distance_of_landscapes(*this, second, power);
    } else {
      return compute_max_norm_distance_of_landscapes(*this, second);
    }
  }

  /**
  * A function to compute scalar product of persistence landscapes.
  * The parameter of this function is a Persistence_landscape.
  * This function is required in Topological_data_with_scalar_product concept.
  **/
  double compute_scalar_product(const Persistence_landscape& second) const {
    return compute_inner_product((*this), second);
  }
  // end of implementation of functions needed for concepts.

  /**
   * This procedure returns y-range of a given level persistence landscape. If a default value is used, the y-range
   * of 0th level landscape is given (and this range contains the ranges of all other landscapes).
  **/
  std::pair<double, double> get_y_range(size_t level = 0) const {
    std::pair<double, double> result;
    if (level < this->land.size()) {
      double maxx = this->compute_maximum();
      double minn = this->compute_minimum();
      result = std::make_pair(minn, maxx);
    } else {
      result = std::make_pair(0, 0);
    }
    return result;
  }

  // a function used to create a gnuplot script for visualization of landscapes
  void plot(const char* filename, double xRangeBegin = std::numeric_limits<double>::max(),
            double xRangeEnd = std::numeric_limits<double>::max(),
            double yRangeBegin = std::numeric_limits<double>::max(),
            double yRangeEnd = std::numeric_limits<double>::max(), int from = std::numeric_limits<int>::max(),
            int to = std::numeric_limits<int>::max());

 protected:
  std::vector<std::vector<std::pair<double, double> > > land;
  size_t number_of_functions_for_vectorization;
  size_t number_of_functions_for_projections_to_reals;

  void construct_persistence_landscape_from_barcode(const std::vector<std::pair<double, double> >& p,
                                                    size_t number_of_levels = std::numeric_limits<size_t>::max());
  Persistence_landscape multiply_lanscape_by_real_number_not_overwrite(double x) const;
  void multiply_lanscape_by_real_number_overwrite(double x);
  friend double compute_maximal_distance_non_symmetric(const Persistence_landscape& pl1,
                                                       const Persistence_landscape& pl2);

  void set_up_numbers_of_functions_for_vectorization_and_projections_to_reals() {
    // warning, this function can be only called after filling in the intervals vector.
    this->number_of_functions_for_vectorization = this->land.size();
    this->number_of_functions_for_projections_to_reals = this->land.size();
  }
};

Persistence_landscape::Persistence_landscape(const char* filename, size_t dimension, size_t number_of_levels) {
  std::vector<std::pair<double, double> > barcode;
  if (dimension < std::numeric_limits<double>::max()) {
    barcode = read_persistence_intervals_in_one_dimension_from_file(filename, dimension);
  } else {
    barcode = read_persistence_intervals_in_one_dimension_from_file(filename);
  }
  this->construct_persistence_landscape_from_barcode(barcode, number_of_levels);
  this->set_up_numbers_of_functions_for_vectorization_and_projections_to_reals();
}

bool operatorEqualDbg = false;
bool Persistence_landscape::operator==(const Persistence_landscape& rhs) const {
  if (this->land.size() != rhs.land.size()) {
    if (operatorEqualDbg) std::clog << "1\n";
    return false;
  }
  for (size_t level = 0; level != this->land.size(); ++level) {
    if (this->land[level].size() != rhs.land[level].size()) {
      if (operatorEqualDbg) std::clog << "this->land[level].size() : " << this->land[level].size() << "\n";
      if (operatorEqualDbg) std::clog << "rhs.land[level].size() : " << rhs.land[level].size() << "\n";
      if (operatorEqualDbg) std::clog << "2\n";
      return false;
    }
    for (size_t i = 0; i != this->land[level].size(); ++i) {
      if (!(almost_equal(this->land[level][i].first, rhs.land[level][i].first) &&
            almost_equal(this->land[level][i].second, rhs.land[level][i].second))) {
        if (operatorEqualDbg)
          std::clog << "this->land[level][i] : " << this->land[level][i].first << " " << this->land[level][i].second
                    << "\n";
        if (operatorEqualDbg)
          std::clog << "rhs.land[level][i] : " << rhs.land[level][i].first << " " << rhs.land[level][i].second << "\n";
        if (operatorEqualDbg) std::clog << "3\n";
        return false;
      }
    }
  }
  return true;
}

Persistence_landscape::Persistence_landscape(const std::vector<std::pair<double, double> >& p,
                                             size_t number_of_levels) {
  this->construct_persistence_landscape_from_barcode(p, number_of_levels);
  this->set_up_numbers_of_functions_for_vectorization_and_projections_to_reals();
}

void Persistence_landscape::construct_persistence_landscape_from_barcode(
    const std::vector<std::pair<double, double> >& p, size_t number_of_levels) {
  bool dbg = false;
  if (dbg) {
    std::clog << "Persistence_landscape::Persistence_landscape( const std::vector< std::pair< double , double > >& p )"
              << std::endl;
  }

  // this is a general algorithm to construct persistence landscapes.
  std::vector<std::pair<double, double> > bars;
  bars.insert(bars.begin(), p.begin(), p.end());
  std::sort(bars.begin(), bars.end(), compare_points_sorting);

  if (dbg) {
    std::clog << "Bars : \n";
    for (size_t i = 0; i != bars.size(); ++i) {
      std::clog << bars[i].first << " " << bars[i].second << "\n";
    }
    getchar();
  }

  std::vector<std::pair<double, double> > characteristicPoints(p.size());
  for (size_t i = 0; i != bars.size(); ++i) {
    characteristicPoints[i] =
        std::make_pair((bars[i].first + bars[i].second) / 2.0, (bars[i].second - bars[i].first) / 2.0);
  }
  std::vector<std::vector<std::pair<double, double> > > Persistence_landscape;
  size_t number_of_levels_in_the_landscape = 0;
  while (!characteristicPoints.empty()) {
    if (dbg) {
      for (size_t i = 0; i != characteristicPoints.size(); ++i) {
        std::clog << "(" << characteristicPoints[i].first << " " << characteristicPoints[i].second << ")\n";
      }
      std::cin.ignore();
    }

    std::vector<std::pair<double, double> > lambda_n;
    lambda_n.push_back(std::make_pair(-std::numeric_limits<int>::max(), 0));
    lambda_n.push_back(std::make_pair(minus_length(characteristicPoints[0]), 0));
    lambda_n.push_back(characteristicPoints[0]);

    if (dbg) {
      std::clog << "1 Adding to lambda_n : (" << -std::numeric_limits<int>::max() << " " << 0 << ") , ("
                << minus_length(characteristicPoints[0]) << " " << 0 << ") , (" << characteristicPoints[0].first << " "
                << characteristicPoints[0].second << ") \n";
    }

    size_t i = 1;
    std::vector<std::pair<double, double> > newCharacteristicPoints;
    while (i < characteristicPoints.size()) {
      size_t p = 1;
      if ((minus_length(characteristicPoints[i]) >= minus_length(lambda_n[lambda_n.size() - 1])) &&
          (birth_plus_deaths(characteristicPoints[i]) > birth_plus_deaths(lambda_n[lambda_n.size() - 1]))) {
        if (minus_length(characteristicPoints[i]) < birth_plus_deaths(lambda_n[lambda_n.size() - 1])) {
          std::pair<double, double> point = std::make_pair(
              (minus_length(characteristicPoints[i]) + birth_plus_deaths(lambda_n[lambda_n.size() - 1])) / 2,
              (birth_plus_deaths(lambda_n[lambda_n.size() - 1]) - minus_length(characteristicPoints[i])) / 2);
          lambda_n.push_back(point);
          if (dbg) {
            std::clog << "2 Adding to lambda_n : (" << point.first << " " << point.second << ")\n";
          }

          if (dbg) {
            std::clog << "characteristicPoints[i+p] : " << characteristicPoints[i + p].first << " "
                      << characteristicPoints[i + p].second << "\n";
            std::clog << "point : " << point.first << " " << point.second << "\n";
            getchar();
          }

          while ((i + p < characteristicPoints.size()) &&
                 (almost_equal(minus_length(point), minus_length(characteristicPoints[i + p]))) &&
                 (birth_plus_deaths(point) <= birth_plus_deaths(characteristicPoints[i + p]))) {
            newCharacteristicPoints.push_back(characteristicPoints[i + p]);
            if (dbg) {
              std::clog << "3.5 Adding to newCharacteristicPoints : (" << characteristicPoints[i + p].first << " "
                        << characteristicPoints[i + p].second << ")\n";
              getchar();
            }
            ++p;
          }

          newCharacteristicPoints.push_back(point);
          if (dbg) {
            std::clog << "4 Adding to newCharacteristicPoints : (" << point.first << " " << point.second << ")\n";
          }

          while ((i + p < characteristicPoints.size()) &&
                 (minus_length(point) <= minus_length(characteristicPoints[i + p])) &&
                 (birth_plus_deaths(point) >= birth_plus_deaths(characteristicPoints[i + p]))) {
            newCharacteristicPoints.push_back(characteristicPoints[i + p]);
            if (dbg) {
              std::clog << "characteristicPoints[i+p] : " << characteristicPoints[i + p].first << " "
                        << characteristicPoints[i + p].second << "\n";
              std::clog << "point : " << point.first << " " << point.second << "\n";
              std::clog << "characteristicPoints[i+p] birth and death : " << minus_length(characteristicPoints[i + p])
                        << " , " << birth_plus_deaths(characteristicPoints[i + p]) << "\n";
              std::clog << "point birth and death : " << minus_length(point) << " , " << birth_plus_deaths(point)
                        << "\n";

              std::clog << "3 Adding to newCharacteristicPoints : (" << characteristicPoints[i + p].first << " "
                        << characteristicPoints[i + p].second << ")\n";
              getchar();
            }
            ++p;
          }

        } else {
          lambda_n.push_back(std::make_pair(birth_plus_deaths(lambda_n[lambda_n.size() - 1]), 0));
          lambda_n.push_back(std::make_pair(minus_length(characteristicPoints[i]), 0));
          if (dbg) {
            std::clog << "5 Adding to lambda_n : (" << birth_plus_deaths(lambda_n[lambda_n.size() - 1]) << " " << 0
                      << ")\n";
            std::clog << "5 Adding to lambda_n : (" << minus_length(characteristicPoints[i]) << " " << 0 << ")\n";
          }
        }
        lambda_n.push_back(characteristicPoints[i]);
        if (dbg) {
          std::clog << "6 Adding to lambda_n : (" << characteristicPoints[i].first << " "
                    << characteristicPoints[i].second << ")\n";
        }
      } else {
        newCharacteristicPoints.push_back(characteristicPoints[i]);
        if (dbg) {
          std::clog << "7 Adding to newCharacteristicPoints : (" << characteristicPoints[i].first << " "
                    << characteristicPoints[i].second << ")\n";
        }
      }
      i = i + p;
    }
    lambda_n.push_back(std::make_pair(birth_plus_deaths(lambda_n[lambda_n.size() - 1]), 0));
    lambda_n.push_back(std::make_pair(std::numeric_limits<int>::max(), 0));

    characteristicPoints = newCharacteristicPoints;

    lambda_n.erase(std::unique(lambda_n.begin(), lambda_n.end()), lambda_n.end());
    this->land.push_back(lambda_n);

    ++number_of_levels_in_the_landscape;
    if (number_of_levels == number_of_levels_in_the_landscape) {
      break;
    }
  }
}

// this function find maximum of lambda_n
double Persistence_landscape::find_max(unsigned lambda) const {
  if (this->land.size() < lambda) return 0;
  double maximum = -std::numeric_limits<int>::max();
  for (size_t i = 0; i != this->land[lambda].size(); ++i) {
    if (this->land[lambda][i].second > maximum) maximum = this->land[lambda][i].second;
  }
  return maximum;
}

double Persistence_landscape::compute_integral_of_landscape() const {
  double result = 0;
  for (size_t i = 0; i != this->land.size(); ++i) {
    for (size_t nr = 2; nr != this->land[i].size() - 1; ++nr) {
      // it suffices to compute every planar integral and then sum them up for each lambda_n
      result += 0.5 * (this->land[i][nr].first - this->land[i][nr - 1].first) *
                (this->land[i][nr].second + this->land[i][nr - 1].second);
    }
  }
  return result;
}

double Persistence_landscape::compute_integral_of_a_level_of_a_landscape(size_t level) const {
  double result = 0;
  if (level >= this->land.size()) {
    // this landscape function is constantly equal 0, so is the integral.
    return result;
  }
  // also negative landscapes are assumed to be zero.
  if (level < 0) return 0;

  for (size_t nr = 2; nr != this->land[level].size() - 1; ++nr) {
    // it suffices to compute every planar integral and then sum them up for each lambda_n
    result += 0.5 * (this->land[level][nr].first - this->land[level][nr - 1].first) *
              (this->land[level][nr].second + this->land[level][nr - 1].second);
  }

  return result;
}

double Persistence_landscape::compute_integral_of_landscape(double p) const {
  bool dbg = false;
  double result = 0;
  for (size_t i = 0; i != this->land.size(); ++i) {
    for (size_t nr = 2; nr != this->land[i].size() - 1; ++nr) {
      if (dbg) std::clog << "nr : " << nr << "\n";
      // In this interval, the landscape has a form f(x) = ax+b. We want to compute integral of (ax+b)^p = 1/a *
      // (ax+b)^{p+1}/(p+1)
      std::pair<double, double> coef = compute_parameters_of_a_line(this->land[i][nr], this->land[i][nr - 1]);
      double a = coef.first;
      double b = coef.second;

      if (dbg)
        std::clog << "(" << this->land[i][nr].first << "," << this->land[i][nr].second << ") , "
                  << this->land[i][nr - 1].first << "," << this->land[i][nr].second << ")" << std::endl;
      if (this->land[i][nr].first == this->land[i][nr - 1].first) continue;
      if (a != 0) {
        result += 1 / (a * (p + 1)) *
                  (pow((a * this->land[i][nr].first + b), p + 1) - pow((a * this->land[i][nr - 1].first + b), p + 1));
      } else {
        result += (this->land[i][nr].first - this->land[i][nr - 1].first) * (pow(this->land[i][nr].second, p));
      }
      if (dbg) {
        std::clog << "a : " << a << " , b : " << b << std::endl;
        std::clog << "result : " << result << std::endl;
      }
    }
  }
  return result;
}

// this is O(log(n)) algorithm, where n is number of points in this->land.
double Persistence_landscape::compute_value_at_a_given_point(unsigned level, double x) const {
  bool compute_value_at_a_given_pointDbg = false;
  // in such a case lambda_level = 0.
  if (level >= this->land.size()) return 0;

  // we know that the points in this->land[level] are ordered according to x coordinate. Therefore, we can find the
  // point by using bisection:
  unsigned coordBegin = 1;
  unsigned coordEnd = this->land[level].size() - 2;

  if (compute_value_at_a_given_pointDbg) {
    std::clog << "Here \n";
    std::clog << "x : " << x << "\n";
    std::clog << "this->land[level][coordBegin].first : " << this->land[level][coordBegin].first << "\n";
    std::clog << "this->land[level][coordEnd].first : " << this->land[level][coordEnd].first << "\n";
  }

  // in this case x is outside the support of the landscape, therefore the value of the landscape is 0.
  if (x <= this->land[level][coordBegin].first) return 0;
  if (x >= this->land[level][coordEnd].first) return 0;

  if (compute_value_at_a_given_pointDbg) std::clog << "Entering to the while loop \n";

  while (coordBegin + 1 != coordEnd) {
    if (compute_value_at_a_given_pointDbg) {
      std::clog << "coordBegin : " << coordBegin << "\n";
      std::clog << "coordEnd : " << coordEnd << "\n";
      std::clog << "this->land[level][coordBegin].first : " << this->land[level][coordBegin].first << "\n";
      std::clog << "this->land[level][coordEnd].first : " << this->land[level][coordEnd].first << "\n";
    }

    unsigned newCord = (unsigned)floor((coordEnd + coordBegin) / 2.0);

    if (compute_value_at_a_given_pointDbg) {
      std::clog << "newCord : " << newCord << "\n";
      std::clog << "this->land[level][newCord].first : " << this->land[level][newCord].first << "\n";
      std::cin.ignore();
    }

    if (this->land[level][newCord].first <= x) {
      coordBegin = newCord;
      if (this->land[level][newCord].first == x) return this->land[level][newCord].second;
    } else {
      coordEnd = newCord;
    }
  }

  if (compute_value_at_a_given_pointDbg) {
    std::clog << "x : " << x << " is between : " << this->land[level][coordBegin].first << " a  "
              << this->land[level][coordEnd].first << "\n";
    std::clog << "the y coords are : " << this->land[level][coordBegin].second << " a  "
              << this->land[level][coordEnd].second << "\n";
    std::clog << "coordBegin : " << coordBegin << "\n";
    std::clog << "coordEnd : " << coordEnd << "\n";
    std::cin.ignore();
  }
  return function_value(this->land[level][coordBegin], this->land[level][coordEnd], x);
}

std::ostream& operator<<(std::ostream& out, Persistence_landscape& land) {
  for (size_t level = 0; level != land.land.size(); ++level) {
    out << "Lambda_" << level << ":" << std::endl;
    for (size_t i = 0; i != land.land[level].size(); ++i) {
      if (land.land[level][i].first == -std::numeric_limits<int>::max()) {
        out << "-inf";
      } else {
        if (land.land[level][i].first == std::numeric_limits<int>::max()) {
          out << "+inf";
        } else {
          out << land.land[level][i].first;
        }
      }
      out << " , " << land.land[level][i].second << std::endl;
    }
  }
  return out;
}

void Persistence_landscape::multiply_lanscape_by_real_number_overwrite(double x) {
  for (size_t dim = 0; dim != this->land.size(); ++dim) {
    for (size_t i = 0; i != this->land[dim].size(); ++i) {
      this->land[dim][i].second *= x;
    }
  }
}

bool AbsDbg = false;
Persistence_landscape Persistence_landscape::abs() {
  Persistence_landscape result;
  for (size_t level = 0; level != this->land.size(); ++level) {
    if (AbsDbg) {
      std::clog << "level: " << level << std::endl;
    }
    std::vector<std::pair<double, double> > lambda_n;
    lambda_n.push_back(std::make_pair(-std::numeric_limits<int>::max(), 0));
    for (size_t i = 1; i != this->land[level].size(); ++i) {
      if (AbsDbg) {
        std::clog << "this->land[" << level << "][" << i << "] : " << this->land[level][i].first << " "
                  << this->land[level][i].second << std::endl;
      }
      // if a line segment between this->land[level][i-1] and this->land[level][i] crosses the x-axis, then we have to
      // add one landscape point t o result
      if ((this->land[level][i - 1].second) * (this->land[level][i].second) < 0) {
        double zero =
            find_zero_of_a_line_segment_between_those_two_points(this->land[level][i - 1], this->land[level][i]);

        lambda_n.push_back(std::make_pair(zero, 0));
        lambda_n.push_back(std::make_pair(this->land[level][i].first, fabs(this->land[level][i].second)));
        if (AbsDbg) {
          std::clog << "Adding pair : (" << zero << ",0)" << std::endl;
          std::clog << "In the same step adding pair : (" << this->land[level][i].first << ","
                    << fabs(this->land[level][i].second) << ") " << std::endl;
          std::cin.ignore();
        }
      } else {
        lambda_n.push_back(std::make_pair(this->land[level][i].first, fabs(this->land[level][i].second)));
        if (AbsDbg) {
          std::clog << "Adding pair : (" << this->land[level][i].first << "," << fabs(this->land[level][i].second)
                    << ") " << std::endl;
          std::cin.ignore();
        }
      }
    }
    result.land.push_back(lambda_n);
  }
  return result;
}

Persistence_landscape* Persistence_landscape::new_abs() {
  Persistence_landscape* result = new Persistence_landscape(*this);
  for (size_t level = 0; level != this->land.size(); ++level) {
    if (AbsDbg) {
      std::clog << "level: " << level << std::endl;
    }
    std::vector<std::pair<double, double> > lambda_n;
    lambda_n.push_back(std::make_pair(-std::numeric_limits<int>::max(), 0));
    for (size_t i = 1; i != this->land[level].size(); ++i) {
      if (AbsDbg) {
        std::clog << "this->land[" << level << "][" << i << "] : " << this->land[level][i].first << " "
                  << this->land[level][i].second << std::endl;
      }
      // if a line segment between this->land[level][i-1] and this->land[level][i] crosses the x-axis, then we have to
      // add one landscape point t o result
      if ((this->land[level][i - 1].second) * (this->land[level][i].second) < 0) {
        double zero =
            find_zero_of_a_line_segment_between_those_two_points(this->land[level][i - 1], this->land[level][i]);

        lambda_n.push_back(std::make_pair(zero, 0));
        lambda_n.push_back(std::make_pair(this->land[level][i].first, fabs(this->land[level][i].second)));
        if (AbsDbg) {
          std::clog << "Adding pair : (" << zero << ",0)" << std::endl;
          std::clog << "In the same step adding pair : (" << this->land[level][i].first << ","
                    << fabs(this->land[level][i].second) << ") " << std::endl;
          std::cin.ignore();
        }
      } else {
        lambda_n.push_back(std::make_pair(this->land[level][i].first, fabs(this->land[level][i].second)));
        if (AbsDbg) {
          std::clog << "Adding pair : (" << this->land[level][i].first << "," << fabs(this->land[level][i].second)
                    << ") " << std::endl;
          std::cin.ignore();
        }
      }
    }
    result->land.push_back(lambda_n);
  }
  return result;
}

Persistence_landscape Persistence_landscape::multiply_lanscape_by_real_number_not_overwrite(double x) const {
  std::vector<std::vector<std::pair<double, double> > > result(this->land.size());
  for (size_t dim = 0; dim != this->land.size(); ++dim) {
    std::vector<std::pair<double, double> > lambda_dim(this->land[dim].size());
    for (size_t i = 0; i != this->land[dim].size(); ++i) {
      lambda_dim[i] = std::make_pair(this->land[dim][i].first, x * this->land[dim][i].second);
    }
    result[dim] = lambda_dim;
  }
  Persistence_landscape res;
  // CHANGE
  // res.land = result;
  res.land.swap(result);
  return res;
}  // multiply_lanscape_by_real_number_overwrite

void Persistence_landscape::print_to_file(const char* filename) const {
  std::ofstream write;
  write.open(filename);
  for (size_t dim = 0; dim != this->land.size(); ++dim) {
    write << "#lambda_" << dim << std::endl;
    for (size_t i = 1; i != this->land[dim].size() - 1; ++i) {
      write << this->land[dim][i].first << "  " << this->land[dim][i].second << std::endl;
    }
  }
  write.close();
}

void Persistence_landscape::load_landscape_from_file(const char* filename) {
  bool dbg = false;
  // removing the current content of the persistence landscape.
  this->land.clear();

  // this constructor reads persistence landscape form a file. This file have to be created by this software before head
  std::ifstream in;
  in.open(filename);
  if (!in.good()) {
    std::cerr << "The file : " << filename << " do not exist. The program will now terminate \n";
    throw "The persistence landscape file do not exist. The program will now terminate \n";
  }

  std::string line;
  std::vector<std::pair<double, double> > landscapeAtThisLevel;

  bool isThisAFirsLine = true;
  while (in.good()) {
    getline(in, line);
    if (!(line.length() == 0 || line[0] == '#')) {
      std::stringstream lineSS;
      lineSS << line;
      double beginn, endd;
      lineSS >> beginn;
      lineSS >> endd;
      landscapeAtThisLevel.push_back(std::make_pair(beginn, endd));
      if (dbg) {
        std::clog << "Reading a point : " << beginn << " , " << endd << std::endl;
      }
    } else {
      if (dbg) {
        std::clog << "IGNORE LINE\n";
        getchar();
      }
      if (!isThisAFirsLine) {
        landscapeAtThisLevel.push_back(std::make_pair(std::numeric_limits<int>::max(), 0));
        this->land.push_back(landscapeAtThisLevel);
        std::vector<std::pair<double, double> > newLevelOdLandscape;
        landscapeAtThisLevel.swap(newLevelOdLandscape);
      }
      landscapeAtThisLevel.push_back(std::make_pair(-std::numeric_limits<int>::max(), 0));
      isThisAFirsLine = false;
    }
  }
  if (landscapeAtThisLevel.size() > 1) {
    // seems that the last line of the file is not finished with the newline sign. We need to put what we have in
    // landscapeAtThisLevel to the constructed landscape.
    landscapeAtThisLevel.push_back(std::make_pair(std::numeric_limits<int>::max(), 0));
    this->land.push_back(landscapeAtThisLevel);
  }

  in.close();
}

template <typename T>
Persistence_landscape operation_on_pair_of_landscapes(const Persistence_landscape& land1,
                                                      const Persistence_landscape& land2) {
  bool operation_on_pair_of_landscapesDBG = false;
  if (operation_on_pair_of_landscapesDBG) {
    std::clog << "operation_on_pair_of_landscapes\n";
    std::cin.ignore();
  }
  Persistence_landscape result;
  std::vector<std::vector<std::pair<double, double> > > land(std::max(land1.land.size(), land2.land.size()));
  result.land = land;
  T oper;

  if (operation_on_pair_of_landscapesDBG) {
    for (size_t i = 0; i != std::min(land1.land.size(), land2.land.size()); ++i) {
      std::clog << "land1.land[" << i << "].size() : " << land1.land[i].size() << std::endl;
      std::clog << "land2.land[" << i << "].size() : " << land2.land[i].size() << std::endl;
    }
    getchar();
  }

  for (size_t i = 0; i != std::min(land1.land.size(), land2.land.size()); ++i) {
    std::vector<std::pair<double, double> > lambda_n;
    size_t p = 0;
    size_t q = 0;
    while ((p + 1 < land1.land[i].size()) && (q + 1 < land2.land[i].size())) {
      if (operation_on_pair_of_landscapesDBG) {
        std::clog << "p : " << p << "\n";
        std::clog << "q : " << q << "\n";
        std::clog << "land1.land.size() : " << land1.land.size() << std::endl;
        std::clog << "land2.land.size() : " << land2.land.size() << std::endl;
        std::clog << "land1.land[" << i << "].size() : " << land1.land[i].size() << std::endl;
        std::clog << "land2.land[" << i << "].size() : " << land2.land[i].size() << std::endl;
        std::clog << "land1.land[i][p].first : " << land1.land[i][p].first << "\n";
        std::clog << "land2.land[i][q].first : " << land2.land[i][q].first << "\n";
      }

      if (land1.land[i][p].first < land2.land[i][q].first) {
        if (operation_on_pair_of_landscapesDBG) {
          std::clog << "first \n";
          std::clog << " function_value(land2.land[i][q-1],land2.land[i][q],land1.land[i][p].first) : "
                    << function_value(land2.land[i][q - 1], land2.land[i][q], land1.land[i][p].first) << "\n";
        }
        lambda_n.push_back(
            std::make_pair(land1.land[i][p].first,
                           oper(static_cast<double>(land1.land[i][p].second),
                                function_value(land2.land[i][q - 1], land2.land[i][q], land1.land[i][p].first))));
        ++p;
        continue;
      }
      if (land1.land[i][p].first > land2.land[i][q].first) {
        if (operation_on_pair_of_landscapesDBG) {
          std::clog << "Second \n";
          std::clog << "function_value(" << land1.land[i][p - 1].first << " " << land1.land[i][p - 1].second << " ,"
                    << land1.land[i][p].first << " " << land1.land[i][p].second << ", " << land2.land[i][q].first
                    << " ) : " << function_value(land1.land[i][p - 1], land1.land[i][p - 1], land2.land[i][q].first)
                    << "\n";
          std::clog << "oper( " << function_value(land1.land[i][p], land1.land[i][p - 1], land2.land[i][q].first) << ","
                    << land2.land[i][q].second << " : "
                    << oper(land2.land[i][q].second,
                            function_value(land1.land[i][p], land1.land[i][p - 1], land2.land[i][q].first))
                    << "\n";
        }
        lambda_n.push_back(std::make_pair(
            land2.land[i][q].first, oper(function_value(land1.land[i][p], land1.land[i][p - 1], land2.land[i][q].first),
                                         land2.land[i][q].second)));
        ++q;
        continue;
      }
      if (land1.land[i][p].first == land2.land[i][q].first) {
        if (operation_on_pair_of_landscapesDBG) std::clog << "Third \n";
        lambda_n.push_back(
            std::make_pair(land2.land[i][q].first, oper(land1.land[i][p].second, land2.land[i][q].second)));
        ++p;
        ++q;
      }
      if (operation_on_pair_of_landscapesDBG) {
        std::clog << "Next iteration \n";
      }
    }
    while ((p + 1 < land1.land[i].size()) && (q + 1 >= land2.land[i].size())) {
      if (operation_on_pair_of_landscapesDBG) {
        std::clog << "New point : " << land1.land[i][p].first
                  << "  oper(land1.land[i][p].second,0) : " << oper(land1.land[i][p].second, 0) << std::endl;
      }
      lambda_n.push_back(std::make_pair(land1.land[i][p].first, oper(land1.land[i][p].second, 0)));
      ++p;
    }
    while ((p + 1 >= land1.land[i].size()) && (q + 1 < land2.land[i].size())) {
      if (operation_on_pair_of_landscapesDBG) {
        std::clog << "New point : " << land2.land[i][q].first
                  << " oper(0,land2.land[i][q].second) : " << oper(0, land2.land[i][q].second) << std::endl;
      }
      lambda_n.push_back(std::make_pair(land2.land[i][q].first, oper(0, land2.land[i][q].second)));
      ++q;
    }
    lambda_n.push_back(std::make_pair(std::numeric_limits<int>::max(), 0));
    // CHANGE
    // result.land[i] = lambda_n;
    result.land[i].swap(lambda_n);
  }
  if (land1.land.size() > std::min(land1.land.size(), land2.land.size())) {
    if (operation_on_pair_of_landscapesDBG) {
      std::clog << "land1.land.size() > std::min( land1.land.size() , land2.land.size() )" << std::endl;
    }
    for (size_t i = std::min(land1.land.size(), land2.land.size()); i != std::max(land1.land.size(), land2.land.size());
         ++i) {
      std::vector<std::pair<double, double> > lambda_n(land1.land[i]);
      for (size_t nr = 0; nr != land1.land[i].size(); ++nr) {
        lambda_n[nr] = std::make_pair(land1.land[i][nr].first, oper(land1.land[i][nr].second, 0));
      }
      // CHANGE
      // result.land[i] = lambda_n;
      result.land[i].swap(lambda_n);
    }
  }
  if (land2.land.size() > std::min(land1.land.size(), land2.land.size())) {
    if (operation_on_pair_of_landscapesDBG) {
      std::clog << "( land2.land.size() > std::min( land1.land.size() , land2.land.size() ) ) " << std::endl;
    }
    for (size_t i = std::min(land1.land.size(), land2.land.size()); i != std::max(land1.land.size(), land2.land.size());
         ++i) {
      std::vector<std::pair<double, double> > lambda_n(land2.land[i]);
      for (size_t nr = 0; nr != land2.land[i].size(); ++nr) {
        lambda_n[nr] = std::make_pair(land2.land[i][nr].first, oper(0, land2.land[i][nr].second));
      }
      // CHANGE
      // result.land[i] = lambda_n;
      result.land[i].swap(lambda_n);
    }
  }
  if (operation_on_pair_of_landscapesDBG) {
    std::clog << "operation_on_pair_of_landscapes END\n";
    std::cin.ignore();
  }
  return result;
}  // operation_on_pair_of_landscapes

double compute_maximal_distance_non_symmetric(const Persistence_landscape& pl1, const Persistence_landscape& pl2) {
  bool dbg = false;
  if (dbg) std::clog << " compute_maximal_distance_non_symmetric \n";
  // this distance is not symmetric. It compute ONLY distance between inflection points of pl1 and pl2.
  double maxDist = 0;
  size_t minimalNumberOfLevels = std::min(pl1.land.size(), pl2.land.size());
  for (size_t level = 0; level != minimalNumberOfLevels; ++level) {
    if (dbg) {
      std::clog << "Level : " << level << std::endl;
      std::clog << "PL1 : \n";
      for (size_t i = 0; i != pl1.land[level].size(); ++i) {
        std::clog << "(" << pl1.land[level][i].first << "," << pl1.land[level][i].second << ") \n";
      }
      std::clog << "PL2 : \n";
      for (size_t i = 0; i != pl2.land[level].size(); ++i) {
        std::clog << "(" << pl2.land[level][i].first << "," << pl2.land[level][i].second << ") \n";
      }
      std::cin.ignore();
    }

    int p2Count = 0;
    // In this case, I consider points at the infinity
    for (size_t i = 1; i != pl1.land[level].size() - 1; ++i) {
      while (true) {
        if ((pl1.land[level][i].first >= pl2.land[level][p2Count].first) &&
            (pl1.land[level][i].first <= pl2.land[level][p2Count + 1].first))
          break;
        p2Count++;
      }
      double val =
          fabs(function_value(pl2.land[level][p2Count], pl2.land[level][p2Count + 1], pl1.land[level][i].first) -
               pl1.land[level][i].second);
      if (maxDist <= val) maxDist = val;

      if (dbg) {
        std::clog << pl1.land[level][i].first << "in [" << pl2.land[level][p2Count].first << ","
                  << pl2.land[level][p2Count + 1].first << "] \n";
        std::clog << "pl1[level][i].second : " << pl1.land[level][i].second << std::endl;
        std::clog << "function_value( pl2[level][p2Count] , pl2[level][p2Count+1] , pl1[level][i].first ) : "
                  << function_value(pl2.land[level][p2Count], pl2.land[level][p2Count + 1], pl1.land[level][i].first)
                  << std::endl;
        std::clog << "val : " << val << std::endl;
        std::cin.ignore();
      }
    }
  }

  if (dbg) std::clog << "minimalNumberOfLevels : " << minimalNumberOfLevels << std::endl;

  if (minimalNumberOfLevels < pl1.land.size()) {
    for (size_t level = minimalNumberOfLevels; level != pl1.land.size(); ++level) {
      for (size_t i = 0; i != pl1.land[level].size(); ++i) {
        if (dbg) std::clog << "pl1[level][i].second  : " << pl1.land[level][i].second << std::endl;
        if (maxDist < pl1.land[level][i].second) maxDist = pl1.land[level][i].second;
      }
    }
  }
  return maxDist;
}

double compute_distance_of_landscapes(const Persistence_landscape& first, const Persistence_landscape& second,
                                      double p) {
  bool dbg = false;
  // This is what we want to compute: (\int_{- \infty}^{+\infty}| first-second |^p)^(1/p). We will do it one step at a
  // time:

  // first-second :
  Persistence_landscape lan = first - second;

  //| first-second |:
  lan = lan.abs();

  if (dbg) {
    std::clog << "Abs of difference ; " << lan << std::endl;
    getchar();
  }

  if (p < std::numeric_limits<double>::max()) {
    // \int_{- \infty}^{+\infty}| first-second |^p
    double result;
    if (p != 1) {
      if (dbg) std::clog << "Power != 1, compute integral to the power p\n";
      result = lan.compute_integral_of_landscape(p);
    } else {
      if (dbg) std::clog << "Power = 1, compute integral \n";
      result = lan.compute_integral_of_landscape();
    }
    // (\int_{- \infty}^{+\infty}| first-second |^p)^(1/p)
    return pow(result, 1.0 / p);
  } else {
    // p == infty
    if (dbg) std::clog << "Power = infty, compute maximum \n";
    return lan.compute_maximum();
  }
}

double compute_max_norm_distance_of_landscapes(const Persistence_landscape& first,
                                               const Persistence_landscape& second) {
  return std::max(compute_maximal_distance_non_symmetric(first, second),
                  compute_maximal_distance_non_symmetric(second, first));
}

bool comparePairsForMerging(std::pair<double, unsigned> first, std::pair<double, unsigned> second) {
  return (first.first < second.first);
}

double compute_inner_product(const Persistence_landscape& l1, const Persistence_landscape& l2) {
  bool dbg = false;
  double result = 0;

  for (size_t level = 0; level != std::min(l1.size(), l2.size()); ++level) {
    if (dbg) {
      std::clog << "Computing inner product for a level : " << level << std::endl;
      getchar();
    }
    auto&& l1_land_level = l1.land[level];
    auto&& l2_land_level = l2.land[level];

    if (l1_land_level.size() * l2_land_level.size() == 0) continue;

    // endpoints of the interval on which we will compute the inner product of two locally linear functions:
    double x1 = -std::numeric_limits<int>::max();
    double x2;
    if (l1_land_level[1].first < l2_land_level[1].first) {
      x2 = l1_land_level[1].first;
    } else {
      x2 = l2_land_level[1].first;
    }

    // iterators for the landscapes l1 and l2
    size_t l1It = 0;
    size_t l2It = 0;

    while ((l1It < l1_land_level.size() - 1) && (l2It < l2_land_level.size() - 1)) {
      // compute the value of a inner product on a interval [x1,x2]

      double a, b, c, d;

      if (l1_land_level[l1It + 1].first != l1_land_level[l1It].first) {
        a = (l1_land_level[l1It + 1].second - l1_land_level[l1It].second) /
            (l1_land_level[l1It + 1].first - l1_land_level[l1It].first);
      } else {
        a = 0;
      }
      b = l1_land_level[l1It].second - a * l1_land_level[l1It].first;
      if (l2_land_level[l2It + 1].first != l2_land_level[l2It].first) {
        c = (l2_land_level[l2It + 1].second - l2_land_level[l2It].second) /
            (l2_land_level[l2It + 1].first - l2_land_level[l2It].first);
      } else {
        c = 0;
      }
      d = l2_land_level[l2It].second - c * l2_land_level[l2It].first;

      double contributionFromThisPart = (a * c * x2 * x2 * x2 / 3 + (a * d + b * c) * x2 * x2 / 2 + b * d * x2) -
                                        (a * c * x1 * x1 * x1 / 3 + (a * d + b * c) * x1 * x1 / 2 + b * d * x1);

      result += contributionFromThisPart;

      if (dbg) {
        std::clog << "[l1_land_level[l1It].first,l1_land_level[l1It+1].first] : " << l1_land_level[l1It].first
                  << " , " << l1_land_level[l1It + 1].first << std::endl;
        std::clog << "[l2_land_level[l2It].first,l2_land_level[l2It+1].first] : " << l2_land_level[l2It].first
                  << " , " << l2_land_level[l2It + 1].first << std::endl;
        std::clog << "a : " << a << ", b : " << b << " , c: " << c << ", d : " << d << std::endl;
        std::clog << "x1 : " << x1 << " , x2 : " << x2 << std::endl;
        std::clog << "contributionFromThisPart : " << contributionFromThisPart << std::endl;
        std::clog << "result : " << result << std::endl;
        getchar();
      }

      // we have two intervals in which functions are constant:
      // [l1_land_level[l1It].first , l1_land_level[l1It+1].first]
      // and
      // [l2_land_level[l2It].first , l2_land_level[l2It+1].first]
      // We also have an interval [x1,x2]. Since the intervals in the landscapes cover the whole R, then it is clear
      // that x2
      // is either l1_land_level[l1It+1].first of l2_land_level[l2It+1].first or both. Lets test it.
      if (x2 == l1_land_level[l1It + 1].first) {
        if (x2 == l2_land_level[l2It + 1].first) {
          // in this case, we increment both:
          ++l2It;
          if (dbg) {
            std::clog << "Incrementing both \n";
          }
        } else {
          if (dbg) {
            std::clog << "Incrementing first \n";
          }
        }
        ++l1It;
      } else {
        // in this case we increment l2It
        ++l2It;
        if (dbg) {
          std::clog << "Incrementing second \n";
        }
      }

      if ( l1It + 1 >= l1_land_level.size()  )break;
      if ( l2It + 1 >= l2_land_level.size()  )break;

      // Now, we shift x1 and x2:
      x1 = x2;
      if (l1_land_level[l1It + 1].first < l2_land_level[l2It + 1].first) {
        x2 = l1_land_level[l1It + 1].first;
      } else {
        x2 = l2_land_level[l2It + 1].first;
      }
    }
  }
  return result;
}

void Persistence_landscape::plot(const char* filename, double xRangeBegin, double xRangeEnd, double yRangeBegin,
                                 double yRangeEnd, int from, int to) {
  // this program create a gnuplot script file that allows to plot persistence diagram.
  std::ofstream out;

  std::ostringstream gnuplot_script;
  gnuplot_script << filename << "_GnuplotScript";
  out.open(gnuplot_script.str().c_str());

  if ((xRangeBegin != std::numeric_limits<double>::max()) || (xRangeEnd != std::numeric_limits<double>::max()) ||
      (yRangeBegin != std::numeric_limits<double>::max()) || (yRangeEnd != std::numeric_limits<double>::max())) {
    out << "set xrange [" << xRangeBegin << " : " << xRangeEnd << "]" << std::endl;
    out << "set yrange [" << yRangeBegin << " : " << yRangeEnd << "]" << std::endl;
  }

  if (from == std::numeric_limits<int>::max()) {
    from = 0;
  }
  if (to == std::numeric_limits<int>::max()) {
    to = this->land.size();
  }

  out << "plot ";
  for (size_t lambda = std::min((size_t)from, this->land.size()); lambda != std::min((size_t)to, this->land.size());
       ++lambda) {
    // out << "     '-' using 1:2 title 'l" << lambda << "' with lp";
    out << "     '-' using 1:2 notitle with lp";
    if (lambda + 1 != std::min((size_t)to, this->land.size())) {
      out << ", \\";
    }
    out << std::endl;
  }

  for (size_t lambda = std::min((size_t)from, this->land.size()); lambda != std::min((size_t)to, this->land.size());
       ++lambda) {
    for (size_t i = 1; i != this->land[lambda].size() - 1; ++i) {
      out << this->land[lambda][i].first << " " << this->land[lambda][i].second << std::endl;
    }
    out << "EOF" << std::endl;
  }
  std::clog << "To visualize, install gnuplot and type the command: gnuplot -persist -e \"load \'"
            << gnuplot_script.str().c_str() << "\'\"" << std::endl;
}

}  // namespace Persistence_representations
}  // namespace Gudhi

#endif  // PERSISTENCE_LANDSCAPE_H_
