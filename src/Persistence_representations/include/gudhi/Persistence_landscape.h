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

#ifndef PERSISTENCE_LANDSCAPE_H_
#define PERSISTENCE_LANDSCAPE_H_

// standard include
#ifdef DEBUG_TRACES
#include <iostream>   // std::cerr, std::clog
#endif
#include <ostream>    // std::ostream
#include <fstream>    // std::ofstream, std::ifstream
#include <sstream>    // std::stringstream, std::ostringstream
#include <stdexcept>  // std::invalid_argument
#include <cstddef>    // std::size_t
#include <limits>     // std::numeric_limits
#include <algorithm>  // std::sort
#include <cmath>      // std::min, std::max, std::pow, std::fabs
#include <utility>    // std::pair
#include <string>
#include <vector>

// gudhi include
#include <gudhi/read_persistence_from_file.h>
#include <gudhi/common_persistence_representations.h>
#include <gudhi/Debug_utils.h>

namespace Gudhi {
namespace Persistence_representations {

// pre declaration needed before C++20 for friends with templates defined inside a class
class Persistence_landscape;
template <typename Operation>
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
 **/
class Persistence_landscape
{
 public:
  /**
   * Default constructor.
   **/
  Persistence_landscape() { this->_set_up_numbers_of_functions_for_vectorization_and_projections_to_reals(); }

  /**
   * Constructor that takes as an input a vector of birth-death pairs.
   **/
  Persistence_landscape(const std::vector<std::pair<double, double> >& p,
                        std::size_t number_of_levels = std::numeric_limits<std::size_t>::max())
  {
    this->_construct_persistence_landscape_from_barcode(p, number_of_levels);
    this->_set_up_numbers_of_functions_for_vectorization_and_projections_to_reals();
  }

  /**
   * Constructor that reads persistence intervals from file and creates persistence landscape. The format of the
   * input file is the following: in each line we put birth-death pair. Last line is assumed to be empty.
   * Even if the points within a line are not ordered, they will be ordered while the input is read.
   **/
  Persistence_landscape(const char* filename,
                        unsigned int dimension = std::numeric_limits<unsigned int>::max(),
                        std::size_t number_of_levels = std::numeric_limits<std::size_t>::max())
  {
    std::vector<std::pair<double, double> > barcode;
    if (dimension < std::numeric_limits<unsigned int>::max()) {
      barcode = read_persistence_intervals_in_one_dimension_from_file(filename, dimension);
    } else {
      barcode = read_persistence_intervals_in_one_dimension_from_file(filename);
    }
    this->_construct_persistence_landscape_from_barcode(barcode, number_of_levels);
    this->_set_up_numbers_of_functions_for_vectorization_and_projections_to_reals();
  }

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
   * functions).
   **/
  double compute_integral_of_landscape() const;

  /**
   * This function compute integral of the 'level'-level of a landscape.
   **/
  double compute_integral_of_a_level_of_a_landscape(std::size_t level) const;

  /**
   * This function compute integral of the landscape p-th power of a landscape (defined formally as sum of integrals
   * on R of p-th powers of all landscape functions).
   **/
  double compute_integral_of_landscape(double p) const;  // this function compute integral of p-th power of landscape.

  /**
   * A function that computes the value of a landscape at a given point.
   * The parameters of the function are: unsigned level and double x.
   * The procedure will compute the value of the level-landscape at the point x.
   **/
  double compute_value_at_a_given_point(unsigned int level, double x) const;

  /**
   * Writing landscape into a stream. A i-th level landscape starts with a string "lambda_i".
   * Then the discontinuity points of the landscapes follows.
   * Shall those points be joined with lines, we will obtain the i-th landscape function.
   **/
  friend std::ostream& operator<<(std::ostream& out, const Persistence_landscape& land)
  {
    for (std::size_t level = 0; level != land.land_.size(); ++level) {
      out << "Lambda_" << level << ":" << std::endl;
      for (std::size_t i = 0; i != land.land_[level].size(); ++i) {
        if (land.land_[level][i].first == -std::numeric_limits<int>::max()) {
          out << "-inf";
        } else {
          if (land.land_[level][i].first == std::numeric_limits<int>::max()) {
            out << "+inf";
          } else {
            out << land.land_[level][i].first;
          }
        }
        out << " , " << land.land_[level][i].second << std::endl;
      }
    }
    return out;
  }

  template <typename Operation>
  friend Persistence_landscape operation_on_pair_of_landscapes(const Persistence_landscape& land1,
                                                               const Persistence_landscape& land2)
  {
#ifdef DEBUG_TRACES
    std::clog << "operation_on_pair_of_landscapes\n";
#endif
    Persistence_landscape result;
    std::vector<std::vector<std::pair<double, double> > > land(std::max(land1.land_.size(), land2.land_.size()));
    result.land_ = land;
    Operation oper;

#ifdef DEBUG_TRACES
    for (std::size_t i = 0; i != std::min(land1.land_.size(), land2.land_.size()); ++i) {
      std::clog << "land1.land[" << i << "].size() : " << land1.land_[i].size() << std::endl;
      std::clog << "land2.land[" << i << "].size() : " << land2.land_[i].size() << std::endl;
    }
#endif

    for (std::size_t i = 0; i != std::min(land1.land_.size(), land2.land_.size()); ++i) {
      std::vector<std::pair<double, double> > lambda_n;
      std::size_t p = 0;
      std::size_t q = 0;
      while ((p + 1 < land1.land_[i].size()) && (q + 1 < land2.land_[i].size())) {
#ifdef DEBUG_TRACES
        std::clog << "p : " << p << "\n";
        std::clog << "q : " << q << "\n";
        std::clog << "land1.land.size() : " << land1.land_.size() << std::endl;
        std::clog << "land2.land.size() : " << land2.land_.size() << std::endl;
        std::clog << "land1.land[" << i << "].size() : " << land1.land_[i].size() << std::endl;
        std::clog << "land2.land[" << i << "].size() : " << land2.land_[i].size() << std::endl;
        std::clog << "land1.land[i][p].first : " << land1.land_[i][p].first << "\n";
        std::clog << "land2.land[i][q].first : " << land2.land_[i][q].first << "\n";
#endif

        if (land1.land_[i][p].first < land2.land_[i][q].first) {
#ifdef DEBUG_TRACES
          std::clog << "first \n";
          std::clog << " function_value(land2.land[i][q-1],land2.land[i][q],land1.land[i][p].first) : "
                    << function_value(land2.land_[i][q - 1], land2.land_[i][q], land1.land_[i][p].first) << "\n";
#endif
          lambda_n.push_back(
              std::make_pair(land1.land_[i][p].first,
                             oper(static_cast<double>(land1.land_[i][p].second),
                                  function_value(land2.land_[i][q - 1], land2.land_[i][q], land1.land_[i][p].first))));
          ++p;
          continue;
        }
        if (land1.land_[i][p].first > land2.land_[i][q].first) {
#ifdef DEBUG_TRACES
          std::clog << "Second \n";
          std::clog << "function_value(" << land1.land_[i][p - 1].first << " " << land1.land_[i][p - 1].second << " ,"
                    << land1.land_[i][p].first << " " << land1.land_[i][p].second << ", " << land2.land_[i][q].first
                    << " ) : " << function_value(land1.land_[i][p - 1], land1.land_[i][p - 1], land2.land_[i][q].first)
                    << "\n";
          std::clog << "oper( " << function_value(land1.land_[i][p], land1.land_[i][p - 1], land2.land_[i][q].first)
                    << "," << land2.land_[i][q].second << " : "
                    << oper(land2.land_[i][q].second,
                            function_value(land1.land_[i][p], land1.land_[i][p - 1], land2.land_[i][q].first))
                    << "\n";
#endif
          lambda_n.push_back(
              std::make_pair(land2.land_[i][q].first,
                             oper(function_value(land1.land_[i][p], land1.land_[i][p - 1], land2.land_[i][q].first),
                                  land2.land_[i][q].second)));
          ++q;
          continue;
        }
        if (land1.land_[i][p].first == land2.land_[i][q].first) {
#ifdef DEBUG_TRACES
          std::clog << "Third \n";
#endif
          lambda_n.push_back(
              std::make_pair(land2.land_[i][q].first, oper(land1.land_[i][p].second, land2.land_[i][q].second)));
          ++p;
          ++q;
        }
#ifdef DEBUG_TRACES
        std::clog << "Next iteration \n";
#endif
      }
      while ((p + 1 < land1.land_[i].size()) && (q + 1 >= land2.land_[i].size())) {
#ifdef DEBUG_TRACES
        std::clog << "New point : " << land1.land_[i][p].first
                  << "  oper(land1.land[i][p].second,0) : " << oper(land1.land_[i][p].second, 0) << std::endl;
#endif
        lambda_n.push_back(std::make_pair(land1.land_[i][p].first, oper(land1.land_[i][p].second, 0)));
        ++p;
      }
      while ((p + 1 >= land1.land_[i].size()) && (q + 1 < land2.land_[i].size())) {
#ifdef DEBUG_TRACES
        std::clog << "New point : " << land2.land_[i][q].first
                  << " oper(0,land2.land[i][q].second) : " << oper(0, land2.land_[i][q].second) << std::endl;
#endif
        lambda_n.push_back(std::make_pair(land2.land_[i][q].first, oper(0, land2.land_[i][q].second)));
        ++q;
      }
      lambda_n.push_back(std::make_pair(std::numeric_limits<int>::max(), 0));
      // CHANGE
      // result.land[i] = lambda_n;
      result.land_[i].swap(lambda_n);
    }
    if (land1.land_.size() > std::min(land1.land_.size(), land2.land_.size())) {
#ifdef DEBUG_TRACES
      std::clog << "land1.land.size() > std::min( land1.land.size() , land2.land.size() )" << std::endl;
#endif
      for (std::size_t i = std::min(land1.land_.size(), land2.land_.size());
           i != std::max(land1.land_.size(), land2.land_.size());
           ++i) {
        std::vector<std::pair<double, double> > lambda_n(land1.land_[i]);
        for (std::size_t nr = 0; nr != land1.land_[i].size(); ++nr) {
          lambda_n[nr] = std::make_pair(land1.land_[i][nr].first, oper(land1.land_[i][nr].second, 0));
        }
        // CHANGE
        // result.land[i] = lambda_n;
        result.land_[i].swap(lambda_n);
      }
    }
    if (land2.land_.size() > std::min(land1.land_.size(), land2.land_.size())) {
#ifdef DEBUG_TRACES
      std::clog << "( land2.land.size() > std::min( land1.land.size() , land2.land.size() ) ) " << std::endl;
#endif
      for (std::size_t i = std::min(land1.land_.size(), land2.land_.size());
           i != std::max(land1.land_.size(), land2.land_.size());
           ++i) {
        std::vector<std::pair<double, double> > lambda_n(land2.land_[i]);
        for (std::size_t nr = 0; nr != land2.land_[i].size(); ++nr) {
          lambda_n[nr] = std::make_pair(land2.land_[i][nr].first, oper(0, land2.land_[i][nr].second));
        }
        // CHANGE
        // result.land[i] = lambda_n;
        result.land_[i].swap(lambda_n);
      }
    }
#ifdef DEBUG_TRACES
    std::clog << "operation_on_pair_of_landscapes END\n";
#endif
    return result;
  }  // operation_on_pair_of_landscapes

  // TODO: `add_two_landscapes` and `subtract_two_landscapes` do not seem to be used outside of resp.
  // `operator+` and `operator-` and their doc is marked private. Are they really necessary?
  // Can `operation_on_pair_of_landscapes` not be directly called in `operator+` and `operator-`?

  /**
   * \private A function that compute sum of two landscapes.
   **/
  friend Persistence_landscape add_two_landscapes(const Persistence_landscape& land1,
                                                  const Persistence_landscape& land2)
  {
    return operation_on_pair_of_landscapes<std::plus<double> >(land1, land2);
  }

  /**
   * \private A function that compute difference of two landscapes.
   **/
  friend Persistence_landscape subtract_two_landscapes(const Persistence_landscape& land1,
                                                       const Persistence_landscape& land2)
  {
    return operation_on_pair_of_landscapes<std::minus<double> >(land1, land2);
  }

  /**
   * An operator +, that compute sum of two landscapes.
   **/
  friend Persistence_landscape operator+(const Persistence_landscape& first, const Persistence_landscape& second)
  {
    return add_two_landscapes(first, second);
  }

  /**
   * An operator -, that compute difference of two landscapes.
   **/
  friend Persistence_landscape operator-(const Persistence_landscape& first, const Persistence_landscape& second)
  {
    return subtract_two_landscapes(first, second);
  }

  /**
   * An operator * that allows multiplication of a landscape by a real number.
   **/
  friend Persistence_landscape operator*(const Persistence_landscape& first, double con)
  {
    return first._multiply_landscape_by_real_number_not_overwrite(con);
  }

  /**
   * An operator * that allows multiplication of a landscape by a real number (order of parameters swapped).
   **/
  friend Persistence_landscape operator*(double con, const Persistence_landscape& first)
  {
    return first._multiply_landscape_by_real_number_not_overwrite(con);
  }

  /**
   * Operator +=. The second parameter is persistence landscape.
   **/
  Persistence_landscape operator+=(const Persistence_landscape& rhs)
  {
    *this = *this + rhs;
    return *this;
  }

  /**
   * Operator -=. The second parameter is a persistence landscape.
   **/
  Persistence_landscape operator-=(const Persistence_landscape& rhs)
  {
    *this = *this - rhs;
    return *this;
  }

  /**
   * Operator *=. The second parameter is a real number by which the y values of all landscape functions are multiplied.
   * The x-values remain unchanged.
   **/
  Persistence_landscape operator*=(double x)
  {
    *this = *this * x;
    return *this;
  }

  /**
   * Operator /=. The second parameter is a real number.
   **/
  Persistence_landscape operator/=(double x)
  {
    if (x == 0) throw std::invalid_argument("In operator /=, division by 0.");
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
  double compute_maximum() const
  {
    double maxValue = 0;
    if (this->land_.size()) {
      maxValue = -std::numeric_limits<int>::max();
      for (std::size_t i = 0; i != this->land_[0].size(); ++i) {
        if (this->land_[0][i].second > maxValue) maxValue = this->land_[0][i].second;
      }
    }
    return maxValue;
  }

  /**
   * \private Computations of minimum (y) value of landscape.
   **/
  double compute_minimum() const
  {
    double minValue = 0;
    if (this->land_.size()) {
      minValue = std::numeric_limits<int>::max();
      for (std::size_t i = 0; i != this->land_[0].size(); ++i) {
        if (this->land_[0][i].second < minValue) minValue = this->land_[0][i].second;
      }
    }
    return minValue;
  }

  /**
   * \private Computations of a \f$L^i\f$ norm of landscape, where i is the input parameter.
   **/
  double compute_norm_of_landscape(double i)
  {
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
  double operator()(unsigned int level, double x) const { return this->compute_value_at_a_given_point(level, x); }

  /**
   * \private Computations of \f$L^{\infty}\f$ distance between two landscapes.
   **/
  friend double compute_max_norm_distance_of_landscapes(const Persistence_landscape& first,
                                                        const Persistence_landscape& second)
  {
    return std::max(compute_maximal_distance_non_symmetric(first, second),
                    compute_maximal_distance_non_symmetric(second, first));
  }

  /**
   * \private Computations of \f$L^{p}\f$ distance between two landscapes. p is the parameter of the procedure.
   **/
  friend double compute_distance_of_landscapes(const Persistence_landscape& first,
                                               const Persistence_landscape& second,
                                               double p)
  {
    // This is what we want to compute: (\int_{- \infty}^{+\infty}| first-second |^p)^(1/p). We will do it one step at a
    // time:

    // first-second :
    Persistence_landscape lan = first - second;

    //| first-second |:
    lan = lan.abs();

#ifdef DEBUG_TRACES
    std::clog << "Abs of difference ; " << lan << std::endl;
#endif

    if (p < std::numeric_limits<double>::max()) {
      // \int_{- \infty}^{+\infty}| first-second |^p
      double result;
      if (p != 1) {
#ifdef DEBUG_TRACES
        std::clog << "Power != 1, compute integral to the power p\n";
#endif
        result = lan.compute_integral_of_landscape(p);
      } else {
#ifdef DEBUG_TRACES
        std::clog << "Power = 1, compute integral \n";
#endif
        result = lan.compute_integral_of_landscape();
      }
      // (\int_{- \infty}^{+\infty}| first-second |^p)^(1/p)
      return std::pow(result, 1.0 / p);
    } else {
      // p == infty
#ifdef DEBUG_TRACES
      std::clog << "Power = infty, compute maximum \n";
#endif
      return lan.compute_maximum();
    }
  }

  /**
   * Function to compute absolute value of a PL function. The representation of persistence landscapes allow to store
   * general PL-function. When computing distance between two landscapes, we compute difference between
   * them. In this case, a general PL-function with negative value can appear as a result. Then in order to compute
   * distance, we need to take its absolute value. This is the purpose of this procedure.
   **/
  Persistence_landscape abs();

  // TODO: it is public but not documented.
  Persistence_landscape* new_abs();

  /**
   * Computes the number of landscape functions.
   **/
  std::size_t size() const { return this->land_.size(); }

  /**
   * Compute maximal value of lambda-level landscape.
   **/
  double find_max(unsigned int lambda) const;

  /**
   * \private Function to compute inner (scalar) product of two landscapes.
   **/
  friend double compute_inner_product(const Persistence_landscape& l1, const Persistence_landscape& l2)
  {
    double result = 0;

    for (std::size_t level = 0; level != std::min(l1.size(), l2.size()); ++level) {
#ifdef DEBUG_TRACES
      std::clog << "Computing inner product for a level : " << level << std::endl;
#endif
      auto&& l1_land_level = l1.land_[level];
      auto&& l2_land_level = l2.land_[level];

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
      std::size_t l1It = 0;
      std::size_t l2It = 0;

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

#ifdef DEBUG_TRACES
        std::clog << "[l1_land_level[l1It].first,l1_land_level[l1It+1].first] : " << l1_land_level[l1It].first << " , "
                  << l1_land_level[l1It + 1].first << std::endl;
        std::clog << "[l2_land_level[l2It].first,l2_land_level[l2It+1].first] : " << l2_land_level[l2It].first << " , "
                  << l2_land_level[l2It + 1].first << std::endl;
        std::clog << "a : " << a << ", b : " << b << " , c: " << c << ", d : " << d << std::endl;
        std::clog << "x1 : " << x1 << " , x2 : " << x2 << std::endl;
        std::clog << "contributionFromThisPart : " << contributionFromThisPart << std::endl;
        std::clog << "result : " << result << std::endl;
#endif

        // we have two intervals in which functions are constant:
        // [l1_land_level[l1It].first , l1_land_level[l1It+1].first]
        // and
        // [l2_land_level[l2It].first , l2_land_level[l2It+1].first]
        // We also have an interval [x1,x2]. Since the intervals in the landscapes cover the whole R, then it is clear
        // that x2 is either l1_land_level[l1It+1].first of l2_land_level[l2It+1].first or both. Lets test it.
        if (x2 == l1_land_level[l1It + 1].first) {
          if (x2 == l2_land_level[l2It + 1].first) {
            // in this case, we increment both:
            ++l2It;
#ifdef DEBUG_TRACES
            std::clog << "Incrementing both \n";
#endif
          } else {
#ifdef DEBUG_TRACES
            std::clog << "Incrementing first \n";
#endif
          }
          ++l1It;
        } else {
          // in this case we increment l2It
          ++l2It;
#ifdef DEBUG_TRACES
          std::clog << "Incrementing second \n";
#endif
        }

        if (l1It + 1 >= l1_land_level.size()) break;
        if (l2It + 1 >= l2_land_level.size()) break;

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

  // Implementations of functions for various concepts.

  /**
   * The number of projections to R is defined to the number of nonzero landscape functions. I-th projection is an
   * integral of i-th landscape function over whole R.
   * This function is required by the Real_valued_topological_data concept.
   * At the moment this function is not tested, since it is quite likely to be changed in the future. Given this, when
   * using it, keep in mind that it will be most likely changed in the next versions.
   **/
  double project_to_R(int number_of_function) const
  {
    return this->compute_integral_of_a_level_of_a_landscape((std::size_t)number_of_function);
  }

  /**
   * The function gives the number of possible projections to R. This function is required by the
   * Real_valued_topological_data concept.
   **/
  std::size_t number_of_projections_to_R() const { return this->number_of_functions_for_projections_to_reals_; }

  /**
   * This function produce a vector of doubles based on a landscape. It is required in a concept
   * Vectorized_topological_data
   */
  std::vector<double> vectorize(int number_of_function) const
  {
    // TODO(PD) think of something smarter over here
    std::vector<double> v;
    if (static_cast<std::size_t>(number_of_function) > this->land_.size()) {
      return v;
    }
    v.reserve(this->land_[number_of_function].size());
    for (std::size_t i = 0; i != this->land_[number_of_function].size(); ++i) {
      v.push_back(this->land_[number_of_function][i].second);
    }
    return v;
  }

  /**
   * This function return the number of functions that allows vectorization of persistence landscape. It is required in
   * a concept Vectorized_topological_data.
   **/
  std::size_t number_of_vectorize_functions() const { return this->number_of_functions_for_vectorization_; }

  /**
   * A function to compute averaged persistence landscape, based on vector of persistence landscapes.
   * This function is required by Topological_data_with_averages concept.
   **/
  void compute_average(const std::vector<Persistence_landscape*>& to_average)
  {
#ifdef DEBUG_TRACES
    std::clog << "to_average.size() : " << to_average.size() << std::endl;
#endif

    std::vector<Persistence_landscape*> nextLevelMerge(to_average.size());
    for (std::size_t i = 0; i != to_average.size(); ++i) {
      nextLevelMerge[i] = to_average[i];
    }

    // In the loop, we will create dynamically a number of intermediate complexes. We have to clean that up, but
    // we cannot erase the initial landscapes we have to average. In this case, we simply check if the nextLevelMerge
    // are the input landscapes or the ones created in that loop by using this extra variable.
    bool is_this_first_level = true;

    while (nextLevelMerge.size() != 1) {
#ifdef DEBUG_TRACES
      std::clog << "nextLevelMerge.size() : " << nextLevelMerge.size() << std::endl;
#endif
      std::vector<Persistence_landscape*> nextNextLevelMerge;
      nextNextLevelMerge.reserve(to_average.size());
      for (std::size_t i = 0; i < nextLevelMerge.size(); i = i + 2) {
#ifdef DEBUG_TRACES
        std::clog << "i : " << i << std::endl;
#endif
        Persistence_landscape* l = new Persistence_landscape;
        if (i + 1 != nextLevelMerge.size()) {
          (*l) = (*nextLevelMerge[i]) + (*nextLevelMerge[i + 1]);
        } else {
          (*l) = *nextLevelMerge[i];
        }
        nextNextLevelMerge.push_back(l);
      }
#ifdef DEBUG_TRACES
      std::clog << "After this iteration \n";
#endif

      if (!is_this_first_level) {
        // deallocate the memory if the vector nextLevelMerge do not consist of the initial landscapes
        for (std::size_t i = 0; i != nextLevelMerge.size(); ++i) {
          delete nextLevelMerge[i];
        }
      }
      is_this_first_level = false;
      nextLevelMerge.swap(nextNextLevelMerge);
    }
    (*this) = (*nextLevelMerge[0]);
    if (!is_this_first_level) delete nextLevelMerge[0];
    (*this) *= 1 / static_cast<double>(to_average.size());
  }

  /**
   * A function to compute distance between persistence landscape.
   * The parameter of this function is a Persistence_landscape.
   * This function is required in Topological_data_with_distances concept.
   * For max norm distance, set power to std::numeric_limits<double>::max()
   **/
  double distance(const Persistence_landscape& second, double power = 1) const
  {
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
  double compute_scalar_product(const Persistence_landscape& second) const
  {
    return compute_inner_product((*this), second);
  }

  // end of implementation of functions needed for concepts.

  /**
   * This procedure returns y-range of a given level persistence landscape. If a default value is used, the y-range
   * of 0th level landscape is given (and this range contains the ranges of all other landscapes).
   **/
  std::pair<double, double> get_y_range(std::size_t level = 0) const
  {
    std::pair<double, double> result;
    if (level < this->land_.size()) {
      double maxx = this->compute_maximum();
      double minn = this->compute_minimum();
      result = std::make_pair(minn, maxx);
    } else {
      result = std::make_pair(0, 0);
    }
    return result;
  }

  // a function used to create a gnuplot script for visualization of landscapes
  void plot(const char* filename,
            double xRangeBegin = std::numeric_limits<double>::max(),
            double xRangeEnd = std::numeric_limits<double>::max(),
            double yRangeBegin = std::numeric_limits<double>::max(),
            double yRangeEnd = std::numeric_limits<double>::max(),
            int from = std::numeric_limits<int>::max(),
            int to = std::numeric_limits<int>::max());

 private:
  friend double compute_maximal_distance_non_symmetric(const Persistence_landscape& pl1,
                                                       const Persistence_landscape& pl2)
  {
#ifdef DEBUG_TRACES
    std::clog << " compute_maximal_distance_non_symmetric \n";
#endif
    // this distance is not symmetric. It compute ONLY distance between inflection points of pl1 and pl2.
    double maxDist = 0;
    std::size_t minimalNumberOfLevels = std::min(pl1.land_.size(), pl2.land_.size());
    for (std::size_t level = 0; level != minimalNumberOfLevels; ++level) {
#ifdef DEBUG_TRACES
      std::clog << "Level : " << level << std::endl;
      std::clog << "PL1 : \n";
      for (std::size_t i = 0; i != pl1.land_[level].size(); ++i) {
        std::clog << "(" << pl1.land_[level][i].first << "," << pl1.land_[level][i].second << ") \n";
      }
      std::clog << "PL2 : \n";
      for (std::size_t i = 0; i != pl2.land_[level].size(); ++i) {
        std::clog << "(" << pl2.land_[level][i].first << "," << pl2.land_[level][i].second << ") \n";
      }
#endif

      int p2Count = 0;
      // In this case, I consider points at the infinity
      for (std::size_t i = 1; i != pl1.land_[level].size() - 1; ++i) {
        while (true) {
          if ((pl1.land_[level][i].first >= pl2.land_[level][p2Count].first) &&
              (pl1.land_[level][i].first <= pl2.land_[level][p2Count + 1].first))
            break;
          p2Count++;
        }
        double val = std::fabs(
            function_value(pl2.land_[level][p2Count], pl2.land_[level][p2Count + 1], pl1.land_[level][i].first) -
            pl1.land_[level][i].second);
        if (maxDist <= val) maxDist = val;

#ifdef DEBUG_TRACES
        std::clog << pl1.land_[level][i].first << "in [" << pl2.land_[level][p2Count].first << ","
                  << pl2.land_[level][p2Count + 1].first << "] \n";
        std::clog << "pl1[level][i].second : " << pl1.land_[level][i].second << std::endl;
        std::clog << "function_value( pl2[level][p2Count] , pl2[level][p2Count+1] , pl1[level][i].first ) : "
                  << function_value(pl2.land_[level][p2Count], pl2.land_[level][p2Count + 1], pl1.land_[level][i].first)
                  << std::endl;
        std::clog << "val : " << val << std::endl;
#endif
      }
    }

#ifdef DEBUG_TRACES
    std::clog << "minimalNumberOfLevels : " << minimalNumberOfLevels << std::endl;
#endif

    if (minimalNumberOfLevels < pl1.land_.size()) {
      for (std::size_t level = minimalNumberOfLevels; level != pl1.land_.size(); ++level) {
        for (std::size_t i = 0; i != pl1.land_[level].size(); ++i) {
#ifdef DEBUG_TRACES
          std::clog << "pl1[level][i].second  : " << pl1.land_[level][i].second << std::endl;
#endif
          if (maxDist < pl1.land_[level][i].second) maxDist = pl1.land_[level][i].second;
        }
      }
    }
    return maxDist;
  }

  std::vector<std::vector<std::pair<double, double> > > land_;
  std::size_t number_of_functions_for_vectorization_;
  std::size_t number_of_functions_for_projections_to_reals_;

  void _construct_persistence_landscape_from_barcode(
      const std::vector<std::pair<double, double> >& p,
      std::size_t number_of_levels = std::numeric_limits<std::size_t>::max());
  Persistence_landscape _multiply_landscape_by_real_number_not_overwrite(double x) const;
  void _multiply_landscape_by_real_number_overwrite(double x);

  // warning, this function can be only called after filling in the intervals vector.
  void _set_up_numbers_of_functions_for_vectorization_and_projections_to_reals()
  {
    this->number_of_functions_for_vectorization_ = this->land_.size();
    this->number_of_functions_for_projections_to_reals_ = this->land_.size();
  }
};

inline bool Persistence_landscape::operator==(const Persistence_landscape& rhs) const
{
  if (this->land_.size() != rhs.land_.size()) {
#ifdef DEBUG_TRACES
    std::clog << "1\n";
#endif
    return false;
  }
  for (std::size_t level = 0; level != this->land_.size(); ++level) {
    if (this->land_[level].size() != rhs.land_[level].size()) {
#ifdef DEBUG_TRACES
      std::clog << "this->land[level].size() : " << this->land_[level].size() << "\n";
      std::clog << "rhs.land[level].size() : " << rhs.land_[level].size() << "\n";
      std::clog << "2\n";
#endif
      return false;
    }
    for (std::size_t i = 0; i != this->land_[level].size(); ++i) {
      if (!(almost_equal(this->land_[level][i].first, rhs.land_[level][i].first) &&
            almost_equal(this->land_[level][i].second, rhs.land_[level][i].second))) {
#ifdef DEBUG_TRACES
        std::clog << "this->land[level][i] : " << this->land_[level][i].first << " " << this->land_[level][i].second
                  << "\n";
        std::clog << "rhs.land[level][i] : " << rhs.land_[level][i].first << " " << rhs.land_[level][i].second << "\n";
        std::clog << "3\n";
#endif
        return false;
      }
    }
  }
  return true;
}

inline void Persistence_landscape::_construct_persistence_landscape_from_barcode(
    const std::vector<std::pair<double, double> >& p,
    std::size_t number_of_levels)
{
#ifdef DEBUG_TRACES
  std::clog << "Persistence_landscape::Persistence_landscape( const std::vector< std::pair< double , double > >& p )"
            << std::endl;
#endif

  // this is a general algorithm to construct persistence landscapes.
  std::vector<std::pair<double, double> > bars;
  bars.insert(bars.begin(), p.begin(), p.end());
  std::sort(bars.begin(), bars.end(), compare_points_sorting);

#ifdef DEBUG_TRACES
  std::clog << "Bars : \n";
  for (std::size_t i = 0; i != bars.size(); ++i) {
    std::clog << bars[i].first << " " << bars[i].second << "\n";
  }
#endif

  std::vector<std::pair<double, double> > characteristicPoints(p.size());
  for (std::size_t i = 0; i != bars.size(); ++i) {
    characteristicPoints[i] =
        std::make_pair((bars[i].first + bars[i].second) / 2.0, (bars[i].second - bars[i].first) / 2.0);
  }
  std::vector<std::vector<std::pair<double, double> > > Persistence_landscape;
  std::size_t number_of_levels_in_the_landscape = 0;
  while (!characteristicPoints.empty()) {
#ifdef DEBUG_TRACES
    for (std::size_t i = 0; i != characteristicPoints.size(); ++i) {
      std::clog << "(" << characteristicPoints[i].first << " " << characteristicPoints[i].second << ")\n";
    }
#endif

    std::vector<std::pair<double, double> > lambda_n;
    lambda_n.push_back(std::make_pair(-std::numeric_limits<int>::max(), 0));
    lambda_n.push_back(std::make_pair(minus_length(characteristicPoints[0]), 0));
    lambda_n.push_back(characteristicPoints[0]);

#ifdef DEBUG_TRACES
    std::clog << "1 Adding to lambda_n : (" << -std::numeric_limits<int>::max() << " " << 0 << ") , ("
              << minus_length(characteristicPoints[0]) << " " << 0 << ") , (" << characteristicPoints[0].first << " "
              << characteristicPoints[0].second << ") \n";
#endif

    std::size_t i = 1;
    std::vector<std::pair<double, double> > newCharacteristicPoints;
    while (i < characteristicPoints.size()) {
      std::size_t p = 1;
      if ((minus_length(characteristicPoints[i]) >= minus_length(lambda_n[lambda_n.size() - 1])) &&
          (birth_plus_deaths(characteristicPoints[i]) > birth_plus_deaths(lambda_n[lambda_n.size() - 1]))) {
        if (minus_length(characteristicPoints[i]) < birth_plus_deaths(lambda_n[lambda_n.size() - 1])) {
          std::pair<double, double> point = std::make_pair(
              (minus_length(characteristicPoints[i]) + birth_plus_deaths(lambda_n[lambda_n.size() - 1])) / 2,
              (birth_plus_deaths(lambda_n[lambda_n.size() - 1]) - minus_length(characteristicPoints[i])) / 2);
          lambda_n.push_back(point);
#ifdef DEBUG_TRACES
          std::clog << "2 Adding to lambda_n : (" << point.first << " " << point.second << ")\n";
#endif

#ifdef DEBUG_TRACES
          std::clog << "characteristicPoints[i+p] : " << characteristicPoints[i + p].first << " "
                    << characteristicPoints[i + p].second << "\n";
          std::clog << "point : " << point.first << " " << point.second << "\n";
#endif

          while ((i + p < characteristicPoints.size()) &&
                 (almost_equal(minus_length(point), minus_length(characteristicPoints[i + p]))) &&
                 (birth_plus_deaths(point) <= birth_plus_deaths(characteristicPoints[i + p]))) {
            newCharacteristicPoints.push_back(characteristicPoints[i + p]);
#ifdef DEBUG_TRACES
            std::clog << "3.5 Adding to newCharacteristicPoints : (" << characteristicPoints[i + p].first << " "
                      << characteristicPoints[i + p].second << ")\n";
#endif
            ++p;
          }

          newCharacteristicPoints.push_back(point);
#ifdef DEBUG_TRACES
          std::clog << "4 Adding to newCharacteristicPoints : (" << point.first << " " << point.second << ")\n";
#endif

          while ((i + p < characteristicPoints.size()) &&
                 (minus_length(point) <= minus_length(characteristicPoints[i + p])) &&
                 (birth_plus_deaths(point) >= birth_plus_deaths(characteristicPoints[i + p]))) {
            newCharacteristicPoints.push_back(characteristicPoints[i + p]);
#ifdef DEBUG_TRACES
            std::clog << "characteristicPoints[i+p] : " << characteristicPoints[i + p].first << " "
                      << characteristicPoints[i + p].second << "\n";
            std::clog << "point : " << point.first << " " << point.second << "\n";
            std::clog << "characteristicPoints[i+p] birth and death : " << minus_length(characteristicPoints[i + p])
                      << " , " << birth_plus_deaths(characteristicPoints[i + p]) << "\n";
            std::clog << "point birth and death : " << minus_length(point) << " , " << birth_plus_deaths(point) << "\n";

            std::clog << "3 Adding to newCharacteristicPoints : (" << characteristicPoints[i + p].first << " "
                      << characteristicPoints[i + p].second << ")\n";
#endif
            ++p;
          }

        } else {
          lambda_n.push_back(std::make_pair(birth_plus_deaths(lambda_n[lambda_n.size() - 1]), 0));
          lambda_n.push_back(std::make_pair(minus_length(characteristicPoints[i]), 0));
#ifdef DEBUG_TRACES
          std::clog << "5 Adding to lambda_n : (" << birth_plus_deaths(lambda_n[lambda_n.size() - 1]) << " " << 0
                    << ")\n";
          std::clog << "5 Adding to lambda_n : (" << minus_length(characteristicPoints[i]) << " " << 0 << ")\n";
#endif
        }
        lambda_n.push_back(characteristicPoints[i]);
#ifdef DEBUG_TRACES
        std::clog << "6 Adding to lambda_n : (" << characteristicPoints[i].first << " "
                  << characteristicPoints[i].second << ")\n";
#endif
      } else {
        newCharacteristicPoints.push_back(characteristicPoints[i]);
#ifdef DEBUG_TRACES
        std::clog << "7 Adding to newCharacteristicPoints : (" << characteristicPoints[i].first << " "
                  << characteristicPoints[i].second << ")\n";
#endif
      }
      i = i + p;
    }
    lambda_n.push_back(std::make_pair(birth_plus_deaths(lambda_n[lambda_n.size() - 1]), 0));
    lambda_n.push_back(std::make_pair(std::numeric_limits<int>::max(), 0));

    characteristicPoints = newCharacteristicPoints;

    lambda_n.erase(std::unique(lambda_n.begin(), lambda_n.end()), lambda_n.end());
    this->land_.push_back(lambda_n);

    ++number_of_levels_in_the_landscape;
    if (number_of_levels == number_of_levels_in_the_landscape) {
      break;
    }
  }
}

// this function find maximum of lambda_n
inline double Persistence_landscape::find_max(unsigned int lambda) const
{
  if (this->land_.size() < lambda) return 0;
  double maximum = -std::numeric_limits<int>::max();
  for (std::size_t i = 0; i != this->land_[lambda].size(); ++i) {
    if (this->land_[lambda][i].second > maximum) maximum = this->land_[lambda][i].second;
  }
  return maximum;
}

inline double Persistence_landscape::compute_integral_of_landscape() const
{
  double result = 0;
  for (std::size_t i = 0; i != this->land_.size(); ++i) {
    for (std::size_t nr = 2; nr != this->land_[i].size() - 1; ++nr) {
      // it suffices to compute every planar integral and then sum them up for each lambda_n
      result += 0.5 * (this->land_[i][nr].first - this->land_[i][nr - 1].first) *
                (this->land_[i][nr].second + this->land_[i][nr - 1].second);
    }
  }
  return result;
}

inline double Persistence_landscape::compute_integral_of_a_level_of_a_landscape(std::size_t level) const
{
  double result = 0;
  if (level >= this->land_.size()) {
    // this landscape function is constantly equal 0, so is the integral.
    return result;
  }
  // also negative landscapes are assumed to be zero.
  if (level < 0) return 0;

  for (std::size_t nr = 2; nr != this->land_[level].size() - 1; ++nr) {
    // it suffices to compute every planar integral and then sum them up for each lambda_n
    result += 0.5 * (this->land_[level][nr].first - this->land_[level][nr - 1].first) *
              (this->land_[level][nr].second + this->land_[level][nr - 1].second);
  }

  return result;
}

inline double Persistence_landscape::compute_integral_of_landscape(double p) const
{
  double result = 0;
  for (std::size_t i = 0; i != this->land_.size(); ++i) {
    for (std::size_t nr = 2; nr != this->land_[i].size() - 1; ++nr) {
#ifdef DEBUG_TRACES
      std::clog << "nr : " << nr << "\n";
#endif
      // In this interval, the landscape has a form f(x) = ax+b. We want to compute integral of (ax+b)^p = 1/a *
      // (ax+b)^{p+1}/(p+1)
      auto [a, b] = compute_parameters_of_a_line(this->land_[i][nr], this->land_[i][nr - 1]);

#ifdef DEBUG_TRACES
      std::clog << "(" << this->land_[i][nr].first << "," << this->land_[i][nr].second << ") , "
                << this->land_[i][nr - 1].first << "," << this->land_[i][nr].second << ")" << std::endl;
#endif
      if (this->land_[i][nr].first == this->land_[i][nr - 1].first) continue;
      if (a != 0) {
        result += 1 / (a * (p + 1)) *
                  (std::pow((a * this->land_[i][nr].first + b), p + 1) -
                   std::pow((a * this->land_[i][nr - 1].first + b), p + 1));
      } else {
        result += (this->land_[i][nr].first - this->land_[i][nr - 1].first) * (std::pow(this->land_[i][nr].second, p));
      }
#ifdef DEBUG_TRACES
      std::clog << "a : " << a << " , b : " << b << std::endl;
      std::clog << "result : " << result << std::endl;
#endif
    }
  }
  return result;
}

// this is O(log(n)) algorithm, where n is number of points in this->land.
inline double Persistence_landscape::compute_value_at_a_given_point(unsigned int level, double x) const
{
  // in such a case lambda_level = 0.
  if (level >= this->land_.size()) return 0;

  // we know that the points in this->land[level] are ordered according to x coordinate. Therefore, we can find the
  // point by using bisection:
  unsigned int coordBegin = 1;
  unsigned int coordEnd = this->land_[level].size() - 2;

#ifdef DEBUG_TRACES
  std::clog << "Here \n";
  std::clog << "x : " << x << "\n";
  std::clog << "this->land[level][coordBegin].first : " << this->land_[level][coordBegin].first << "\n";
  std::clog << "this->land[level][coordEnd].first : " << this->land_[level][coordEnd].first << "\n";
#endif

  // in this case x is outside the support of the landscape, therefore the value of the landscape is 0.
  if (x <= this->land_[level][coordBegin].first) return 0;
  if (x >= this->land_[level][coordEnd].first) return 0;

#ifdef DEBUG_TRACES
  std::clog << "Entering to the while loop \n";
#endif

  while (coordBegin + 1 != coordEnd) {
#ifdef DEBUG_TRACES
    std::clog << "coordBegin : " << coordBegin << "\n";
    std::clog << "coordEnd : " << coordEnd << "\n";
    std::clog << "this->land[level][coordBegin].first : " << this->land_[level][coordBegin].first << "\n";
    std::clog << "this->land[level][coordEnd].first : " << this->land_[level][coordEnd].first << "\n";
#endif

    unsigned int newCord = (unsigned int)floor((coordEnd + coordBegin) / 2.0);

#ifdef DEBUG_TRACES
    std::clog << "newCord : " << newCord << "\n";
    std::clog << "this->land[level][newCord].first : " << this->land_[level][newCord].first << "\n";
#endif

    if (this->land_[level][newCord].first <= x) {
      coordBegin = newCord;
      if (this->land_[level][newCord].first == x) return this->land_[level][newCord].second;
    } else {
      coordEnd = newCord;
    }
  }

#ifdef DEBUG_TRACES
  std::clog << "x : " << x << " is between : " << this->land_[level][coordBegin].first << " a  "
            << this->land_[level][coordEnd].first << "\n";
  std::clog << "the y coords are : " << this->land_[level][coordBegin].second << " a  "
            << this->land_[level][coordEnd].second << "\n";
  std::clog << "coordBegin : " << coordBegin << "\n";
  std::clog << "coordEnd : " << coordEnd << "\n";
#endif
  return function_value(this->land_[level][coordBegin], this->land_[level][coordEnd], x);
}

inline void Persistence_landscape::_multiply_landscape_by_real_number_overwrite(double x)
{
  for (std::size_t dim = 0; dim != this->land_.size(); ++dim) {
    for (std::size_t i = 0; i != this->land_[dim].size(); ++i) {
      this->land_[dim][i].second *= x;
    }
  }
}

inline Persistence_landscape Persistence_landscape::abs()
{
  Persistence_landscape result;
  for (std::size_t level = 0; level != this->land_.size(); ++level) {
#ifdef DEBUG_TRACES
    std::clog << "level: " << level << std::endl;
#endif
    std::vector<std::pair<double, double> > lambda_n;
    lambda_n.push_back(std::make_pair(-std::numeric_limits<int>::max(), 0));
    for (std::size_t i = 1; i != this->land_[level].size(); ++i) {
#ifdef DEBUG_TRACES
      std::clog << "this->land[" << level << "][" << i << "] : " << this->land_[level][i].first << " "
                << this->land_[level][i].second << std::endl;
#endif
      // if a line segment between this->land[level][i-1] and this->land[level][i] crosses the x-axis, then we have to
      // add one landscape point t o result
      if ((this->land_[level][i - 1].second) * (this->land_[level][i].second) < 0) {
        double zero =
            find_zero_of_a_line_segment_between_those_two_points(this->land_[level][i - 1], this->land_[level][i]);

        lambda_n.push_back(std::make_pair(zero, 0));
        lambda_n.push_back(std::make_pair(this->land_[level][i].first, std::fabs(this->land_[level][i].second)));
#ifdef DEBUG_TRACES
        std::clog << "Adding pair : (" << zero << ",0)" << std::endl;
        std::clog << "In the same step adding pair : (" << this->land_[level][i].first << ","
                  << std::fabs(this->land_[level][i].second) << ") " << std::endl;
#endif
      } else {
        lambda_n.push_back(std::make_pair(this->land_[level][i].first, std::fabs(this->land_[level][i].second)));
#ifdef DEBUG_TRACES
        std::clog << "Adding pair : (" << this->land_[level][i].first << "," << std::fabs(this->land_[level][i].second)
                  << ") " << std::endl;
#endif
      }
    }
    result.land_.push_back(lambda_n);
  }
  return result;
}

inline Persistence_landscape* Persistence_landscape::new_abs()
{
  Persistence_landscape* result = new Persistence_landscape(*this);
  for (std::size_t level = 0; level != this->land_.size(); ++level) {
#ifdef DEBUG_TRACES
    std::clog << "level: " << level << std::endl;
#endif
    std::vector<std::pair<double, double> > lambda_n;
    lambda_n.push_back(std::make_pair(-std::numeric_limits<int>::max(), 0));
    for (std::size_t i = 1; i != this->land_[level].size(); ++i) {
#ifdef DEBUG_TRACES
      std::clog << "this->land[" << level << "][" << i << "] : " << this->land_[level][i].first << " "
                << this->land_[level][i].second << std::endl;
#endif
      // if a line segment between this->land[level][i-1] and this->land[level][i] crosses the x-axis, then we have to
      // add one landscape point t o result
      if ((this->land_[level][i - 1].second) * (this->land_[level][i].second) < 0) {
        double zero =
            find_zero_of_a_line_segment_between_those_two_points(this->land_[level][i - 1], this->land_[level][i]);

        lambda_n.push_back(std::make_pair(zero, 0));
        lambda_n.push_back(std::make_pair(this->land_[level][i].first, std::fabs(this->land_[level][i].second)));
#ifdef DEBUG_TRACES
        std::clog << "Adding pair : (" << zero << ",0)" << std::endl;
        std::clog << "In the same step adding pair : (" << this->land_[level][i].first << ","
                  << std::fabs(this->land_[level][i].second) << ") " << std::endl;
#endif
      } else {
        lambda_n.push_back(std::make_pair(this->land_[level][i].first, std::fabs(this->land_[level][i].second)));
#ifdef DEBUG_TRACES
        std::clog << "Adding pair : (" << this->land_[level][i].first << "," << std::fabs(this->land_[level][i].second)
                  << ") " << std::endl;
#endif
      }
    }
    result->land_.push_back(lambda_n);
  }
  return result;
}

inline Persistence_landscape Persistence_landscape::_multiply_landscape_by_real_number_not_overwrite(double x) const
{
  std::vector<std::vector<std::pair<double, double> > > result(this->land_.size());
  for (std::size_t dim = 0; dim != this->land_.size(); ++dim) {
    std::vector<std::pair<double, double> > lambda_dim(this->land_[dim].size());
    for (std::size_t i = 0; i != this->land_[dim].size(); ++i) {
      lambda_dim[i] = std::make_pair(this->land_[dim][i].first, x * this->land_[dim][i].second);
    }
    result[dim] = lambda_dim;
  }
  Persistence_landscape res;
  // CHANGE
  // res.land = result;
  res.land_.swap(result);
  return res;
}

inline void Persistence_landscape::print_to_file(const char* filename) const
{
  std::ofstream write;
  write.open(filename);
  for (std::size_t dim = 0; dim != this->land_.size(); ++dim) {
    write << "#lambda_" << dim << std::endl;
    for (std::size_t i = 1; i != this->land_[dim].size() - 1; ++i) {
      write << this->land_[dim][i].first << "  " << this->land_[dim][i].second << std::endl;
    }
  }
  write.close();
}

inline void Persistence_landscape::load_landscape_from_file(const char* filename)
{
  // removing the current content of the persistence landscape.
  this->land_.clear();

  // this constructor reads persistence landscape form a file. This file have to be created by this software before head
  std::ifstream in;
  in.open(filename);
  if (!in.good()) {
#ifdef DEBUG_TRACES
    std::cerr << "The file : " << filename << " do not exist. The program will now terminate \n";
#endif
    throw std::invalid_argument("The persistence landscape file do not exist.");
  }

  std::string line;
  std::vector<std::pair<double, double> > landscapeAtThisLevel;

  bool isThisAFirsLine = true;
  while (in.good()) {
    getline(in, line);
    if (!(line.length() == 0 || line[0] == '#')) {
      std::stringstream lineSS;
      lineSS << line;
      double begin, end;
      lineSS >> begin;
      lineSS >> end;
      landscapeAtThisLevel.push_back(std::make_pair(begin, end));
#ifdef DEBUG_TRACES
      std::clog << "Reading a point : " << begin << " , " << end << std::endl;
#endif
    } else {
#ifdef DEBUG_TRACES
      std::clog << "IGNORE LINE\n";
#endif
      if (!isThisAFirsLine) {
        landscapeAtThisLevel.push_back(std::make_pair(std::numeric_limits<int>::max(), 0));
        this->land_.push_back(landscapeAtThisLevel);
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
    this->land_.push_back(landscapeAtThisLevel);
  }

  in.close();
}

inline void Persistence_landscape::plot(const char* filename,
                                        double xRangeBegin,
                                        double xRangeEnd,
                                        double yRangeBegin,
                                        double yRangeEnd,
                                        int from,
                                        int to)
{
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
    to = this->land_.size();
  }

  out << "plot ";
  for (std::size_t lambda = std::min((std::size_t)from, this->land_.size());
       lambda != std::min((std::size_t)to, this->land_.size());
       ++lambda) {
    // out << "     '-' using 1:2 title 'l" << lambda << "' with lp";
    out << "     '-' using 1:2 notitle with lp";
    if (lambda + 1 != std::min((std::size_t)to, this->land_.size())) {
      out << ", \\";
    }
    out << std::endl;
  }

  for (std::size_t lambda = std::min((std::size_t)from, this->land_.size());
       lambda != std::min((std::size_t)to, this->land_.size());
       ++lambda) {
    for (std::size_t i = 1; i != this->land_[lambda].size() - 1; ++i) {
      out << this->land_[lambda][i].first << " " << this->land_[lambda][i].second << std::endl;
    }
    out << "EOF" << std::endl;
  }
#ifdef DEBUG_TRACES
  std::clog << "To visualize, install gnuplot and type the command: gnuplot -persist -e \"load \'"
            << gnuplot_script.str().c_str() << "\'\"" << std::endl;
#endif
}

}  // namespace Persistence_representations
}  // namespace Gudhi

#endif  // PERSISTENCE_LANDSCAPE_H_
