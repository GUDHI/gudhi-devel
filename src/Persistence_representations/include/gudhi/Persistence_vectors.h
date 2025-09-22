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

#ifndef PERSISTENCE_VECTORS_H_
#define PERSISTENCE_VECTORS_H_

// standard include
#ifdef DEBUG_TRACES
#include <iostream>   // std::cerr, std::clog
#endif
#include <ostream>    // std::ostream
#include <fstream>    // std::ofstream, std::ifstream
#include <sstream>    // std::stringstream
#include <stdexcept>  // std::invalid_argument
#include <cstddef>    // std::size_t
#include <cmath>      // std::min, std::max, std::pow, std::fabs
#include <algorithm>  // std::make_heap, std::pop_heap, std::push_heap, std::sort
#include <limits>     // std::numeric_limits
#include <utility>    // std::pair
#include <vector>

// gudhi include
#include <gudhi/read_persistence_from_file.h>
#include <gudhi/common_persistence_representations.h>
#include <gudhi/Debug_utils.h>

namespace Gudhi {
namespace Persistence_representations {

// TODO: the template F is only used in the constructor to build `sorted_vector_of_distances_`. Once done it is
// used nowhere anymore or did I miss something?
// Would it not be better to pass the method as argument to the constructor?
// But I guess, that won't be retro-compatible anymore...
/**
 * \class Vector_distances_in_diagram Persistence_vectors.h gudhi/Persistence_vectors.h
 * \brief A class implementing persistence vectors.
 *
 * \ingroup Persistence_representations
 *
 * \details
 * This is an implementation of idea presented in the paper <i>Stable Topological Signatures for Points on 3D
 * Shapes</i> \cite Carriere_Oudot_Ovsjanikov_top_signatures_3d .<br>
 * The parameter of the class is the class that computes distance used to construct the vectors. The typical function
 * is either Euclidean of maximum (Manhattan) distance.
 *
 * This class implements the following concepts: Vectorized_topological_data, Topological_data_with_distances,
 * Real_valued_topological_data, Topological_data_with_averages, Topological_data_with_scalar_product
 **/
template <typename F>
class Vector_distances_in_diagram
{
 public:
  /**
   * The default constructor.
   **/
  Vector_distances_in_diagram() {}

  /**
   * The constructor that takes as an input a multiset of persistence intervals (given as vector of birth-death
   * pairs). The second parameter is the desired length of the output vectors.
   **/
  Vector_distances_in_diagram(const std::vector<std::pair<double, double> >& intervals, std::size_t where_to_cut)
      : intervals_(intervals), where_to_cut_(where_to_cut)
  {
    this->_compute_sorted_vector_of_distances_via_vector_sorting();
    this->_set_up_numbers_of_functions_for_vectorization_and_projections_to_reals();
  }

  /**
   * The constructor taking as an input a file with birth-death pairs. The second parameter is the desired length of
   * the output vectors.
   **/
  Vector_distances_in_diagram(const char* filename,
                              std::size_t where_to_cut,
                              unsigned int dimension = std::numeric_limits<unsigned int>::max())
      : intervals_(read_persistence_intervals_in_one_dimension_from_file(
            filename,
            dimension == std::numeric_limits<unsigned int>::max() ? -1 : static_cast<int>(dimension))),
        where_to_cut_(where_to_cut)
  {
    this->_compute_sorted_vector_of_distances_via_heap();
    _set_up_numbers_of_functions_for_vectorization_and_projections_to_reals();
  }

  /**
   * Writing to a stream.
   **/
  friend std::ostream& operator<<(std::ostream& out, const Vector_distances_in_diagram& d)
  {
    for (std::size_t i = 0; i != std::min(d.sorted_vector_of_distances_.size(), d.where_to_cut_); ++i) {
      out << d.sorted_vector_of_distances_[i] << " ";
    }
    return out;
  }

  /**
   * This procedure gives the value of a vector on a given position.
   **/
  double vector_in_position(std::size_t position) const
  {
    if (position >= this->sorted_vector_of_distances_.size())
      throw("Wrong position in accessing Vector_distances_in_diagram::sorted_vector_of_distances_\n");
    return this->sorted_vector_of_distances_[position];
  }

  /**
   * Return a size of a vector.
   **/
  std::size_t size() const { return this->sorted_vector_of_distances_.size(); }

  /**
   * Write a vector to a file.
   **/
  void write_to_file(const char* filename) const;

  /**
   * Write a vector to a file.
   **/
  void print_to_file(const char* filename) const { this->write_to_file(filename); }

  /**
   * Loading a vector to a file.
   **/
  void load_from_file(const char* filename);

  /**
   * Comparison operators:
   **/
  bool operator==(const Vector_distances_in_diagram& second) const
  {
    if (this->sorted_vector_of_distances_.size() != second.sorted_vector_of_distances_.size()) return false;
    for (std::size_t i = 0; i != this->sorted_vector_of_distances_.size(); ++i) {
      if (!almost_equal(this->sorted_vector_of_distances_[i], second.sorted_vector_of_distances_[i])) return false;
    }
    return true;
  }

  bool operator!=(const Vector_distances_in_diagram& second) const { return !(*this == second); }

  // Implementations of functions for various concepts.
  /**
   * Compute projection to real numbers of persistence vector. This function is required by the
   * Real_valued_topological_data concept
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
   * Compute a vectorization of a persistent vectors. It is required in a concept Vectorized_topological_data.
   **/
  std::vector<double> vectorize(int number_of_function) const;

  /**
   * This function return the number of functions that allows vectorization of a persistence vector. It is required
   * in a concept Vectorized_topological_data.
   **/
  std::size_t number_of_vectorize_functions() const { return this->number_of_functions_for_vectorization_; }

  /**
   * Compute a average of two persistent vectors. This function is required by Topological_data_with_averages concept.
   **/
  void compute_average(const std::vector<Vector_distances_in_diagram*>& to_average);

  /**
   * Compute a distance of two persistent vectors. This function is required in Topological_data_with_distances concept.
   * For max norm distance, set power to std::numeric_limits<double>::max().
   **/
  double distance(const Vector_distances_in_diagram& second, double power = 1) const;

  /**
   * Compute a scalar product of two persistent vectors. This function is required in
   * Topological_data_with_scalar_product concept.
   **/
  double compute_scalar_product(const Vector_distances_in_diagram& second) const;

  // end of implementation of functions needed for concepts.

  /**
   * For visualization use output from vectorize and build histograms.
   **/
  const std::vector<double>& output_for_visualization() const { return this->sorted_vector_of_distances_; }

  /**
   * Create a gnuplot script to visualize the data structure.
   **/
  void plot(const char* filename) const
  {
    std::stringstream gnuplot_script;
    gnuplot_script << filename << "_GnuplotScript";
    std::ofstream out;
    out.open(gnuplot_script.str().c_str());
    out << "set style data histogram" << std::endl;
    out << "set style histogram cluster gap 1" << std::endl;
    out << "set style fill solid border -1" << std::endl;
    out << "plot '-' notitle" << std::endl;
    for (std::size_t i = 0; i != this->sorted_vector_of_distances_.size(); ++i) {
      out << this->sorted_vector_of_distances_[i] << std::endl;
    }
    out << std::endl;
    out.close();
#ifdef DEBUG_TRACES
    std::clog << "To visualize, install gnuplot and type the command: gnuplot -persist -e \"load \'"
              << gnuplot_script.str().c_str() << "\'\"" << std::endl;
#endif
  }

  /**
   * The x-range of the persistence vector.
   **/
  std::pair<double, double> get_x_range() const { return std::make_pair(0, this->sorted_vector_of_distances_.size()); }

  /**
   * The y-range of the persistence vector.
   **/
  std::pair<double, double> get_y_range() const
  {
    if (this->sorted_vector_of_distances_.size() == 0) return std::make_pair(0, 0);
    return std::make_pair(this->sorted_vector_of_distances_[0], 0);
  }

  // arithmetic operations:
  template <typename Operation>
  friend Vector_distances_in_diagram operation_on_pair_of_vectors(const Vector_distances_in_diagram& first,
                                                                  const Vector_distances_in_diagram& second,
                                                                  Operation&& operation)
  {
    Vector_distances_in_diagram result;
    result.sorted_vector_of_distances_.reserve(
        std::max(first.sorted_vector_of_distances_.size(), second.sorted_vector_of_distances_.size()));
    for (std::size_t i = 0;
         i != std::min(first.sorted_vector_of_distances_.size(), second.sorted_vector_of_distances_.size());
         ++i) {
      result.sorted_vector_of_distances_.push_back(
          operation(first.sorted_vector_of_distances_[i], second.sorted_vector_of_distances_[i]));
    }
    if (first.sorted_vector_of_distances_.size() ==
        std::min(first.sorted_vector_of_distances_.size(), second.sorted_vector_of_distances_.size())) {
      for (std::size_t i =
               std::min(first.sorted_vector_of_distances_.size(), second.sorted_vector_of_distances_.size());
           i != std::max(first.sorted_vector_of_distances_.size(), second.sorted_vector_of_distances_.size());
           ++i) {
        result.sorted_vector_of_distances_.push_back(operation(0, second.sorted_vector_of_distances_[i]));
      }
    } else {
      for (std::size_t i =
               std::min(first.sorted_vector_of_distances_.size(), second.sorted_vector_of_distances_.size());
           i != std::max(first.sorted_vector_of_distances_.size(), second.sorted_vector_of_distances_.size());
           ++i) {
        result.sorted_vector_of_distances_.push_back(operation(first.sorted_vector_of_distances_[i], 0));
      }
    }
    return result;
  }

  /**
   * This function implements an operation of multiplying Vector_distances_in_diagram by a scalar.
   **/
  Vector_distances_in_diagram multiply_by_scalar(double scalar) const
  {
    Vector_distances_in_diagram result;
    result.sorted_vector_of_distances_.reserve(this->sorted_vector_of_distances_.size());
    for (std::size_t i = 0; i != this->sorted_vector_of_distances_.size(); ++i) {
      result.sorted_vector_of_distances_.push_back(scalar * this->sorted_vector_of_distances_[i]);
    }
    return result;
  }

  /**
   * This function computes a sum of two objects of a type Vector_distances_in_diagram.
   **/
  friend Vector_distances_in_diagram operator+(const Vector_distances_in_diagram& first,
                                               const Vector_distances_in_diagram& second)
  {
    return operation_on_pair_of_vectors(first, second, std::plus<double>());
  }

  /**
   * This function computes a difference of two objects of a type Vector_distances_in_diagram.
   **/
  friend Vector_distances_in_diagram operator-(const Vector_distances_in_diagram& first,
                                               const Vector_distances_in_diagram& second)
  {
    return operation_on_pair_of_vectors(first, second, std::minus<double>());
  }

  /**
   * This function computes a product of an object of a type Vector_distances_in_diagram with real number.
   **/
  friend Vector_distances_in_diagram operator*(double scalar, const Vector_distances_in_diagram& A)
  {
    return A.multiply_by_scalar(scalar);
  }

  /**
   * This function computes a product of an object of a type Vector_distances_in_diagram with real number.
   **/
  friend Vector_distances_in_diagram operator*(const Vector_distances_in_diagram& A, double scalar)
  {
    return A.multiply_by_scalar(scalar);
  }

  /**
   * This function computes a product of an object of a type Vector_distances_in_diagram with real number.
   **/
  Vector_distances_in_diagram operator*(double scalar) { return this->multiply_by_scalar(scalar); }

  /**
   * += operator for Vector_distances_in_diagram.
   **/
  Vector_distances_in_diagram operator+=(const Vector_distances_in_diagram& rhs)
  {
    *this = *this + rhs;
    return *this;
  }

  /**
   * -= operator for Vector_distances_in_diagram.
   **/
  Vector_distances_in_diagram operator-=(const Vector_distances_in_diagram& rhs)
  {
    *this = *this - rhs;
    return *this;
  }

  /**
   * *= operator for Vector_distances_in_diagram.
   **/
  Vector_distances_in_diagram operator*=(double x)
  {
    *this = *this * x;
    return *this;
  }

  /**
   * /= operator for Vector_distances_in_diagram.
   **/
  Vector_distances_in_diagram operator/=(double x)
  {
    if (x == 0) throw("In operator /=, division by 0. Program terminated.");
    *this = *this * (1 / x);
    return *this;
  }

 private:
  std::vector<std::pair<double, double> > intervals_;
  std::vector<double> sorted_vector_of_distances_;
  std::size_t number_of_functions_for_vectorization_;
  std::size_t number_of_functions_for_projections_to_reals_;
  std::size_t where_to_cut_;

  void _compute_sorted_vector_of_distances_via_heap();
  void _compute_sorted_vector_of_distances_via_vector_sorting();

  Vector_distances_in_diagram(const std::vector<double>& sorted_vector_of_distances)
      : sorted_vector_of_distances_(sorted_vector_of_distances)
  {
    this->_set_up_numbers_of_functions_for_vectorization_and_projections_to_reals();
  }

  // warning, this function can be only called after filling in the intervals vector.
  void _set_up_numbers_of_functions_for_vectorization_and_projections_to_reals()
  {
    this->number_of_functions_for_vectorization_ = this->sorted_vector_of_distances_.size();
    this->number_of_functions_for_projections_to_reals_ = this->sorted_vector_of_distances_.size();
  }
};

template <typename F>
void Vector_distances_in_diagram<F>::_compute_sorted_vector_of_distances_via_heap()
{
#ifdef DEBUG_TRACES
  std::clog << "Here are the intervals : \n";
  for (std::size_t i = 0; i != this->intervals_.size(); ++i) {
    std::clog << this->intervals_[i].first << " , " << this->intervals_[i].second << std::endl;
  }
#endif
  // TODO: should this also modify where_to_cut_ ?
  auto where_to_cut =
      std::min(where_to_cut_,
               (std::size_t)(0.5 * this->intervals_.size() * (this->intervals_.size() - 1) + this->intervals_.size()));

  std::vector<double> heap(where_to_cut, std::numeric_limits<int>::max());
  std::make_heap(heap.begin(), heap.end());
  F f;

  // for every pair of points in the diagram, compute the minimum of their distance, and distance of those points from
  // diagonal
  for (std::size_t i = 0; i < this->intervals_.size(); ++i) {
    for (std::size_t j = i + 1; j < this->intervals_.size(); ++j) {
      double value =
          std::min(f(this->intervals_[i], this->intervals_[j]),
                   std::min(f(this->intervals_[i],
                              std::make_pair(0.5 * (this->intervals_[i].first + this->intervals_[i].second),
                                             0.5 * (this->intervals_[i].first + this->intervals_[i].second))),
                            f(this->intervals_[j],
                              std::make_pair(0.5 * (this->intervals_[j].first + this->intervals_[j].second),
                                             0.5 * (this->intervals_[j].first + this->intervals_[j].second)))));

#ifdef DEBUG_TRACES
      std::clog << "Value : " << value << std::endl;
      std::clog << "heap.front() : " << heap.front() << std::endl;
#endif

      if (-value < heap.front()) {
#ifdef DEBUG_TRACES
        std::clog << "Replacing : " << heap.front() << " with : " << -value << std::endl;
#endif
        // remove the first element from the heap
        std::pop_heap(heap.begin(), heap.end());
        // and put value there instead:
        heap[where_to_cut - 1] = -value;
        std::push_heap(heap.begin(), heap.end());
      }
    }
  }

  // now add distances of all points from diagonal
  for (std::size_t i = 0; i < this->intervals_.size(); ++i) {
    double value = f(this->intervals_[i],
                     std::make_pair(0.5 * (this->intervals_[i].first + this->intervals_[i].second),
                                    0.5 * (this->intervals_[i].first + this->intervals_[i].second)));
    if (-value < heap.front()) {
      // remove the first element from the heap
      std::pop_heap(heap.begin(), heap.end());
      // and put value there instead:
      heap[where_to_cut - 1] = -value;
      std::push_heap(heap.begin(), heap.end());
    }
  }

  std::sort_heap(heap.begin(), heap.end());
  for (std::size_t i = 0; i != heap.size(); ++i) {
    if (heap[i] == std::numeric_limits<int>::max()) {
      heap[i] = 0;
    } else {
      heap[i] *= -1;
    }
  }

#ifdef DEBUG_TRACES
  std::clog << "This is the heap after all the operations :\n";
  for (std::size_t i = 0; i != heap.size(); ++i) {
    std::clog << heap[i] << " ";
  }
  std::clog << std::endl;
#endif

  this->sorted_vector_of_distances_.swap(heap);
}

template <typename F>
void Vector_distances_in_diagram<F>::_compute_sorted_vector_of_distances_via_vector_sorting()
{
  std::vector<double> distances;
  distances.reserve(
      (std::size_t)(0.5 * this->intervals_.size() * (this->intervals_.size() - 1) + this->intervals_.size()));
  F f;

  // for every pair of points in the diagram, compute the minimum of their distance, and distance of those points from
  // diagonal
  for (std::size_t i = 0; i < this->intervals_.size(); ++i) {
    // add distance of i-th point in the diagram from the diagonal to the distances vector
    distances.push_back(f(this->intervals_[i],
                          std::make_pair(0.5 * (this->intervals_[i].first + this->intervals_[i].second),
                                         0.5 * (this->intervals_[i].first + this->intervals_[i].second))));
    for (std::size_t j = i + 1; j < this->intervals_.size(); ++j) {
      double value =
          std::min(f(this->intervals_[i], this->intervals_[j]),
                   std::min(f(this->intervals_[i],
                              std::make_pair(0.5 * (this->intervals_[i].first + this->intervals_[i].second),
                                             0.5 * (this->intervals_[i].first + this->intervals_[i].second))),
                            f(this->intervals_[j],
                              std::make_pair(0.5 * (this->intervals_[j].first + this->intervals_[j].second),
                                             0.5 * (this->intervals_[j].first + this->intervals_[j].second)))));
      distances.push_back(value);
    }
  }
  std::sort(distances.begin(), distances.end(), std::greater<double>());
  if (distances.size() > where_to_cut_) distances.resize(where_to_cut_);

  this->sorted_vector_of_distances_.swap(distances);
}

// Implementations of functions for various concepts.
template <typename F>
double Vector_distances_in_diagram<F>::project_to_R(int number_of_function) const
{
  // TODO: why not type `number_of_function` directly as `std::size_t`? Are there cases where the user wants this
  // method to throw when `number_of_function` < 0 ?
  if (static_cast<std::size_t>(number_of_function) > this->number_of_functions_for_projections_to_reals_)
    throw std::invalid_argument("Wrong index of a function in a method Vector_distances_in_diagram<F>::project_to_R");
  if (number_of_function < 0)
    throw std::invalid_argument("Wrong index of a function in a method Vector_distances_in_diagram<F>::project_to_R");

  double result = 0;
  for (std::size_t i = 0; i != static_cast<std::size_t>(number_of_function); ++i) {
    result += sorted_vector_of_distances_[i];
  }
  return result;
}

template <typename F>
void Vector_distances_in_diagram<F>::compute_average(const std::vector<Vector_distances_in_diagram*>& to_average)
{
  if (to_average.size() == 0) {
    (*this) = Vector_distances_in_diagram<F>();
    return;
  }

  std::size_t maximal_length_of_vector = 0;
  for (std::size_t i = 0; i != to_average.size(); ++i) {
    if (to_average[i]->sorted_vector_of_distances_.size() > maximal_length_of_vector) {
      maximal_length_of_vector = to_average[i]->sorted_vector_of_distances_.size();
    }
  }

  std::vector<double> av(maximal_length_of_vector, 0);
  for (std::size_t i = 0; i != to_average.size(); ++i) {
    for (std::size_t j = 0; j != to_average[i]->sorted_vector_of_distances_.size(); ++j) {
      av[j] += to_average[i]->sorted_vector_of_distances_[j];
    }
  }

  for (std::size_t i = 0; i != maximal_length_of_vector; ++i) {
    av[i] /= static_cast<double>(to_average.size());
  }
  this->sorted_vector_of_distances_.swap(av);
  this->where_to_cut_ = av.size();
}

template <typename F>
double Vector_distances_in_diagram<F>::distance(const Vector_distances_in_diagram& second_, double power) const
{
#ifdef DEBUG_TRACES
  std::clog << "Entering double Vector_distances_in_diagram<F>::distance( const Abs_Topological_data_with_distances* "
               "second , double power ) procedure \n";
  std::clog << "Power : " << power << std::endl;
  std::clog << "This : " << *this << std::endl;
  std::clog << "second : " << second_ << std::endl;
#endif

  double result = 0;
  for (std::size_t i = 0;
       i != std::min(this->sorted_vector_of_distances_.size(), second_.sorted_vector_of_distances_.size());
       ++i) {
    if (power == 1) {
#ifdef DEBUG_TRACES
      std::clog << "|" << this->sorted_vector_of_distances_[i] << " -  " << second_.sorted_vector_of_distances_[i]
                << " |  : " << std::fabs(this->sorted_vector_of_distances_[i] - second_.sorted_vector_of_distances_[i])
                << std::endl;
#endif
      result += std::fabs(this->sorted_vector_of_distances_[i] - second_.sorted_vector_of_distances_[i]);
    } else {
      if (power < std::numeric_limits<double>::max()) {
        result +=
            std::pow(std::fabs(this->sorted_vector_of_distances_[i] - second_.sorted_vector_of_distances_[i]), power);
      } else {
        // max norm
        if (result < std::fabs(this->sorted_vector_of_distances_[i] - second_.sorted_vector_of_distances_[i]))
          result = std::fabs(this->sorted_vector_of_distances_[i] - second_.sorted_vector_of_distances_[i]);
      }
#ifdef DEBUG_TRACES
      std::clog << "| " << this->sorted_vector_of_distances_[i] << " - " << second_.sorted_vector_of_distances_[i]
                << " : " << std::fabs(this->sorted_vector_of_distances_[i] - second_.sorted_vector_of_distances_[i])
                << std::endl;
#endif
    }
  }
  if (this->sorted_vector_of_distances_.size() != second_.sorted_vector_of_distances_.size()) {
    if (this->sorted_vector_of_distances_.size() > second_.sorted_vector_of_distances_.size()) {
      for (std::size_t i = second_.sorted_vector_of_distances_.size(); i != this->sorted_vector_of_distances_.size();
           ++i) {
        result += std::fabs(this->sorted_vector_of_distances_[i]);
      }
    } else {
      // this->sorted_vector_of_distances_.size() < second_.sorted_vector_of_distances_.size()
      for (std::size_t i = this->sorted_vector_of_distances_.size(); i != second_.sorted_vector_of_distances_.size();
           ++i) {
        result += std::fabs(second_.sorted_vector_of_distances_[i]);
      }
    }
  }

  if (power != 1) {
    result = std::pow(result, (1.0 / power));
  }
  return result;
}

template <typename F>
std::vector<double> Vector_distances_in_diagram<F>::vectorize(int number_of_function) const
{
  // TODO: why not type `number_of_function` directly as `std::size_t`? Are there cases where the user wants this
  // method to throw when `number_of_function` < 0 ?
  if (static_cast<std::size_t>(number_of_function) > this->number_of_functions_for_vectorization_)
    throw std::invalid_argument("Wrong index of a function in a method Vector_distances_in_diagram<F>::vectorize");
  if (number_of_function < 0)
    throw std::invalid_argument("Wrong index of a function in a method Vector_distances_in_diagram<F>::vectorize");

  std::vector<double> result(
      std::min(static_cast<std::size_t>(number_of_function), this->sorted_vector_of_distances_.size()));
  for (std::size_t i = 0;
       i != std::min(static_cast<std::size_t>(number_of_function), this->sorted_vector_of_distances_.size());
       ++i) {
    result[i] = this->sorted_vector_of_distances_[i];
  }
  return result;
}

template <typename F>
void Vector_distances_in_diagram<F>::write_to_file(const char* filename) const
{
  std::ofstream out;
  out.open(filename);

  for (std::size_t i = 0; i != this->sorted_vector_of_distances_.size(); ++i) {
    out << this->sorted_vector_of_distances_[i] << " ";
  }

  out.close();
}

template <typename F>
void Vector_distances_in_diagram<F>::load_from_file(const char* filename)
{
  std::ifstream in;
  in.open(filename);
  // check if the file exist.
  if (!in.good()) {
#ifdef DEBUG_TRACES
    std::cerr << "The file : " << filename << " do not exist. The program will now terminate \n";
#endif
    throw std::invalid_argument("The persistence landscape file do not exist.");
  }

  double number;
  while (in >> number) {
    this->sorted_vector_of_distances_.push_back(number);
  }
  in.close();
}

template <typename F>
double Vector_distances_in_diagram<F>::compute_scalar_product(const Vector_distances_in_diagram& second_vector) const
{
  double result = 0;
  for (std::size_t i = 0;
       i != std::min(this->sorted_vector_of_distances_.size(), second_vector.sorted_vector_of_distances_.size());
       ++i) {
    result += this->sorted_vector_of_distances_[i] * second_vector.sorted_vector_of_distances_[i];
  }
  return result;
}

}  // namespace Persistence_representations
}  // namespace Gudhi

#endif  // PERSISTENCE_VECTORS_H_
