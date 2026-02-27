/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Pawel Dlotko and Mathieu Carriere
 *
 *    Copyright (C) 2019 Inria
 *
 *    Modification(s):
 *      - 2018/04 MC: Add discrete/non-discrete mechanism and non-discrete version
 *      - 2025/06 Hannah Schreiber: Various small bug fixes (missing `inline`s, `DEBUG_TRACES`s etc.)
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PERSISTENCE_HEAT_MAPS_H_
#define PERSISTENCE_HEAT_MAPS_H_

// standard include
#ifdef DEBUG_TRACES
#include <iostream>   // std::clog, std::cerr
#endif
#include <fstream>    // std::ofstream, std::ifstream
#include <sstream>    // std::stringstream
#include <cstddef>    // std::size_t
#include <stdexcept>  // std::invalid_argument
#include <cmath>      // std::sqrt, std::pow, std::atan, std::exp
#include <limits>     // std::numeric_limits
#include <algorithm>  // std::nth_element
#include <functional> // std::function
#include <utility>    // std::pair
#include <vector>
#include <string>

// gudhi include
#include <gudhi/read_persistence_from_file.h>
#include <gudhi/common_persistence_representations.h>
#include <gudhi/Debug_utils.h>

namespace Gudhi {
namespace Persistence_representations {

/**
 * This is a simple procedure to create n by n (or 2*pixel_radius times 2*pixel_radius) cubical approximation of a
 * Gaussian kernel.
 *
 * @ingroup Persistence_representations
 **/
inline std::vector<std::vector<double> > create_Gaussian_filter(std::size_t pixel_radius, double sigma, int approximation=2)
{
  // we are computing the kernel mask to 2 standard deviations away from the center. We discretize it in a grid of a
  // size 2*pixel_radius times 2*pixel_radius.

  double r = 0;
  double sigma_sqr = sigma * sigma;

  // sum is for normalization
  double sum = 0;

  // initialization of a kernel:
  std::vector<std::vector<double> > kernel(2 * pixel_radius + 1, std::vector<double>(2 * pixel_radius + 1, 0));

#ifdef DEBUG_TRACES
  std::clog << "Kernel initialize \n";
  std::clog << "pixel_radius : " << pixel_radius << std::endl;
  std::clog << "kernel.size() : " << kernel.size() << std::endl;
#endif

  for (int x = -pixel_radius; x <= static_cast<int>(pixel_radius); x++) {
    for (int y = -pixel_radius; y <= static_cast<int>(pixel_radius); y++) {
      double real_x = approximation * sigma * x / pixel_radius;
      double real_y = approximation * sigma * y / pixel_radius;
      r = std::sqrt(real_x * real_x + real_y * real_y);
      kernel[x + pixel_radius][y + pixel_radius] = (std::exp(-(r * r) / (2 * sigma_sqr))) / (2 * 3.141592 * sigma_sqr);
      sum += kernel[x + pixel_radius][y + pixel_radius];
    }
  }

  // normalize the kernel
  for (std::size_t i = 0; i != kernel.size(); ++i) {
    for (std::size_t j = 0; j != kernel[i].size(); ++j) {
      kernel[i][j] /= sum;
    }
  }

#ifdef DEBUG_TRACES
  std::clog << "Here is the kernel : \n";
  for (std::size_t i = 0; i != kernel.size(); ++i) {
    for (std::size_t j = 0; j != kernel[i].size(); ++j) {
      std::clog << kernel[i][j] << " ";
    }
    std::clog << std::endl;
  }
#endif
  return kernel;
}

// There are various options to scale the points depending on their location. One can for instance:
// (1) do nothing (scale all of them with the weight 1), as with constant_scaling_function
// (2) Scale them by the distance to the diagonal, as with distance_from_diagonal_scaling
// (3) Scale them with the square of their distance to diagonal, as with squared_distance_from_diagonal_scaling
// (4) Scale them with arc_tan_of_persistence_of_point
// (5) Scale them with weight_by_setting_maximal_interval_to_have_length_one

// TODO: replace all std::pow(x, 2) with x * x ?

/**
 * This is one of a scaling functions used to weight points depending on their persistence and/or location in the
 * diagram.
 * This particular functionality is a function which always assign value 1 to a point in the diagram.
 *
 * \ingroup Persistence_representations
 **/
class constant_scaling_function
{
 public:
  double operator()(const std::pair<double, double>& point_in_diagram) { return 1; }
};

/**
 * This is one of a scaling functions used to weight points depending on their persistence and/or location in the
 * diagram.
 * The scaling given by this function to a point (b,d) is Euclidean distance of (b,d) from diagonal.
 *
 * \ingroup Persistence_representations
 **/
class distance_from_diagonal_scaling
{
 public:
  double operator()(const std::pair<double, double>& point_in_diagram)
  {
    // (point_in_diagram.first+point_in_diagram.second)/2.0
    return std::sqrt(std::pow((point_in_diagram.first - (point_in_diagram.first + point_in_diagram.second) / 2.0), 2) +
                     std::pow((point_in_diagram.second - (point_in_diagram.first + point_in_diagram.second) / 2.0), 2));
  }
};

/**
 * This is one of a scaling functions used to weight points depending on their persistence and/or location in the
 * diagram.
 * The scaling given by this function to a point (b,d) is a square of Euclidean distance of (b,d) from diagonal.
 *
 * \ingroup Persistence_representations
 **/
class squared_distance_from_diagonal_scaling
{
 public:
  double operator()(const std::pair<double, double>& point_in_diagram)
  {
    return std::pow((point_in_diagram.first - (point_in_diagram.first + point_in_diagram.second) / 2.0), 2) +
           std::pow((point_in_diagram.second - (point_in_diagram.first + point_in_diagram.second) / 2.0), 2);
  }
};

/**
 * This is one of a scaling functions used to weight points depending on their persistence and/or location in the
 * diagram.
 * The scaling given by this function to a point \f$ (b,d) \f$ is the arc tangent of the length of the bar
 * corresponding to this point (i.e. \f$ arctangent(b-d) \f$).
 *
 * \ingroup Persistence_representations
 **/
class arc_tan_of_persistence_of_point
{
 public:
  double operator()(const std::pair<double, double>& point_in_diagram)
  {
    return std::atan(point_in_diagram.second - point_in_diagram.first);
  }
};

/**
 * This is one of a scaling functions used to weight points depending on their persistence and/or location in the
 * diagram.
 * This scaling function does not only depend on a point in the diagram, but on the whole diagram.
 * The longest persistence pair gets a scaling 1. Any other pair gets a scaling between 0 and 1, which is proportional
 * to its length compared to the longest.
 *
 * \ingroup Persistence_representations
 **/
class weight_by_setting_maximal_interval_to_have_length_one
{
 public:
  weight_by_setting_maximal_interval_to_have_length_one(double len) : length_of_maximal_interval(len) {}

  double operator()(const std::pair<double, double>& point_in_diagram)
  {
    return (point_in_diagram.second - point_in_diagram.first) / this->length_of_maximal_interval;
  }

 private:
  double length_of_maximal_interval;
};

/**
 * \class Persistence_heat_maps Persistence_heat_maps.h gudhi/Persistence_heat_maps.h
 * \brief A class implementing persistence heat maps.
 *
 * \ingroup Persistence_representations
 *
 * This class implements the following concepts: Vectorized_topological_data, Topological_data_with_distances,
 * Real_valued_topological_data, Topological_data_with_averages, Topological_data_with_scalar_product
 **/
template <typename Scaling_of_kernels = constant_scaling_function>
class Persistence_heat_maps
{
 public:
  /**
   * The default constructor. A scaling function from the diagonal is set up to a constant function. The image is not
   * erased below the diagonal. The Gaussians have diameter 5.
   **/
  Persistence_heat_maps()
      : min_(0),
        max_(0),
        number_of_functions_for_vectorization_(1),
        number_of_functions_for_projections_to_reals_(1),
        f_(),
        erase_below_diagonal_(false)
  {}

  /**
   * @brief Constructor. All parameters except `interval` are optional.
   *
   * @param interval A vector of pairs of doubles (representing persistence intervals).
   * @param filter A Gaussian filter generated by @ref create_Gaussian_filter filter "". Default value: a Gaussian
   * filter of a radius 5.
   * @param erase_below_diagonal A boolean which determines if the area of image below the diagonal should be erased
   * or not. Default value: false (it will not be erased).
   * @param number_of_pixels The number of pixels in each direction. Default value: 1000.
   * @param min The minimal x and y value of points that are to be taken into account. If set to the same value than
   * `max`, its value is re-inferred from the minimal point in the given interval. Assumed to be less or equal to
   * `max`. Default value: `std::numeric_limits<double>::max()`.
   * @param max The maximal x and y value of points that are to be taken into account. If set to the same value than
   * `min`, its value is re-inferred from the maximal point in the given interval. Assumed to be greater or equal to
   * `min`. Default value: `std::numeric_limits<double>::max()`.
   */
  Persistence_heat_maps(const std::vector<std::pair<double, double> >& interval,
                        const std::vector<std::vector<double> >& filter = create_Gaussian_filter(5, 1),
                        bool erase_below_diagonal = false,
                        std::size_t number_of_pixels = 1000,
                        double min = std::numeric_limits<double>::max(),
                        double max = std::numeric_limits<double>::max())
      : number_of_functions_for_vectorization_(1), number_of_functions_for_projections_to_reals_(1)
  {
    this->_construct(interval, filter, erase_below_diagonal, number_of_pixels, min, max);
  }

  /**
   * @brief Constructor from a file. All parameters except `interval` are optional.
   *
   * @param filename The name of a file with persistence intervals. The file should be readable by the function
   * @ref read_persistence_intervals_in_one_dimension_from_file "".
   * @param filter A Gaussian filter generated by @ref create_Gaussian_filter filter "". Default value: a Gaussian
   * filter of a radius 5.
   * @param erase_below_diagonal A boolean which determines if the area of image below the diagonal should be erased
   * or not. Default value: false (it will not be erased).
   * @param number_of_pixels The number of pixels in each direction. Default value: 1000.
   * @param min The minimal x and y value of points that are to be taken into account. If set to the same value than
   * `max`, its value is re-inferred from the minimal point in the given interval. Assumed to be less or equal to
   * `max`. Default value: `std::numeric_limits<double>::max()`.
   * @param max The maximal x and y value of points that are to be taken into account. If set to the same value than
   * `min`, its value is re-inferred from the maximal point in the given interval. Assumed to be greater or equal to
   * `min`. Default value: `std::numeric_limits<double>::max()`.
   * @param dimension If anything other than `std::numeric_limits<unsigned int>::max()`, only the intervals in
   * this given dimension are token into account. Default value: `std::numeric_limits<double>::max()`.
   */
  Persistence_heat_maps(const char* filename,
                        const std::vector<std::vector<double> >& filter = create_Gaussian_filter(5, 1),
                        bool erase_below_diagonal = false,
                        std::size_t number_of_pixels = 1000,
                        double min = std::numeric_limits<double>::max(),
                        double max = std::numeric_limits<double>::max(),
                        unsigned int dimension = std::numeric_limits<unsigned int>::max())
      : number_of_functions_for_vectorization_(1), number_of_functions_for_projections_to_reals_(1)
  {
    int dim = dimension == std::numeric_limits<unsigned int>::max() ? -1 : static_cast<int>(dimension);
    auto intervals_ = read_persistence_intervals_in_one_dimension_from_file(filename, dim);
    this->_construct(intervals_, filter, erase_below_diagonal, number_of_pixels, min, max);
  }

  /**
   * Construction that takes as inputs (1) the diagram, and (2) the grid parameters (min, max and number of samples for
   * x and y axes).
   **/
  Persistence_heat_maps(
      const std::vector<std::pair<double, double> >& interval,
      const std::function<double(const std::pair<double, double>&, const std::pair<double, double>&)>& kernel,
      std::size_t number_of_x_pixels,
      std::size_t number_of_y_pixels,
      double min_x = 0,
      double max_x = 1,
      double min_y = 0,
      double max_y = 1)
      : min_(min_x),
        max_(max_x),
        heat_map_(number_of_y_pixels),
        discrete_(true),
        number_of_functions_for_vectorization_(1),
        number_of_functions_for_projections_to_reals_(1)
  {
    double step_x = (max_x - min_x) / (number_of_x_pixels - 1);
    double step_y = (max_y - min_y) / (number_of_y_pixels - 1);

    int num_pts = interval.size();
    for (std::size_t i = 0; i < number_of_y_pixels; i++) {
      double y = min_y + i * step_y;
      this->heat_map_[i].reserve(number_of_x_pixels);
      for (std::size_t j = 0; j < number_of_x_pixels; j++) {
        double x = min_x + j * step_x;
        std::pair<double, double> grid_point(x, y);
        double pixel_value = 0;
        for (int k = 0; k < num_pts; k++) pixel_value += this->f_(interval[k]) * kernel(interval[k], grid_point);
        this->heat_map_[i].push_back(pixel_value);
      }
    }
  }

  /**
   * Construction that takes only the diagram as input (weight and 2D kernel are template parameters)
   **/
  Persistence_heat_maps(
      const std::vector<std::pair<double, double> >& interval,
      const std::function<double(const std::pair<double, double>&, const std::pair<double, double>&)>& kernel)
      : discrete_(false),
        number_of_functions_for_vectorization_(1),
        number_of_functions_for_projections_to_reals_(1),
        kernel_(kernel),
        interval_(interval)
  {
    int num_pts = this->interval_.size();
    for (int i = 0; i < num_pts; i++) this->weights_.push_back(this->f_(this->interval_[i]));
  }

  /**
   * Compute a mean value of a collection of heat maps and store it in the current object. Note that all the persistence
   * maps send in a vector to this procedure need to have the same parameters.
   * If this is not the case, the program will throw an exception.
   **/
  void compute_mean(const std::vector<Persistence_heat_maps*>& maps);

  /**
   * Compute a median value of a collection of heat maps and store it in the current object. Note that all the
   * persistence maps send in a vector to this procedure need to have the same parameters.
   * If this is not the case, the program will throw an exception.
   **/
  void compute_median(const std::vector<Persistence_heat_maps*>& maps);

  /**
   * Compute a percentage of active (i.e) values above the cutoff of a collection of heat maps.
   **/
  void compute_percentage_of_active(const std::vector<Persistence_heat_maps*>& maps, std::size_t cutoff = 1);

  // put to file subroutine
  /**
   * The function outputs the persistence image to a text file. The format as follow:
   * In the first line, the values min and max of the image are stored
   * In the next lines, we have the persistence images in a form of a bitmap image.
   **/
  void print_to_file(const char* filename) const;

  /**
   * A function that load a heat map from file to the current object (and erase whatever was stored in the current
   * object before).
   **/
  void load_from_file(const char* filename);

  // TODO: replace all use of this method with a GUDHI_CHECK?
  /**
   * The procedure checks if min_, max_ and this->heat_maps_ sizes are the same.
   **/
  bool check_if_the_same(const Persistence_heat_maps& second) const
  {
    if (this->heat_map_.size() != second.heat_map_.size()) {
#ifdef DEBUG_TRACES
      std::clog << "this->heat_map_.size() : " << this->heat_map_.size()
                << " \n second.heat_map_.size() : " << second.heat_map_.size() << std::endl;
#endif
      return false;
    }
    if (this->min_ != second.min_) {
#ifdef DEBUG_TRACES
      std::clog << "this->min_ : " << this->min_ << ", second.min_ : " << second.min_ << std::endl;
#endif
      return false;
    }
    if (this->max_ != second.max_) {
#ifdef DEBUG_TRACES
      std::clog << "this->max_ : " << this->max_ << ", second.max_ : " << second.max_ << std::endl;
#endif
      return false;
    }
    // in the other case we may assume that the persistence images are defined on the same domain.
    return true;
  }

  /**
   * Return minimal range value of persistent image.
   **/
  double get_min() const { return this->min_; }

  /**
   * Return maximal range value of persistent image.
   **/
  double get_max() const { return this->max_; }

  /**
   * Operator == to check if to persistence heat maps are the same.
   **/
  bool operator==(const Persistence_heat_maps& rhs) const
  {
    if (!this->check_if_the_same(rhs)) {
#ifdef DEBUG_TRACES
      std::clog << "The domains are not the same \n";
#endif
      return false;  // in this case, the domains are not the same, so the maps cannot be the same.
    }
    for (std::size_t i = 0; i != this->heat_map_.size(); ++i) {
      for (std::size_t j = 0; j != this->heat_map_[i].size(); ++j) {
        if (!almost_equal(this->heat_map_[i][j], rhs.heat_map_[i][j])) {
#ifdef DEBUG_TRACES
          std::clog << "this->heat_map_[" << i << "][" << j << "] = " << this->heat_map_[i][j] << std::endl;
          std::clog << "rhs.heat_map_[" << i << "][" << j << "] = " << rhs.heat_map_[i][j] << std::endl;
#endif
          return false;
        }
      }
    }
    return true;
  }

  /**
   * Operator != to check if to persistence heat maps are different.
   **/
  bool operator!=(const Persistence_heat_maps& rhs) const { return !((*this) == rhs); }

  /**
   * A function to generate a gnuplot script to visualize the persistent image.
   **/
  void plot(const char* filename) const;

  template <typename Operation_type>
  friend Persistence_heat_maps operation_on_pair_of_heat_maps(const Persistence_heat_maps& first,
                                                              const Persistence_heat_maps& second,
                                                              Operation_type&& operation)
  {
    // first check if the heat maps are compatible
    if (!first.check_if_the_same(second)) {
#ifdef DEBUG_TRACES
      std::cerr << "Sizes of the heat maps are not compatible. The program will now terminate \n";
#endif
      throw std::invalid_argument("Sizes of the heat maps are not compatible. The program will now terminate \n");
    }
    Persistence_heat_maps result;
    result.min_ = first.min_;
    result.max_ = first.max_;
    result.heat_map_.reserve(first.heat_map_.size());
    for (std::size_t i = 0; i != first.heat_map_.size(); ++i) {
      std::vector<double> v;
      v.reserve(first.heat_map_[i].size());
      for (std::size_t j = 0; j != first.heat_map_[i].size(); ++j) {
        v.push_back(operation(first.heat_map_[i][j], second.heat_map_[i][j]));
      }
      result.heat_map_.push_back(v);
    }
    return result;
  }  // operation_on_pair_of_heat_maps

  // TODO: would this method not make more sense as a friend to mirror `operation_on_pair_of_heat_maps`?
  // Or as modifying the heat map itself instead of a copy (and invert the hierarchy of operator* and operator*=) ?
  /**
   * Multiplication of Persistence_heat_maps by scalar (so that all values of the heat map gets multiplied by that
   * scalar).
   **/
  Persistence_heat_maps multiply_by_scalar(double scalar) const
  {
    Persistence_heat_maps result;
    result.min_ = this->min_;
    result.max_ = this->max_;
    result.heat_map_.reserve(this->heat_map_.size());
    for (std::size_t i = 0; i != this->heat_map_.size(); ++i) {
      std::vector<double> v;
      v.reserve(this->heat_map_[i].size());
      for (std::size_t j = 0; j != this->heat_map_[i].size(); ++j) {
        v.push_back(this->heat_map_[i][j] * scalar);
      }
      result.heat_map_.push_back(v);
    }
    return result;
  }

  /**
   * This function computes a sum of two objects of a type Persistence_heat_maps.
   **/
  friend Persistence_heat_maps operator+(const Persistence_heat_maps& first, const Persistence_heat_maps& second)
  {
    return operation_on_pair_of_heat_maps(first, second, std::plus<double>());
  }

  /**
   * This function computes a difference of two objects of a type Persistence_heat_maps.
   **/
  friend Persistence_heat_maps operator-(const Persistence_heat_maps& first, const Persistence_heat_maps& second)
  {
    return operation_on_pair_of_heat_maps(first, second, std::minus<double>());
  }

  /**
   * This function computes a product of an object of a type Persistence_heat_maps with real number.
   **/
  friend Persistence_heat_maps operator*(double scalar, const Persistence_heat_maps& A)
  {
    return A.multiply_by_scalar(scalar);
  }

  /**
   * This function computes a product of an object of a type Persistence_heat_maps with real number.
   **/
  friend Persistence_heat_maps operator*(const Persistence_heat_maps& A, double scalar)
  {
    return A.multiply_by_scalar(scalar);
  }

  /**
   * This function computes a product of an object of a type Persistence_heat_maps with real number.
   **/
  Persistence_heat_maps operator*(double scalar) { return this->multiply_by_scalar(scalar); }

  /**
   * += operator for Persistence_heat_maps.
   **/
  Persistence_heat_maps& operator+=(const Persistence_heat_maps& rhs)
  {
    *this = *this + rhs;
    return *this;
  }

  /**
   * -= operator for Persistence_heat_maps.
   **/
  Persistence_heat_maps& operator-=(const Persistence_heat_maps& rhs)
  {
    *this = *this - rhs;
    return *this;
  }

  /**
   * *= operator for Persistence_heat_maps.
   **/
  Persistence_heat_maps& operator*=(double x)
  {
    *this = *this * x;
    return *this;
  }

  /**
   * /= operator for Persistence_heat_maps.
   **/
  Persistence_heat_maps& operator/=(double x)
  {
    if (x == 0) throw std::invalid_argument("In operator /=, division by 0. Program terminated.");
    *this = *this * (1 / x);
    return *this;
  }

  // Implementations of functions for various concepts.

  /**
   * This function produce a vector of doubles based on a persistence heat map. It is required in a concept
   * Vectorized_topological_data
   */
  std::vector<double> vectorize(int number_of_function) const;

  /**
   * This function return the number of functions that allows vectorization of persistence heat map. It is required
   * in a concept Vectorized_topological_data.
   **/
  std::size_t number_of_vectorize_functions() const { return this->number_of_functions_for_vectorization_; }

  /**
   * This function is required by the Real_valued_topological_data concept. It returns various projections on the
   * persistence heat map to a real line.
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
   * A function to compute distance between persistence heat maps.
   * The parameter of this function is a const reference to an object of a class Persistence_heat_maps.
   * This function is required in Topological_data_with_distances concept.
   * For max norm distance, set power to std::numeric_limits<double>::max()
   **/
  double distance(const Persistence_heat_maps& second_, double power = 1) const;

  /**
   * A function to compute averaged persistence heat map, based on vector of persistence heat maps.
   * This function is required by Topological_data_with_averages concept.
   **/
  void compute_average(const std::vector<Persistence_heat_maps*>& to_average);

  /**
   * A function to compute scalar product of persistence heat maps.
   * The parameter of this function is a const reference to an object of a class Persistence_heat_maps.
   * This function is required in Topological_data_with_scalar_product concept.
   **/
  double compute_scalar_product(const Persistence_heat_maps& second_) const;

  // end of implementation of functions needed for concepts.

  /**
   * The x-range of the persistence heat map.
   **/
  std::pair<double, double> get_x_range() const { return std::make_pair(this->min_, this->max_); }

  /**
   * The y-range of the persistence heat map.
   **/
  std::pair<double, double> get_y_range() const { return this->get_x_range(); }

 protected:
  double min_;
  double max_;
  std::vector<std::vector<double> > heat_map_;

 private:
  std::vector<std::vector<double> > _check_and_initialize_maps(const std::vector<Persistence_heat_maps*>& maps);
  void _construct(const std::vector<std::pair<double, double> >& intervals,
                 const std::vector<std::vector<double> >& filter = create_Gaussian_filter(5, 1),
                 bool erase_below_diagonal = false,
                 std::size_t number_of_pixels = 1000,
                 double min = std::numeric_limits<double>::max(),
                 double max = std::numeric_limits<double>::max());

  // Boolean indicating if we are computing persistence image (true) or persistence weighted gaussian kernel (false)
  // TODO: constexpr ?
  const bool discrete_ = true;
  std::size_t number_of_functions_for_vectorization_;
  std::size_t number_of_functions_for_projections_to_reals_;
  std::function<double(std::pair<double, double>, std::pair<double, double>)> kernel_;
  std::vector<std::pair<double, double> > interval_;
  std::vector<double> weights_;

  // data
  Scaling_of_kernels f_;
  bool erase_below_diagonal_;
};

// if min == max, then the program is requested to set up the values itself based on persistence intervals
template <typename Scaling_of_kernels>
void Persistence_heat_maps<Scaling_of_kernels>::_construct(
    const std::vector<std::pair<double, double> >& intervals,
    const std::vector<std::vector<double> >& filter,
    bool erase_below_diagonal,
    std::size_t number_of_pixels,
    double min,
    double max)
{
#ifdef DEBUG_TRACES
  std::clog << "Entering construct procedure \n";
#endif

  Scaling_of_kernels f;
  this->f_ = f;

#ifdef DEBUG_TRACES
  std::clog << "min and max passed to construct() procedure: " << min << " " << max << std::endl;
#endif

  if (min == max) {
#ifdef DEBUG_TRACES
    std::clog << "min and max parameters will be determined based on intervals \n";
#endif
    // in this case, we want the program to set up the min and max values by itself.
    min = std::numeric_limits<double>::max();
    max = -std::numeric_limits<double>::max();

    for (std::size_t i = 0; i != intervals.size(); ++i) {
      if (intervals[i].first < min) min = intervals[i].first;
      if (intervals[i].second > max) max = intervals[i].second;
    }
    // now we have the structure filled in, and moreover we know min and max values of the interval, so we know the
    // range.

    // add some more space:
    min -= fabs(max - min) / 100;
    max += fabs(max - min) / 100;
  }

#ifdef DEBUG_TRACES
  std::clog << "min_ : " << min << std::endl;
  std::clog << "max_ : " << max << std::endl;
  std::clog << "number_of_pixels : " << number_of_pixels << std::endl;
#endif

  this->min_ = min;
  this->max_ = max;

  // initialization of the structure heat_map_
  std::vector<std::vector<double> > heat_map_;
  for (std::size_t i = 0; i != number_of_pixels; ++i) {
    std::vector<double> v(number_of_pixels, 0);
    heat_map_.push_back(v);
  }
  this->heat_map_ = heat_map_;

#ifdef DEBUG_TRACES
  std::clog << "Done creating of the heat map, now we will fill in the structure \n";
#endif

  for (std::size_t pt_nr = 0; pt_nr != intervals.size(); ++pt_nr) {
    // compute the value of intervals_[pt_nr] in the grid:
    int x_grid =
        static_cast<int>((intervals[pt_nr].first - this->min_) / (this->max_ - this->min_) * number_of_pixels);
    int y_grid =
        static_cast<int>(((intervals[pt_nr].second - intervals[pt_nr].first) - this->min_) / (this->max_ - this->min_) * number_of_pixels);

#ifdef DEBUG_TRACES
    std::clog << "point : " << intervals[pt_nr].first << " , " << intervals[pt_nr].second << std::endl;
    std::clog << "x_grid : " << x_grid << std::endl;
    std::clog << "y_grid : " << y_grid << std::endl;
#endif

    // x_grid and y_grid gives a center of the kernel. We want to have its lower left corner. To get this, we need to
    // shift x_grid and y_grid by a grid diameter.
    x_grid -= filter.size() / 2;
    y_grid -= filter.size() / 2;
    // note that the numbers x_grid and y_grid may be negative.

#ifdef DEBUG_TRACES
    std::clog << "After shift : \n";
    std::clog << "x_grid : " << x_grid << std::endl;
    std::clog << "y_grid : " << y_grid << std::endl;
#endif

    double scaling_value = this->f_(intervals[pt_nr]);

    for (std::size_t i = 0; i != filter.size(); ++i) {
      for (std::size_t j = 0; j != filter.size(); ++j) {
        // if the point (x_grid+i,y_grid+j) is the correct point in the grid.
        if (((x_grid + i) >= 0) && (x_grid + i < this->heat_map_.size()) && ((y_grid + j) >= 0) &&
            (y_grid + j < this->heat_map_.size())) {
#ifdef DEBUG_TRACES
          std::clog << y_grid + j << " " << x_grid + i << std::endl;
#endif
          this->heat_map_[y_grid + j][x_grid + i] += scaling_value * filter[i][j];
#ifdef DEBUG_TRACES
          std::clog << "Position : (" << x_grid + i << "," << y_grid + j
                    << ") got increased by the value : " << filter[i][j] << std::endl;
#endif
        }
      }
    }
  }

  // now it remains to cut everything below diagonal if the user wants us to.
  if (erase_below_diagonal) {
    for (std::size_t i = 0; i != this->heat_map_.size(); ++i) {
      for (std::size_t j = i; j != this->heat_map_.size(); ++j) {
        this->heat_map_[i][j] = 0;
      }
    }
  }
}  // construct

template <typename Scaling_of_kernels>
std::vector<std::vector<double> > Persistence_heat_maps<Scaling_of_kernels>::_check_and_initialize_maps(
    const std::vector<Persistence_heat_maps*>& maps)
{
  // checking if all the heat maps are of the same size:
  GUDHI_CHECK_code(
    for (std::size_t i = 0; i != maps.size(); ++i) {
      GUDHI_CHECK(maps[i]->heat_map_.size() == maps[0]->heat_map_.size(),
                  "Sizes of Persistence_heat_maps are not compatible.");
      GUDHI_CHECK(maps[i]->heat_map_[0].size() == maps[0]->heat_map_[0].size(),
                  "Sizes of Persistence_heat_maps are not compatible.");
    }
  );

  const auto& map = maps[0]->heat_map_;
  return std::vector<std::vector<double> >(map.size(), std::vector<double>(map[0].size(), 0));
}

template <typename Scaling_of_kernels>
void Persistence_heat_maps<Scaling_of_kernels>::compute_median(const std::vector<Persistence_heat_maps*>& maps)
{
  std::vector<std::vector<double> > heat_maps = this->_check_and_initialize_maps(maps);

  std::vector<double> to_compute_median(maps.size());
  for (std::size_t i = 0; i != heat_maps.size(); ++i) {
    for (std::size_t j = 0; j != heat_maps[i].size(); ++j) {
      for (std::size_t map_no = 0; map_no != maps.size(); ++map_no) {
        to_compute_median[map_no] = maps[map_no]->heat_map_[i][j];
      }
      std::nth_element(
          to_compute_median.begin(), to_compute_median.begin() + to_compute_median.size() / 2, to_compute_median.end());
      heat_maps[i][j] = to_compute_median[to_compute_median.size() / 2];
    }
  }
  this->heat_map_ = heat_maps;
  this->min_ = maps[0]->min_;
  this->max_ = maps[0]->max_;
}

template <typename Scaling_of_kernels>
void Persistence_heat_maps<Scaling_of_kernels>::compute_mean(const std::vector<Persistence_heat_maps*>& maps)
{
  std::vector<std::vector<double> > heat_maps = this->_check_and_initialize_maps(maps);
  for (std::size_t i = 0; i != heat_maps.size(); ++i) {
    for (std::size_t j = 0; j != heat_maps[i].size(); ++j) {
      double mean = 0;
      for (std::size_t map_no = 0; map_no != maps.size(); ++map_no) {
        mean += maps[map_no]->heat_map_[i][j];
      }
      heat_maps[i][j] = mean / static_cast<double>(maps.size());
    }
  }
  this->heat_map_ = heat_maps;
  this->min_ = maps[0]->min_;
  this->max_ = maps[0]->max_;
}

template <typename Scaling_of_kernels>
void Persistence_heat_maps<Scaling_of_kernels>::compute_percentage_of_active(
    const std::vector<Persistence_heat_maps*>& maps,
    std::size_t cutoff)
{
  std::vector<std::vector<double> > heat_maps = this->_check_and_initialize_maps(maps);

  for (std::size_t i = 0; i != heat_maps.size(); ++i) {
    for (std::size_t j = 0; j != heat_maps[i].size(); ++j) {
      std::size_t number_of_active_levels = 0;
      for (std::size_t map_no = 0; map_no != maps.size(); ++map_no) {
        if (maps[map_no]->heat_map_[i][j]) number_of_active_levels++;
      }
      if (number_of_active_levels > cutoff) {
        heat_maps[i][j] = number_of_active_levels;
      } else {
        heat_maps[i][j] = 0;
      }
    }
  }
  this->heat_map_ = heat_maps;
  this->min_ = maps[0]->min_;
  this->max_ = maps[0]->max_;
}

template <typename Scaling_of_kernels>
void Persistence_heat_maps<Scaling_of_kernels>::plot(const char* filename) const
{
  std::ofstream out;
  std::stringstream gnuplot_script;
  gnuplot_script << filename << "_GnuplotScript";

  out.open(gnuplot_script.str().c_str());
  out << "plot      '-' matrix with image" << std::endl;
  for (std::size_t i = 0; i != this->heat_map_.size(); ++i) {
    for (std::size_t j = 0; j != this->heat_map_[i].size(); ++j) {
      out << this->heat_map_[i][j] << " ";
    }
    out << std::endl;
  }
  out.close();
#ifdef DEBUG_TRACES
  std::clog << "To visualize, install gnuplot and type the command: gnuplot -persist -e \"load \'"
            << gnuplot_script.str().c_str() << "\'\"" << std::endl;
#endif
}

template <typename Scaling_of_kernels>
void Persistence_heat_maps<Scaling_of_kernels>::print_to_file(const char* filename) const
{
  std::ofstream out;
  out.open(filename);

  // First we store this->min_ and this->max_ values:
  out << this->min_ << " " << this->max_ << std::endl;
  for (std::size_t i = 0; i != this->heat_map_.size(); ++i) {
    for (std::size_t j = 0; j != this->heat_map_[i].size(); ++j) {
      out << this->heat_map_[i][j] << " ";
    }
    out << std::endl;
  }
  out.close();
}

template <typename Scaling_of_kernels>
void Persistence_heat_maps<Scaling_of_kernels>::load_from_file(const char* filename)
{
  std::ifstream in;
  in.open(filename);

  // checking if the file exist / if it was open.
  if (!in.good()) {
#ifdef DEBUG_TRACES
    std::cerr << "The file : " << filename << " do not exist. The program will now terminate \n";
#endif
    throw std::invalid_argument("The persistence landscape file do not exist.");
  }

  // now we read the file one by one.

  in >> this->min_ >> this->max_;
#ifdef DEBUG_TRACES
  std::clog << "Reading the following values of min and max : " << this->min_ << " , " << this->max_ << std::endl;
#endif

  std::string temp;
  std::getline(in, temp);
  while (in.good()) {
    std::getline(in, temp);
    std::stringstream lineSS;
    lineSS << temp;

    std::vector<double> line_of_heat_map;
    while (lineSS.good()) {
      double point;

      lineSS >> point;
      line_of_heat_map.push_back(point);
#ifdef DEBUG_TRACES
      std::clog << point << " ";
#endif
    }
#ifdef DEBUG_TRACES
    std::clog << std::endl;
#endif

    if (in.good()) this->heat_map_.push_back(line_of_heat_map);
  }
  in.close();
#ifdef DEBUG_TRACES
  std::clog << "Done \n";
#endif
}

// Implementation of virtual methods:
template <typename Scaling_of_kernels>
std::vector<double> Persistence_heat_maps<Scaling_of_kernels>::vectorize(int number_of_function) const
{
  std::vector<double> result;

  if (!discrete_) {
#ifdef DEBUG_TRACES
    std::clog << "No vectorize method in case of infinite dimensional vectorization" << std::endl;
#endif
    return result;
  }

  // convert this->heat_map_ into one large vector:
  std::size_t size_of_result = 0;
  for (std::size_t i = 0; i != this->heat_map_.size(); ++i) {
    size_of_result += this->heat_map_[i].size();
  }

  result.reserve(size_of_result);

  for (std::size_t i = 0; i != this->heat_map_.size(); ++i) {
    for (std::size_t j = 0; j != this->heat_map_[i].size(); ++j) {
      result.push_back(this->heat_map_[i][j]);
    }
  }

  return result;
}

template <typename Scaling_of_kernels>
double Persistence_heat_maps<Scaling_of_kernels>::distance(const Persistence_heat_maps& second,
                                                                   double power) const
{
  if (this->discrete_) {
    // first we need to check if (*this) and second are defined on the same domain and have the same dimensions:
    if (!this->check_if_the_same(second)) {
#ifdef DEBUG_TRACES
      std::cerr << "The persistence images are of non compatible sizes. We cannot therefore compute distance between "
                   "them. The program will now terminate";
#endif
      throw std::invalid_argument("The persistence images are of non compatible sizes.");
    }

    // if we are here, we know that the two persistence images are defined on the same domain, so we can start
    // computing their distances:

    double distance = 0;
    if (power < std::numeric_limits<double>::max()) {
      for (std::size_t i = 0; i != this->heat_map_.size(); ++i) {
        for (std::size_t j = 0; j != this->heat_map_[i].size(); ++j) {
          distance += std::pow(std::fabs(this->heat_map_[i][j] - second.heat_map_[i][j]), power);
        }
      }
    } else {
      // in this case, we compute max norm distance
      for (std::size_t i = 0; i != this->heat_map_.size(); ++i) {
        for (std::size_t j = 0; j != this->heat_map_[i].size(); ++j) {
          auto diff = std::fabs(this->heat_map_[i][j] - second.heat_map_[i][j]);
          if (distance < diff) distance = diff;
        }
      }
    }
    return distance;
  } else {
    return std::sqrt(this->compute_scalar_product(*this) + second.compute_scalar_product(second) -
                     2 * this->compute_scalar_product(second));
  }
}

template <typename Scaling_of_kernels>
double Persistence_heat_maps<Scaling_of_kernels>::project_to_R(int number_of_function) const
{
  double result = 0;
  for (std::size_t i = 0; i != this->heat_map_.size(); ++i) {
    for (std::size_t j = 0; j != this->heat_map_[i].size(); ++j) {
      result += this->heat_map_[i][j];
    }
  }
  return result;
}

template <typename Scaling_of_kernels>
void Persistence_heat_maps<Scaling_of_kernels>::compute_average(
    const std::vector<Persistence_heat_maps*>& to_average)
{
  this->compute_mean(to_average);
}

template <typename Scaling_of_kernels>
double Persistence_heat_maps<Scaling_of_kernels>::compute_scalar_product(
    const Persistence_heat_maps& second) const
{
  if (discrete_) {
    // first we need to check if (*this) and second are defined on the same domain and have the same dimensions:
    if (!this->check_if_the_same(second)) {
#ifdef DEBUG_TRACES
      std::cerr << "The persistence images are of non compatible sizes. We cannot therefore compute distance between "
                   "them. The program will now terminate";
#endif
      throw std::invalid_argument("The persistence images are of non compatible sizes.");
    }

    // if we are here, we know that the two persistence images are defined on the same domain, so we can start computing
    // their scalar product:
    double scalar_prod = 0;
    for (std::size_t i = 0; i != this->heat_map_.size(); ++i) {
      for (std::size_t j = 0; j != this->heat_map_[i].size(); ++j) {
        scalar_prod += this->heat_map_[i][j] * second.heat_map_[i][j];
      }
    }
    return scalar_prod;
  } else {
    int num_pts1 = this->interval_.size();
    int num_pts2 = second.interval_.size();
    double kernel_val = 0;
    for (int i = 0; i < num_pts1; i++) {
      std::pair<double, double> pi = this->interval_[i];
      for (int j = 0; j < num_pts2; j++) {
        std::pair<double, double> pj = second.interval_[j];
        kernel_val += this->weights_[i] * second.weights_[j] * this->kernel_(pi, pj);
      }
    }
    return kernel_val;
  }
}

}  // namespace Persistence_representations
}  // namespace Gudhi

#endif  // PERSISTENCE_HEAT_MAPS_H_
