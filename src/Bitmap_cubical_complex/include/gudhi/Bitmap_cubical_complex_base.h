/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2015 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef BITMAP_CUBICAL_COMPLEX_BASE_H_
#define BITMAP_CUBICAL_COMPLEX_BASE_H_

#include <gudhi/Debug_utils.h>

#include <boost/config.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/iterator_range.hpp>

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <limits>
#include <utility>
#include <stdexcept>
#include <cstddef>
#include <numeric>
#include <functional>

namespace Gudhi {

namespace cubical_complex {

/**
 * @brief Cubical complex represented as a bitmap, class with basic implementation.
 * @ingroup cubical_complex
 * @details This is a class implementing a basic bitmap data structure to store cubical complexes.
 * It implements only the most basic subroutines.
 * The idea of the bitmap is the following. Our aim is to have a memory efficient
 * data structure to store d-dimensional cubical complex
 * C being a cubical decomposition
 * of a rectangular region of a space. This is achieved by storing C as a
 * vector of bits (this is where the name 'bitmap' came from).
 * Each cell is represented by a single
 * bit (in case of black and white bitmaps, or by a single element of a type T
 * (here T is a filtration type of a bitmap, typically a double).
 * All the information needed for homology and
 * persistent homology computations (like dimension of a cell, boundary and
 * coboundary elements of a cell, are then obtained from the
 * position of the element in C.
 * The default filtration used in this implementation is the lower star filtration.
 */
template <typename T>
class Bitmap_cubical_complex_base {
 public:
  typedef T filtration_type;

  /**
   *Default constructor
   **/
  Bitmap_cubical_complex_base() {}
  /**
   * There are a few constructors of a Bitmap_cubical_complex_base class.
   * First one, that takes vector<unsigned>, creates an empty bitmap of a dimension equal
   * the number of elements in the
   * input vector and size in the i-th dimension equal the number in the position i-of the input vector.
   */
  explicit Bitmap_cubical_complex_base(const std::vector<unsigned>& sizes);
  /**
   * The second constructor takes as a input a Perseus style file. For more details,
   * please consult the documentations of
   * Perseus software as well as examples attached to this
   * implementation.
   **/
  explicit Bitmap_cubical_complex_base(const char* perseus_style_file);
  /**
   * The last constructor of a Bitmap_cubical_complex_base class accepts vector of dimensions (as the first one)
   * with vector of filtration values of vertices or top dimensional cells depending on the input_top_cells flag.
   **/
  Bitmap_cubical_complex_base(const std::vector<unsigned>& dimensions, const std::vector<T>& cells, bool input_top_cells = true);

  /**
   * Destructor of the Bitmap_cubical_complex_base class.
   **/
  virtual ~Bitmap_cubical_complex_base() {}

  /**
   * The functions get_boundary_of_a_cell, get_coboundary_of_a_cell, get_dimension_of_a_cell
   * and get_cell_data are the basic
   * functions that compute boundary / coboundary / dimension and the filtration
   * value form a position of a cell in the structure of a bitmap. The input parameter of all of those function is a
   * non-negative integer, indicating a position of a cube in the data structure.
   * In the case of functions that compute (co)boundary, the output is a vector if non-negative integers pointing to
   * the positions of (co)boundary element of the input cell.
   * The boundary elements are guaranteed to be returned so that the
   * incidence coefficients of boundary elements are alternating.
   */
  virtual inline std::vector<std::size_t> get_boundary_of_a_cell(std::size_t cell) const;
  /**
   * The functions get_coboundary_of_a_cell, get_coboundary_of_a_cell,
   * get_dimension_of_a_cell and get_cell_data are the basic
   * functions that compute boundary / coboundary / dimension and the filtration
   * value form a position of a cell in the structure of a bitmap.
   * The input parameter of all of those function is a non-negative integer,
   * indicating a position of a cube in the data structure.
   * In the case of functions that compute (co)boundary, the output is a vector if
   * non-negative integers pointing to the
   * positions of (co)boundary element of the input cell.
   * Note that unlike in the case of boundary, over here the elements are
   * not guaranteed to be returned with alternating incidence numbers.
   *
   **/
  virtual inline std::vector<std::size_t> get_coboundary_of_a_cell(std::size_t cell) const;

  /**
   * This function finds a top-dimensional cell that is incident to the input cell and has
   * the same filtration value. In case several cells are suitable, an arbitrary one is
   * returned. Note that the input parameter can be a cell of any dimension (vertex, edge, etc).
   * On the other hand, the output is always indicating the position of
   * a top-dimensional cube in the data structure.
   * \pre The filtration values are assigned as per `impose_lower_star_filtration()`.
   **/
  inline size_t get_top_dimensional_coface_of_a_cell(size_t splx);

  /**
   * This function finds a vertex that is incident to the input cell and has
   * the same filtration value. In case several cells are suitable, an arbitrary one is
   * returned. Note that the input parameter can be a cell of any dimension (vertex, edge, etc).
   * On the other hand, the output is always indicating the position of
   * a vertex in the data structure.
   * \pre The filtration values are assigned as per `impose_lower_star_filtration_from_vertices()`.
   **/
  inline size_t get_vertex_of_a_cell(size_t splx);

  /**
  * This procedure compute incidence numbers between cubes. For a cube \f$A\f$ of
  * dimension n and a cube \f$B \subset A\f$ of dimension n-1, an incidence
  * between \f$A\f$ and \f$B\f$ is the integer with which \f$B\f$ appears in the boundary of \f$A\f$.
  * Note that first parameter is a cube of dimension n,
  * and the second parameter is an adjusted cube in dimension n-1.
  * Given \f$A = [b_1,e_1] \times \ldots \ [b_{j-1},e_{j-1}] \times [b_{j},e_{j}] \times [b_{j+1},e_{j+1}] \times \ldots
  *\times [b_{n},e_{n}] \f$
  * such that \f$ b_{j} \neq e_{j} \f$
  * and \f$B = [b_1,e_1] \times \ldots \ [b_{j-1},e_{j-1}] \times [a,a] \times [b_{j+1},e_{j+1}] \times \ldots \times
  *[b_{n},e_{n}] \f$
  * where \f$ a = b_{j}\f$ or \f$ a = e_{j}\f$, the incidence between \f$A\f$ and \f$B\f$
  * computed by this procedure is given by formula:
  * \f$ c\ (-1)^{\sum_{i=1}^{j-1} dim [b_{i},e_{i}]}  \f$
  * Where \f$ dim [b_{i},e_{i}] = 0 \f$ if \f$ b_{i}=e_{i} \f$ and 1 in other case.
  * c is -1 if \f$ a = b_{j}\f$ and 1 if \f$ a = e_{j}\f$.
  * @exception std::logic_error In case when the cube \f$B\f$ is not n-1
  * dimensional face of a cube \f$A\f$.
  **/
  virtual int compute_incidence_between_cells(std::size_t coface, std::size_t face) const {
    // first get the counters for coface and face:
    std::vector<unsigned> coface_counter = this->compute_counter_for_given_cell(coface);
    std::vector<unsigned> face_counter = this->compute_counter_for_given_cell(face);

    // coface_counter and face_counter should agree at all positions except from one:
    int number_of_position_in_which_counters_do_not_agree = -1;
    std::size_t number_of_full_faces_that_comes_before = 0;
    for (std::size_t i = 0; i != coface_counter.size(); ++i) {
      if ((coface_counter[i] % 2 == 1) && (number_of_position_in_which_counters_do_not_agree == -1)) {
        ++number_of_full_faces_that_comes_before;
      }
      if (coface_counter[i] != face_counter[i]) {
        if (number_of_position_in_which_counters_do_not_agree != -1) {
          std::cerr << "Cells given to compute_incidence_between_cells procedure do not form a pair of coface-face.\n";
          throw std::logic_error(
              "Cells given to compute_incidence_between_cells procedure do not form a pair of coface-face.");
        }
        number_of_position_in_which_counters_do_not_agree = i;
      }
    }

    int incidence = 1;
    if (number_of_full_faces_that_comes_before % 2) incidence = -1;
    // if the face cell is on the right from coface cell:
    if (coface_counter[number_of_position_in_which_counters_do_not_agree] + 1 ==
        face_counter[number_of_position_in_which_counters_do_not_agree]) {
      incidence *= -1;
    }

    return incidence;
  }

  /**
* In the case of get_dimension_of_a_cell function, the output is a non-negative integer
* indicating the dimension of a cell.
* Note that unlike in the case of boundary, over here the elements are
* not guaranteed to be returned with alternating incidence numbers.
* To compute incidence between cells use compute_incidence_between_cells
* procedure
**/
  inline unsigned get_dimension_of_a_cell(std::size_t cell) const;

  /**
   * In the case of get_cell_data, the output parameter is a reference to the value of a cube in a given position.
   * This allows reading and changing the value of filtration. Note that if the value of a filtration is changed, the
   * code do not check if we have a filtration or not. i.e. it do not check if the value of a filtration of a cell is
   * not smaller than the value of a filtration of its boundary and not greater than the value of its coboundary.
   **/
  inline T& get_cell_data(std::size_t cell);

  /**
   * Typical input used to construct a baseBitmap class is a filtration given at the top dimensional cells.
   * Then, there are a few ways one can pick the filtration of lower dimensional
   * cells. The most typical one is by so called lower star filtration. This function is always called by any
   * constructor which takes the top dimensional cells. If you use such a constructor,
   * then there is no need to call this function. Call it only if you are putting the filtration
   * of the cells by your own (for instance by using Top_dimensional_cells_iterator).
   **/
  void impose_lower_star_filtration();  // assume that top dimensional cells are already set.

  /**
   * Set cells filtrations given those of the vertices, and based on lower star filtration.
   * This is already called by the relevant constructors.
   **/
  void impose_lower_star_filtration_from_vertices();  // assume that vertices are already set.

  /**
   * Returns dimension of a complex.
   **/
  inline unsigned dimension() const { return sizes.size(); }

  /**
   * Returns number of all cubes in the data structure.
   **/
  inline std::size_t size() const { return this->data.size(); }

  /**
   * Function that put the input data to bins. By putting data to bins we mean rounding them to a sequence of values
   * equally distributed in the range of data.
   * Sometimes if most of the cells have different birth-death times, the performance of the algorithms to compute
   * persistence gets worst. When dealing with this type of data, one may want to put different values on cells to
   * some number of bins. The function put_data_to_bins( std::size_t number_of_bins ) is designed for that purpose.
   * The parameter of the function is the number of bins (distinct values) we want to have in the cubical complex.
   **/
  void put_data_to_bins(std::size_t number_of_bins);

  /**
   * Function that put the input data to bins. By putting data to bins we mean rounding them to a sequence of values
   * equally distributed in the range of data.
   * Sometimes if most of the cells have different birth-death times, the performance of the algorithms to compute
   * persistence gets worst. When dealing with this type of data, one may want to put different values on cells to
   * some number of bins. The function put_data_to_bins( T diameter_of_bin ) is designed for that purpose.
   * The parameter of it is the diameter of each bin. Note that the bottleneck distance between the persistence
   * diagram of the cubical complex before and after using such a function will be bounded by the parameter
   * diameter_of_bin.
   **/
  void put_data_to_bins(T diameter_of_bin);

  /**
   * Functions to find min and max values of filtration.
   **/
  std::pair<T, T> min_max_filtration();

  // ITERATORS

  /**
   * @brief Iterator through all cells in the complex (in order they appear in the structure -- i.e.
   * in lexicographical order).
   **/
  typedef boost::counting_iterator<std::size_t> All_cells_iterator;

  /**
   * Function returning a All_cells_iterator to the first cell of the bitmap.
   **/
  All_cells_iterator all_cells_iterator_begin() const { return All_cells_iterator(0); }

  /**
   * Function returning a All_cells_iterator beyond the last cell of the bitmap.
   **/
  All_cells_iterator all_cells_iterator_end() const { return All_cells_iterator(data.size()); }

  /**
   * @brief Range corresponding to All_cells_iterator
   **/
  typedef boost::iterator_range<All_cells_iterator> All_cells_range;

  /** Returns a range over all cells. */
  All_cells_range all_cells_range() const { return All_cells_range(all_cells_iterator_begin(), all_cells_iterator_end()); }

  /**
   * Boundary_range class provides ranges for boundary iterators.
   **/
  typedef typename std::vector<std::size_t>::const_iterator Boundary_iterator;
  typedef typename std::vector<std::size_t> Boundary_range;

  /**
   * boundary_simplex_range creates an object of a Boundary_simplex_range class
   * that provides ranges for the Boundary_simplex_iterator.
   **/
  Boundary_range boundary_range(std::size_t sh) { return this->get_boundary_of_a_cell(sh); }

  /**
   * Coboundary_range class provides ranges for boundary iterators.
   **/
  typedef typename std::vector<std::size_t>::const_iterator Coboundary_iterator;
  typedef typename std::vector<std::size_t> Coboundary_range;

  /**
   * boundary_simplex_range creates an object of a Boundary_simplex_range class
   * that provides ranges for the Boundary_simplex_iterator.
   **/
  Coboundary_range coboundary_range(std::size_t sh) { return this->get_coboundary_of_a_cell(sh); }

  /**
   * @brief Iterator through top dimensional cells of the complex. The cells appear in order they are stored
   * in the structure (i.e. in lexicographical order)
   **/
  class Top_dimensional_cells_iterator {
   public:
    typedef std::input_iterator_tag iterator_category;
    typedef std::size_t value_type;
    typedef std::ptrdiff_t difference_type;
    typedef value_type* pointer;
    typedef value_type reference;

    Top_dimensional_cells_iterator(Bitmap_cubical_complex_base* b) : counter(b->dimension()), b(b) {}

    Top_dimensional_cells_iterator operator++() {
      // first find first element of the counter that can be increased:
      std::size_t dim = 0;
      while ((dim != this->b->dimension()) && (this->counter[dim] == this->b->sizes[dim] - 1)) ++dim;

      if (dim != this->b->dimension()) {
        ++this->counter[dim];
        for (std::size_t i = 0; i != dim; ++i) {
          this->counter[i] = 0;
        }
      } else {
        ++this->counter[0];
      }
      return *this;
    }

    Top_dimensional_cells_iterator operator++(int) {
      Top_dimensional_cells_iterator result = *this;
      ++(*this);
      return result;
    }

    bool operator==(const Top_dimensional_cells_iterator& rhs) const {
      if (this->b != rhs.b) return false;
      if (this->counter.size() != rhs.counter.size()) return false;
      for (std::size_t i = 0; i != this->counter.size(); ++i) {
        if (this->counter[i] != rhs.counter[i]) return false;
      }
      return true;
    }

    bool operator!=(const Top_dimensional_cells_iterator& rhs) const { return !(*this == rhs); }

    /*
     * The operator * returns position of a cube in the structure of cubical complex. This position can be then used as
     * an argument of the following functions:
     * get_boundary_of_a_cell, get_coboundary_of_a_cell, get_dimension_of_a_cell to get information about the cell
     * boundary and coboundary and dimension
     * and in function get_cell_data to get a filtration of a cell.
     */
    std::size_t operator*() { return this->compute_index_in_bitmap(); }

    std::size_t compute_index_in_bitmap() const {
      std::size_t index = 0;
      for (std::size_t i = 0; i != this->counter.size(); ++i) {
        index += (2 * this->counter[i] + 1) * this->b->multipliers[i];
      }
      return index;
    }

    void print_counter() const {
      for (std::size_t i = 0; i != this->counter.size(); ++i) {
        std::clog << this->counter[i] << " ";
      }
    }
    friend class Bitmap_cubical_complex_base;

   protected:
    std::vector<std::size_t> counter;
    Bitmap_cubical_complex_base* b;
  };

  /**
   * Function returning a Top_dimensional_cells_iterator to the first top dimensional cell of the bitmap.
   **/
  Top_dimensional_cells_iterator top_dimensional_cells_iterator_begin() {
    Top_dimensional_cells_iterator a(this);
    return a;
  }

  /**
   * Function returning a Top_dimensional_cells_iterator beyond the last top dimensional cell of the bitmap.
   **/
  Top_dimensional_cells_iterator top_dimensional_cells_iterator_end() {
    Top_dimensional_cells_iterator a(this);
    for (std::size_t i = 0; i != this->dimension(); ++i) {
      a.counter[i] = this->sizes[i] - 1;
    }
    a.counter[0]++;
    return a;
  }

  /**
   * @brief Range corresponding to Top_dimensional_cells_iterator
   **/
  class Top_dimensional_cells_range {
   public:
    Top_dimensional_cells_range(Bitmap_cubical_complex_base* b) : b(b) {}

    Top_dimensional_cells_iterator begin() { return b->top_dimensional_cells_iterator_begin(); }

    Top_dimensional_cells_iterator end() { return b->top_dimensional_cells_iterator_end(); }

   private:
    Bitmap_cubical_complex_base<T>* b;
  };

  /** Returns a range over all top-dimensional cells. */
  Top_dimensional_cells_range top_dimensional_cells_range() { return Top_dimensional_cells_range(this); }

  //****************************************************************************************************************//
  //****************************************************************************************************************//
  //****************************************************************************************************************//
  //****************************************************************************************************************//
  class Vertices_iterator {
   public:
    typedef std::input_iterator_tag iterator_category;
    typedef std::size_t value_type;
    typedef std::ptrdiff_t difference_type;
    typedef value_type* pointer;
    typedef value_type reference;

    Vertices_iterator(Bitmap_cubical_complex_base* b) : counter(b->dimension()), b(b) {}

    Vertices_iterator operator++() {
      // first find first element of the counter that can be increased:
      std::size_t dim = 0;
      while ((dim != this->b->dimension()) && (this->counter[dim] == this->b->sizes[dim])) ++dim;

      if (dim != this->b->dimension()) {
        ++this->counter[dim];
        for (std::size_t i = 0; i != dim; ++i) {
          this->counter[i] = 0;
        }
      } else {
        ++this->counter[0];
      }
      return *this;
    }

    Vertices_iterator operator++(int) {
      Vertices_iterator result = *this;
      ++(*this);
      return result;
    }

    bool operator==(const Vertices_iterator& rhs) const {
      if (this->b != rhs.b) return false;
      GUDHI_CHECK(this->counter.size() == rhs.counter.size(), "impossible");
      for (std::size_t i = 0; i != this->counter.size(); ++i) {
        if (this->counter[i] != rhs.counter[i]) return false;
      }
      return true;
    }

    bool operator!=(const Vertices_iterator& rhs) const { return !(*this == rhs); }

    /*
     * The operator * returns position of a cube in the structure of cubical complex. This position can be then used as
     * an argument of the following functions:
     * get_boundary_of_a_cell, get_coboundary_of_a_cell, get_dimension_of_a_cell to get information about the cell
     * boundary and coboundary and dimension
     * and in function get_cell_data to get a filtration of a cell.
     */
    std::size_t operator*() const { return this->compute_index_in_bitmap(); }

    std::size_t compute_index_in_bitmap() const {
      std::size_t index = 0;
      for (std::size_t i = 0; i != this->counter.size(); ++i) {
        index += 2 * this->counter[i] * this->b->multipliers[i];
      }
      return index;
    }

    void print_counter() const {
      for (std::size_t i = 0; i != this->counter.size(); ++i) {
        std::clog << this->counter[i] << " ";
      }
    }
    friend class Bitmap_cubical_complex_base;

   protected:
    std::vector<std::size_t> counter;
    Bitmap_cubical_complex_base* b;
  };

  /*
   * Function returning a Vertices_iterator to the first vertex of the bitmap.
   */
  Vertices_iterator vertices_iterator_begin() {
    Vertices_iterator a(this);
    return a;
  }

  /*
   * Function returning a Vertices_iterator to the last vertex of the bitmap.
   */
  Vertices_iterator vertices_iterator_end() {
    Vertices_iterator a(this);
    for (std::size_t i = 0; i != this->dimension(); ++i) {
      a.counter[i] = this->sizes[i];
    }
    a.counter[0]++;
    return a;
  }

  /*
   * @brief Vertices_iterator_range class provides ranges for Vertices_iterator_range
   */
  class Vertices_range {
   public:
    Vertices_range(Bitmap_cubical_complex_base* b) : b(b) {}

    Vertices_iterator begin() { return b->vertices_iterator_begin(); }

    Vertices_iterator end() { return b->vertices_iterator_end(); }

   private:
    Bitmap_cubical_complex_base* b;
  };

  /* Returns a range over all vertices. */
  Vertices_range vertices_range() { return Vertices_range(this); }

  //****************************************************************************************************************//
  //****************************************************************************************************************//
  //****************************************************************************************************************//
  //****************************************************************************************************************//

  inline std::size_t number_cells() const { return this->data.size(); }

  //****************************************************************************************************************//
  //****************************************************************************************************************//
  //****************************************************************************************************************//
  //****************************************************************************************************************//

  /// @private @brief execute the function on each cell id corresponding to a vertex, in increasing order.
  template <class F> void for_each_vertex(F&&f) {
    for_each_vertex_rec(f, 0, multipliers.size()-1);
  }

 protected:
  std::vector<unsigned> sizes;
  std::vector<unsigned> multipliers;
  std::vector<T> data;

  template <class F> void for_each_vertex_rec(F&&f, std::size_t base, int dim);
  void propagate_from_vertices_rec(int special_dim, int current_dim, std::size_t base);

  void set_up_containers(const std::vector<unsigned>& sizes, bool is_pos_inf) {
    // The fact that multipliers[0]=1 is relied on by optimizations in other functions
    unsigned multiplier = 1;
    for (std::size_t i = 0; i != sizes.size(); ++i) {
      this->sizes.push_back(sizes[i]);
      this->multipliers.push_back(multiplier);
      multiplier *= 2 * sizes[i] + 1;
    }
    if(is_pos_inf)
      this->data = std::vector<T>(multiplier, std::numeric_limits<T>::infinity());
    else
      this->data = std::vector<T>(multiplier, -std::numeric_limits<T>::infinity());
  }

  std::size_t compute_position_in_bitmap(const std::vector<unsigned>& counter) {
    std::size_t position = 0;
    for (std::size_t i = 0; i != this->multipliers.size(); ++i) {
      position += this->multipliers[i] * counter[i];
    }
    return position;
  }

  std::vector<unsigned> compute_counter_for_given_cell(std::size_t cell) const {
    std::vector<unsigned> counter;
    counter.reserve(this->sizes.size());
    for (std::size_t dim = this->sizes.size(); dim > 1; --dim) {
      std::size_t quot = cell / this->multipliers[dim - 1];
      cell = cell % this->multipliers[dim - 1];
      counter.push_back(quot);
    }
    // Split out the last iteration to save a costly division by multipliers[0]=1
    counter.push_back(cell);
    std::reverse(counter.begin(), counter.end());
    return counter;
  }
  void read_perseus_style_file(const char* perseus_style_file);
  void setup_bitmap_based_on_top_dimensional_cells_list(const std::vector<unsigned>& sizes_in_following_directions,
                                                        const std::vector<T>& top_dimensional_cells);
  void setup_bitmap_based_on_vertices(const std::vector<unsigned>& sizes_in_following_directions,
                                      const std::vector<T>& vertices);
  Bitmap_cubical_complex_base(const char* perseus_style_file, const std::vector<bool>& directions);
  Bitmap_cubical_complex_base(const std::vector<unsigned>& sizes, const std::vector<bool>& directions);
  Bitmap_cubical_complex_base(const std::vector<unsigned>& dimensions, const std::vector<T>& cells,
                              const std::vector<bool>& directions, bool input_top_cells);
};

template <typename T>
void Bitmap_cubical_complex_base<T>::put_data_to_bins(std::size_t number_of_bins) {

  std::pair<T, T> min_max = this->min_max_filtration();
  T dx = (min_max.second - min_max.first) / (T)number_of_bins;

  // now put the data into the appropriate bins:
  for (std::size_t i = 0; i != this->data.size(); ++i) {
#ifdef DEBUG_TRACES
    std::clog << "Before binning : " << this->data[i] << std::endl;
#endif
    this->data[i] = min_max.first + dx * (this->data[i] - min_max.first) / number_of_bins;
#ifdef DEBUG_TRACES
    std::clog << "After binning : " << this->data[i] << std::endl;
#endif
  }
}

template <typename T>
void Bitmap_cubical_complex_base<T>::put_data_to_bins(T diameter_of_bin) {
  std::pair<T, T> min_max = this->min_max_filtration();

  std::size_t number_of_bins = (min_max.second - min_max.first) / diameter_of_bin;
  // now put the data into the appropriate bins:
  for (std::size_t i = 0; i != this->data.size(); ++i) {
#ifdef DEBUG_TRACES
    std::clog << "Before binning : " << this->data[i] << std::endl;
#endif
    this->data[i] = min_max.first + diameter_of_bin * (this->data[i] - min_max.first) / number_of_bins;
#ifdef DEBUG_TRACES
    std::clog << "After binning : " << this->data[i] << std::endl;
#endif
  }
}

template <typename T>
std::pair<T, T> Bitmap_cubical_complex_base<T>::min_max_filtration() {
  std::pair<T, T> min_max(std::numeric_limits<T>::infinity(), -std::numeric_limits<T>::infinity());
  for (std::size_t i = 0; i != this->data.size(); ++i) {
    if (this->data[i] < min_max.first) min_max.first = this->data[i];
    if (this->data[i] > min_max.second) min_max.second = this->data[i];
  }
  return min_max;
}

template <typename T>
Bitmap_cubical_complex_base<T>::Bitmap_cubical_complex_base(const std::vector<unsigned>& sizes) {
  this->set_up_containers(sizes, true);
}

template <typename T>
void Bitmap_cubical_complex_base<T>::setup_bitmap_based_on_top_dimensional_cells_list(
    const std::vector<unsigned>& sizes_in_following_directions, const std::vector<T>& top_dimensional_cells) {
  this->set_up_containers(sizes_in_following_directions, true);
  std::size_t number_of_top_dimensional_elements = std::accumulate(std::begin(sizes_in_following_directions),
                                                   std::end(sizes_in_following_directions), std::size_t(1),
                                                   std::multiplies<std::size_t>());
  if (number_of_top_dimensional_elements != top_dimensional_cells.size()) {
    std::cerr << "Error in constructor Bitmap_cubical_complex_base ( std::vector<unsigned> "
              << "sizes_in_following_directions, std::vector<T> top_dimensional_cells ). Number of top dimensional "
              << "elements that follow from sizes_in_following_directions vector is different from the size of "
              << "top_dimensional_cells vector."
              << std::endl;
    throw std::invalid_argument(
        "Error in constructor Bitmap_cubical_complex_base( std::vector<unsigned> sizes_in_following_directions,"
        "std::vector<T> top_dimensional_cells ). Number of top dimensional elements that follow from "
        "sizes_in_following_directions vector is different from the size of top_dimensional_cells vector.");
  }

  std::size_t index = 0;
  for (auto it = this->top_dimensional_cells_iterator_begin();
       it != this->top_dimensional_cells_iterator_end(); ++it) {
    this->get_cell_data(*it) = top_dimensional_cells[index];
    ++index;
  }
  this->impose_lower_star_filtration();
}

template <typename T>
void Bitmap_cubical_complex_base<T>::setup_bitmap_based_on_vertices(const std::vector<unsigned>& sizes_in_following_directions,
                                                                    const std::vector<T>& vertices) {
  std::vector<unsigned> top_cells_sizes;
  std::transform (sizes_in_following_directions.begin(), sizes_in_following_directions.end(), std::back_inserter(top_cells_sizes),
               [](int i){ return i-1;});
  this->set_up_containers(top_cells_sizes, false);
  std::size_t number_of_vertices = std::accumulate(std::begin(sizes_in_following_directions),
                                                   std::end(sizes_in_following_directions), (std::size_t)1,
                                                   std::multiplies<std::size_t>());
  if (number_of_vertices != vertices.size()) {
    std::cerr << "Error in constructor Bitmap_cubical_complex_base ( std::vector<unsigned> "
              << "sizes_in_following_directions, std::vector<T> vertices ). Number of vertices "
              << "elements that follow from sizes_in_following_directions vector is different from the size of "
              << "vertices vector."
              << std::endl;
    throw std::invalid_argument(
        "Error in constructor Bitmap_cubical_complex_base( std::vector<unsigned> sizes_in_following_directions,"
        "std::vector<T> vertices ). Number of vertices elements that follow from "
        "sizes_in_following_directions vector is different from the size of vertices vector.");
  }

  for_each_vertex([this, &vertices, index=(std::size_t)0] (auto cell) mutable { get_cell_data(cell) = vertices[index++]; });
  this->impose_lower_star_filtration_from_vertices();
}

template <typename T>
size_t Bitmap_cubical_complex_base<T>::get_top_dimensional_coface_of_a_cell(size_t splx) {
  if (this->get_dimension_of_a_cell(splx) == this->dimension()){return splx;}
  else {
    for (auto v : this->get_coboundary_of_a_cell(splx)){
      if(this->get_cell_data(v) == this->get_cell_data(splx)){
        return this->get_top_dimensional_coface_of_a_cell(v);
      }
    }
  }
  BOOST_UNREACHABLE_RETURN(-2);
}

template <typename T>
size_t Bitmap_cubical_complex_base<T>::get_vertex_of_a_cell(size_t splx) {
  if (this->get_dimension_of_a_cell(splx) == 0){return splx;}
  else {
    for (auto v : this->get_boundary_of_a_cell(splx)){
      if(this->get_cell_data(v) == this->get_cell_data(splx)){
        return this->get_vertex_of_a_cell(v);
      }
    }
  }
  BOOST_UNREACHABLE_RETURN(-2);
}

template <typename T>
Bitmap_cubical_complex_base<T>::Bitmap_cubical_complex_base(const std::vector<unsigned>& sizes_in_following_directions,
                                                            const std::vector<T>& cells, bool input_top_cells) {
  if (input_top_cells) {
    this->setup_bitmap_based_on_top_dimensional_cells_list(sizes_in_following_directions, cells);
  } else {
    this->setup_bitmap_based_on_vertices(sizes_in_following_directions, cells);
  }
}

template <typename T>
void Bitmap_cubical_complex_base<T>::read_perseus_style_file(const char* perseus_style_file) {
  std::ifstream inFiltration(perseus_style_file);
  if(!inFiltration) throw std::ios_base::failure(std::string("Could not open the file ") + perseus_style_file);
  unsigned dimensionOfData;
  inFiltration >> dimensionOfData;

#ifdef DEBUG_TRACES
  std::clog << "dimensionOfData : " << dimensionOfData << std::endl;
#endif

  std::vector<unsigned> sizes;
  sizes.reserve(dimensionOfData);
  // all dimensions multiplied
  std::size_t dimensions = 1;
  for (std::size_t i = 0; i != dimensionOfData; ++i) {
    unsigned size_in_this_dimension;
    inFiltration >> size_in_this_dimension;
    sizes.push_back(size_in_this_dimension);
    dimensions *= size_in_this_dimension;
#ifdef DEBUG_TRACES
    std::clog << "size_in_this_dimension : " << size_in_this_dimension << std::endl;
#endif
  }
  this->set_up_containers(sizes, true);

  Bitmap_cubical_complex_base<T>::Top_dimensional_cells_iterator it = this->top_dimensional_cells_iterator_begin();

  T filtrationLevel = 0.;
  std::size_t filtration_counter = 0;
  while (!inFiltration.eof()) {
    std::string line;
    getline(inFiltration, line);
    if (line.length() != 0) {
      int n = sscanf(line.c_str(), "%lf", &filtrationLevel);
      if (n != 1) {
        std::string perseus_error("Bad Perseus file format. This line is incorrect : " + line);
        throw std::ios_base::failure(perseus_error.c_str());
      }

#ifdef DEBUG_TRACES
      std::clog << "Cell of an index : " << it.compute_index_in_bitmap()
                << " and dimension: " << this->get_dimension_of_a_cell(it.compute_index_in_bitmap())
                << " get the value : " << filtrationLevel << std::endl;
#endif
      this->get_cell_data(*it) = filtrationLevel;
      ++it;
      ++filtration_counter;
    }
  }

  if (filtration_counter != dimensions) {
    std::string perseus_error("Bad Perseus file format. Read " + std::to_string(filtration_counter) + " expected " + \
      std::to_string(dimensions) + " values");
    throw std::ios_base::failure(perseus_error);
  }

  inFiltration.close();
  this->impose_lower_star_filtration();
}

template <typename T>
Bitmap_cubical_complex_base<T>::Bitmap_cubical_complex_base(const char* perseus_style_file,
                                                            const std::vector<bool>& directions) {
  // this constructor is here just for compatibility with a class that creates cubical complexes with periodic boundary
  // conditions.
  // It ignores the last parameter of the function.
  this->read_perseus_style_file(perseus_style_file);
}

template <typename T>
Bitmap_cubical_complex_base<T>::Bitmap_cubical_complex_base(const std::vector<unsigned>& sizes,
                                                            const std::vector<bool>& directions) {
  // this constructor is here just for compatibility with a class that creates cubical complexes with periodic boundary
  // conditions.
  // It ignores the last parameter of the function.
  this->set_up_containers(sizes, true);
}

template <typename T>
Bitmap_cubical_complex_base<T>::Bitmap_cubical_complex_base(const std::vector<unsigned>& dimensions,
                                                            const std::vector<T>& cells,
                                                            const std::vector<bool>& directions,
                                                            bool input_top_cells) {
  // this constructor is here just for compatibility with a class that creates cubical complexes with periodic boundary
  // conditions.
  // It ignores the last parameter of the function.
  if (input_top_cells) {
    this->setup_bitmap_based_on_top_dimensional_cells_list(dimensions, cells);
  } else {
    this->setup_bitmap_based_on_vertices(dimensions, cells);
  }
}

template <typename T>
Bitmap_cubical_complex_base<T>::Bitmap_cubical_complex_base(const char* perseus_style_file) {
  this->read_perseus_style_file(perseus_style_file);
}

template <typename T>
std::vector<std::size_t> Bitmap_cubical_complex_base<T>::get_boundary_of_a_cell(std::size_t cell) const {
  std::vector<std::size_t> boundary_elements;

  // Speed traded of for memory. Check if it is better in practice.
  boundary_elements.reserve(this->dimension() * 2);

  std::size_t sum_of_dimensions = 0;
  std::size_t cell1 = cell;
  for (std::size_t i = this->multipliers.size(); i > 1; --i) {
    unsigned position = cell1 / this->multipliers[i - 1];
    cell1 = cell1 % this->multipliers[i - 1];
    if (position % 2 == 1) {
      if (sum_of_dimensions % 2) {
        boundary_elements.push_back(cell + this->multipliers[i - 1]);
        boundary_elements.push_back(cell - this->multipliers[i - 1]);
      } else {
        boundary_elements.push_back(cell - this->multipliers[i - 1]);
        boundary_elements.push_back(cell + this->multipliers[i - 1]);
      }
      ++sum_of_dimensions;
    }
  }
  // Split out the last iteration to save a costly division by multipliers[0]=1
  if (cell1 % 2 == 1) {
    if (sum_of_dimensions % 2) {
      boundary_elements.push_back(cell + 1);
      boundary_elements.push_back(cell - 1);
    } else {
      boundary_elements.push_back(cell - 1);
      boundary_elements.push_back(cell + 1);
    }
    ++sum_of_dimensions;
  }

  return boundary_elements;
}

template <typename T>
std::vector<std::size_t> Bitmap_cubical_complex_base<T>::get_coboundary_of_a_cell(std::size_t cell) const {
  std::vector<unsigned> counter = this->compute_counter_for_given_cell(cell);
  std::vector<std::size_t> coboundary_elements;
  std::size_t cell1 = cell;
  for (std::size_t i = this->multipliers.size(); i > 1; --i) {
    // It is a bit sad to recompute those divisions when we just did them in compute_counter_for_given_cell.
    unsigned position = cell1 / this->multipliers[i - 1];
    cell1 = cell1 % this->multipliers[i - 1];
    if (position % 2 == 0) {
      if ((cell > this->multipliers[i - 1]) && (counter[i - 1] != 0)) {
        coboundary_elements.push_back(cell - this->multipliers[i - 1]);
      }
      if ((cell + this->multipliers[i - 1] < this->data.size()) && (counter[i - 1] != 2 * this->sizes[i - 1])) {
        coboundary_elements.push_back(cell + this->multipliers[i - 1]);
      }
    }
  }
  if (cell1 % 2 == 0) {
    if ((cell > 1) && (counter[0] != 0)) {
      coboundary_elements.push_back(cell - 1);
    }
    if ((cell + 1 < this->data.size()) && (counter[0] != 2 * this->sizes[0])) {
      coboundary_elements.push_back(cell + 1);
    }
  }
  return coboundary_elements;
}

template <typename T>
unsigned Bitmap_cubical_complex_base<T>::get_dimension_of_a_cell(std::size_t cell) const {
#ifdef DEBUG_TRACES
  std::clog << "\n\n\n Computing position o a cell of an index : " << cell << std::endl;
#endif
  unsigned dimension = 0;
  for (std::size_t i = this->multipliers.size(); i > 1; --i) {
    unsigned position = cell / this->multipliers[i - 1];
    std::size_t newcell = cell % this->multipliers[i - 1];

#ifdef DEBUG_TRACES
    std::clog << "i-1 :" << i - 1 << std::endl;
    std::clog << "cell : " << cell << std::endl;
    std::clog << "position : " << position << std::endl;
    std::clog << "multipliers[" << i - 1 << "] = " << this->multipliers[i - 1] << std::endl;
#endif

    if (position % 2 == 1) {
#ifdef DEBUG_TRACES
      std::clog << "Nonzero length in this direction \n";
#endif
      dimension++;
    }
    cell = newcell;
  }
  // Split out the last iteration to save a costly division by multipliers[0]=1
#ifdef DEBUG_TRACES
  std::clog << "i-1 :" << 0 << std::endl;
  std::clog << "cell : " << cell << std::endl;
  std::clog << "position : " << cell << std::endl;
  std::clog << "multipliers[" << 0 << "] = " << 1 << std::endl;
#endif

  if (cell % 2 == 1) {
#ifdef DEBUG_TRACES
    std::clog << "Nonzero length in this direction \n";
#endif
    dimension++;
  }

  return dimension;
}

template <typename T>
inline T& Bitmap_cubical_complex_base<T>::get_cell_data(std::size_t cell) {
  return this->data[cell];
}

template <typename T>
void Bitmap_cubical_complex_base<T>::impose_lower_star_filtration() {
  // this vector will be used to check which elements have already been taken care of in imposing lower star filtration
  std::vector<bool> is_this_cell_considered(this->data.size(), false);

  std::vector<std::size_t> indices_to_consider;
  // we assume here that we already have a filtration on the top dimensional cells and
  // we have to extend it to lower ones.
  for (auto it = this->top_dimensional_cells_iterator_begin();
       it != this->top_dimensional_cells_iterator_end(); ++it) {
    indices_to_consider.push_back(it.compute_index_in_bitmap());
  }

  while (indices_to_consider.size()) {
#ifdef DEBUG_TRACES
    std::clog << "indices_to_consider in this iteration \n";
    for (auto index : indices_to_consider) {
      std::clog << index << "  ";
    }
#endif
    std::vector<std::size_t> new_indices_to_consider;
    for (auto index : indices_to_consider) {
      std::vector<std::size_t> bd = this->get_boundary_of_a_cell(index);
      for (auto boundary : bd) {
#ifdef DEBUG_TRACES
        std::clog << "filtration of a cell : " << boundary << " is : " << this->data[boundary]
                  << " while of a cell: " << index << " is: " << this->data[index]
                  << std::endl;
#endif
        if (this->data[boundary] > this->data[index]) {
          this->data[boundary] = this->data[index];
#ifdef DEBUG_TRACES
          std::clog << "Setting the value of a cell : " << boundary
                    << " to : " << this->data[index] << std::endl;
#endif
        }
        if (is_this_cell_considered[boundary] == false) {
          new_indices_to_consider.push_back(boundary);
          is_this_cell_considered[boundary] = true;
        }
      }
    }
    indices_to_consider.swap(new_indices_to_consider);
  }
}

template <typename T>
template <class F>
void Bitmap_cubical_complex_base<T>::for_each_vertex_rec(F&&f, std::size_t base, int dim) {
  if (dim > 0)
    for(std::size_t i = 0; i < sizes[dim] + 1; ++i)
      for_each_vertex_rec(f, base + multipliers[dim] * 2 * i, dim - 1);
  else
    for(std::size_t i = 0; i < sizes[0] + 1; ++i)
      f(base + 2 * i);
}

template <typename T>
void Bitmap_cubical_complex_base<T>::impose_lower_star_filtration_from_vertices() {
  int max_dim = multipliers.size()-1;
  for (int dim = max_dim; dim >= 0; --dim)
    propagate_from_vertices_rec(dim, max_dim, 0);
}

template <typename T>
void Bitmap_cubical_complex_base<T>::propagate_from_vertices_rec (int special_dim, int current_dim, std::size_t base) {
  if (special_dim == current_dim) {
    propagate_from_vertices_rec(special_dim, current_dim - 1, base);
    return;
  }
  if (current_dim < 0) {
    std::size_t step = multipliers[special_dim];
    for(std::size_t i = 0; i < sizes[special_dim]; ++i) {
      std::size_t ref = base + step * 2 * i;
      data[ref + step] = std::max(data[ref], data[ref + 2 * step]);
    }
    return;
  }
  if (current_dim < special_dim)
    for(std::size_t i = 0; i < sizes[current_dim] + 1; ++i)
      propagate_from_vertices_rec(special_dim, current_dim - 1, base + multipliers[current_dim] * 2 * i);
  else
    for(std::size_t i = 0; i < 2 * sizes[current_dim] + 1; ++i)
      propagate_from_vertices_rec(special_dim, current_dim - 1, base + multipliers[current_dim] * i);
}

template <typename T>
bool compareFirstElementsOfTuples(const std::pair<std::pair<T, std::size_t>, char>& first,
                                  const std::pair<std::pair<T, std::size_t>, char>& second) {
  if (first.first.first < second.first.first) {
    return true;
  } else {
    if (first.first.first > second.first.first) {
      return false;
    }
    // in this case first.first.first == second.first.first, so we need to compare dimensions
    return first.second < second.second;
  }
}

}  // namespace cubical_complex

namespace Cubical_complex = cubical_complex;

}  // namespace Gudhi

#endif  // BITMAP_CUBICAL_COMPLEX_BASE_H_
