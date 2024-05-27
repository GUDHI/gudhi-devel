/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2015 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef BITMAP_CUBICAL_COMPLEX_PERIODIC_BOUNDARY_CONDITIONS_BASE_H_
#define BITMAP_CUBICAL_COMPLEX_PERIODIC_BOUNDARY_CONDITIONS_BASE_H_

#include <gudhi/Bitmap_cubical_complex_base.h>
#include <gudhi/Debug_utils.h>

#include <cmath>
#include <limits>  // for numeric_limits<>
#include <vector>
#include <stdexcept>
#include <cstddef>

namespace Gudhi {

namespace cubical_complex {

// in this class, we are storing all the elements which are in normal bitmap (i.e. the bitmap without the periodic
// boundary conditions). But, we set up the iterators and the procedures to compute boundary and coboundary in the way
// that it is all right. We assume here that all the cells that are on the left / bottom and so on remains, while all
// the cells on the right / top are not in the Bitmap_cubical_complex_periodic_boundary_conditions_base

/**
 * @brief Cubical complex with periodic boundary conditions represented as a bitmap.
 * @ingroup cubical_complex
 * @details This is a class implementing a bitmap data structure with periodic boundary conditions. Most of the
 * functions are
 * identical to the functions from Bitmap_cubical_complex_base.
 * The ones that needed to be updated are the constructors and get_boundary_of_a_cell and get_coboundary_of_a_cell.
 */
template <typename T>
class Bitmap_cubical_complex_periodic_boundary_conditions_base : public Bitmap_cubical_complex_base<T> {
 public:
  // constructors that take an extra parameter:

  /**
   * Default constructor of Bitmap_cubical_complex_periodic_boundary_conditions_base class.
   */
  Bitmap_cubical_complex_periodic_boundary_conditions_base() {}
  /**
   * A constructor of Bitmap_cubical_complex_periodic_boundary_conditions_base class that takes the following
   * parameters: (1) vector with numbers of top dimensional cells in all dimensions and (2) vector of booleans. If
   * at i-th position of this vector there is true value, that means that periodic boundary conditions are to be
   * imposed in this direction. In case of false, the periodic boundary conditions will not be imposed in the direction
   * i.
   */
  Bitmap_cubical_complex_periodic_boundary_conditions_base(
      const std::vector<unsigned>& sizes,
      const std::vector<bool>& directions_in_which_periodic_b_cond_are_to_be_imposed);
  /**
   * A constructor of Bitmap_cubical_complex_periodic_boundary_conditions_base class that takes the name of Perseus
   * style file as an input. Please consult the documentation about the specification of the file.
   */
  explicit Bitmap_cubical_complex_periodic_boundary_conditions_base(const char* perseusStyleFile);
  /**
   * A constructor of Bitmap_cubical_complex_periodic_boundary_conditions_base class that takes the following
   * parameters: (1) vector with numbers of top dimensional cells in all dimensions and (2) vector of top dimensional
   * cells (ordered lexicographically) and (3) vector of booleans. If at i-th position of this vector there is true
   * value, that means that periodic boundary conditions are to be imposed in this direction. In case of false, the
   * periodic boundary conditions will not be imposed in the direction i.
   */
  Bitmap_cubical_complex_periodic_boundary_conditions_base(
      const std::vector<unsigned>& dimensions, const std::vector<T>& cells,
      const std::vector<bool>& directions_in_which_periodic_b_cond_are_to_be_imposed,
      bool input_top_cells);

  /**
   * Destructor of the Bitmap_cubical_complex_periodic_boundary_conditions_base class.
   **/
  virtual ~Bitmap_cubical_complex_periodic_boundary_conditions_base() {}

  // overwritten methods co compute boundary and coboundary
  /**
   * A version of a function that return boundary of a given cell for an object of
   * Bitmap_cubical_complex_periodic_boundary_conditions_base class.
   * The boundary elements are guaranteed to be returned so that the
   * incidence coefficients are alternating.
   */
  virtual std::vector<std::size_t> get_boundary_of_a_cell(std::size_t cell) const override;

  /**
   * A version of a function that return coboundary of a given cell for an object of
   * Bitmap_cubical_complex_periodic_boundary_conditions_base class.
   * Note that unlike in the case of boundary, over here the elements are
   * not guaranteed to be returned with alternating incidence numbers.
   * To compute incidence between cells use compute_incidence_between_cells
   * procedure
   */
  virtual std::vector<std::size_t> get_coboundary_of_a_cell(std::size_t cell) const override;

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
  *[b_{n},e_{n}]s \f$
  * where \f$ a = b_{j}\f$ or \f$ a = e_{j}\f$, the incidence between \f$A\f$ and \f$B\f$
  * computed by this procedure is given by formula:
  * \f$ c\ (-1)^{\sum_{i=1}^{j-1} dim [b_{i},e_{i}]}  \f$
  * Where \f$ dim [b_{i},e_{i}] = 0 \f$ if \f$ b_{i}=e_{i} \f$ and 1 in other case.
  * c is -1 if \f$ a = b_{j}\f$ and 1 if \f$ a = e_{j}\f$.
  * @exception std::logic_error In case when the cube \f$B\f$ is not n-1
  * dimensional face of a cube \f$A\f$.
  **/
  virtual int compute_incidence_between_cells(std::size_t coface, std::size_t face) const override {
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
    if ((coface_counter[number_of_position_in_which_counters_do_not_agree] + 1 ==
         face_counter[number_of_position_in_which_counters_do_not_agree]) ||
        ((coface_counter[number_of_position_in_which_counters_do_not_agree] != 1) &&
         (face_counter[number_of_position_in_which_counters_do_not_agree] == 0))) {
      incidence *= -1;
    }

    return incidence;
  }

  // The non-periodic code works for top dimensional cells, but not vertices.
  class Vertices_iterator {
   public:
    typedef std::input_iterator_tag iterator_category;
    typedef std::size_t value_type;
    typedef std::ptrdiff_t difference_type;
    typedef value_type* pointer;
    typedef value_type reference;

    Vertices_iterator(Bitmap_cubical_complex_periodic_boundary_conditions_base* b) : counter(b->dimension()), b(b) {}

    Vertices_iterator operator++() {
      // first find first element of the counter that can be increased:
      std::size_t dim = 0;
      while ((dim != this->b->dimension()) && (this->counter[dim] == this->b->sizes[dim] - this->b->directions_in_which_periodic_b_cond_are_to_be_imposed[dim])) ++dim;

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
    friend class Bitmap_cubical_complex_periodic_boundary_conditions_base;

   protected:
    std::vector<std::size_t> counter;
    Bitmap_cubical_complex_periodic_boundary_conditions_base* b;
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
      a.counter[i] = this->sizes[i] - this->directions_in_which_periodic_b_cond_are_to_be_imposed[i];
    }
    a.counter[0]++;
    return a;
  }

  /*
   * @brief Vertices_iterator_range class provides ranges for Vertices_iterator_range
   */
  class Vertices_range {
   public:
    Vertices_range(Bitmap_cubical_complex_periodic_boundary_conditions_base* b) : b(b) {}

    Vertices_iterator begin() { return b->vertices_iterator_begin(); }

    Vertices_iterator end() { return b->vertices_iterator_end(); }

   private:
    Bitmap_cubical_complex_periodic_boundary_conditions_base* b;
  };

  /* Returns a range over all vertices. */
  Vertices_range vertices_range() { return Vertices_range(this); }

  /**
   * Set cells filtrations given those of the vertices, and based on lower star filtration.
   * This is already called by the relevant constructors.
   **/
  void impose_lower_star_filtration_from_vertices();

 protected:
  std::vector<bool> directions_in_which_periodic_b_cond_are_to_be_imposed;

  void set_up_containers(const std::vector<unsigned>& sizes, bool is_pos_inf) {
    // The fact that multipliers[0]=1 is relied on by optimizations in other functions
    unsigned multiplier = 1;
    for (std::size_t i = 0; i != sizes.size(); ++i) {
      this->sizes.push_back(sizes[i]);
      this->multipliers.push_back(multiplier);

      if (directions_in_which_periodic_b_cond_are_to_be_imposed[i]) {
        multiplier *= 2 * sizes[i];
      } else {
        multiplier *= 2 * sizes[i] + 1;
      }
    }
    // std::reverse( this->sizes.begin() , this->sizes.end() );
    if (is_pos_inf)
      this->data = std::vector<T>(multiplier, std::numeric_limits<T>::infinity());
    else
      this->data = std::vector<T>(multiplier, -std::numeric_limits<T>::infinity());
  }
  Bitmap_cubical_complex_periodic_boundary_conditions_base(const std::vector<unsigned>& sizes);
  Bitmap_cubical_complex_periodic_boundary_conditions_base(const std::vector<unsigned>& dimensions,
                                                           const std::vector<T>& cells,
                                                           bool input_top_cells);

  /**
   * Procedures used to construct the data structures in the class.
  **/
  void construct_complex_based_on_top_dimensional_cells(
      const std::vector<unsigned>& dimensions, const std::vector<T>& topDimensionalCells,
      const std::vector<bool>& directions_in_which_periodic_b_cond_are_to_be_imposed);
  void construct_complex_based_on_vertices(const std::vector<unsigned>& dimensions, const std::vector<T>& vertices,
      const std::vector<bool>& directions_in_which_periodic_b_cond_are_to_be_imposed);
};

template <typename T>
void Bitmap_cubical_complex_periodic_boundary_conditions_base<T>::construct_complex_based_on_top_dimensional_cells(
    const std::vector<unsigned>& dimensions, const std::vector<T>& topDimensionalCells,
    const std::vector<bool>& directions_in_which_periodic_b_cond_are_to_be_imposed) {
  this->directions_in_which_periodic_b_cond_are_to_be_imposed = directions_in_which_periodic_b_cond_are_to_be_imposed;
  this->set_up_containers(dimensions, true);

  std::size_t i = 0;
  for (auto it = this->top_dimensional_cells_iterator_begin(); it != this->top_dimensional_cells_iterator_end(); ++it) {
    this->get_cell_data(*it) = topDimensionalCells[i];
    ++i;
  }
  this->impose_lower_star_filtration();
}

template <typename T>
void Bitmap_cubical_complex_periodic_boundary_conditions_base<T>::construct_complex_based_on_vertices(
    const std::vector<unsigned>& dimensions, const std::vector<T>& vertices,
    const std::vector<bool>& directions_in_which_periodic_b_cond_are_to_be_imposed) {
  this->directions_in_which_periodic_b_cond_are_to_be_imposed = directions_in_which_periodic_b_cond_are_to_be_imposed;

  std::vector<unsigned> top_cells_sizes;
  std::transform (dimensions.begin(), dimensions.end(), directions_in_which_periodic_b_cond_are_to_be_imposed.begin(),
      std::back_inserter(top_cells_sizes), [](unsigned i, bool b){ return i - !b;});
  this->set_up_containers(top_cells_sizes, false);

  std::size_t i = 0;
  for (auto it = this->vertices_iterator_begin(); it != this->vertices_iterator_end(); ++it) {
    this->get_cell_data(*it) = vertices[i];
    ++i;
  }
  this->impose_lower_star_filtration_from_vertices();
}

template <typename T>
Bitmap_cubical_complex_periodic_boundary_conditions_base<T>::Bitmap_cubical_complex_periodic_boundary_conditions_base(
    const std::vector<unsigned>& sizes,
    const std::vector<bool>& directions_in_which_periodic_b_cond_are_to_be_imposed) 
    : directions_in_which_periodic_b_cond_are_to_be_imposed(directions_in_which_periodic_b_cond_are_to_be_imposed) {
  this->set_up_containers(sizes, true);
}

template <typename T>
Bitmap_cubical_complex_periodic_boundary_conditions_base<T>::Bitmap_cubical_complex_periodic_boundary_conditions_base(
    const char* perseus_style_file) {
  // for Perseus style files:

  std::ifstream inFiltration(perseus_style_file);
  if(!inFiltration) throw std::ios_base::failure(std::string("Could not open the file ") + perseus_style_file);
  unsigned dimensionOfData;
  inFiltration >> dimensionOfData;

  this->directions_in_which_periodic_b_cond_are_to_be_imposed = std::vector<bool>(dimensionOfData, false);

  std::vector<unsigned> sizes;
  sizes.reserve(dimensionOfData);
  for (std::size_t i = 0; i != dimensionOfData; ++i) {
    int size_in_this_dimension;
    inFiltration >> size_in_this_dimension;
    if (size_in_this_dimension < 0) {
      this->directions_in_which_periodic_b_cond_are_to_be_imposed[i] = true;
    }
    sizes.push_back(abs(size_in_this_dimension));
  }
  this->set_up_containers(sizes, true);

  typename Bitmap_cubical_complex_periodic_boundary_conditions_base<T>::Top_dimensional_cells_iterator
  it = this->top_dimensional_cells_iterator_begin();

  while (!inFiltration.eof()) {
    double filtrationLevel;
    inFiltration >> filtrationLevel;
    if (inFiltration.eof()) break;

#ifdef DEBUG_TRACES
    std::clog << "Cell of an index : " << it.compute_index_in_bitmap()
              << " and dimension: " << this->get_dimension_of_a_cell(it.compute_index_in_bitmap())
              << " get the value : " << filtrationLevel << std::endl;
#endif
    this->get_cell_data(*it) = filtrationLevel;
    ++it;
  }
  inFiltration.close();
  this->impose_lower_star_filtration();
}

template <typename T>
Bitmap_cubical_complex_periodic_boundary_conditions_base<T>::Bitmap_cubical_complex_periodic_boundary_conditions_base(
    const std::vector<unsigned>& sizes) {
  this->directions_in_which_periodic_b_cond_are_to_be_imposed = std::vector<bool>(sizes.size(), false);
  this->set_up_containers(sizes, true);
}

template <typename T>
Bitmap_cubical_complex_periodic_boundary_conditions_base<T>::Bitmap_cubical_complex_periodic_boundary_conditions_base(
    const std::vector<unsigned>& dimensions, const std::vector<T>& cells, bool input_top_cells) {
  std::vector<bool> directions_in_which_periodic_b_cond_are_to_be_imposed = std::vector<bool>(dimensions.size(), false);
  if (input_top_cells) {
    this->construct_complex_based_on_top_dimensional_cells(dimensions, cells,
                                                           directions_in_which_periodic_b_cond_are_to_be_imposed);
  } else {
    this->construct_complex_based_on_vertices(dimensions, cells,
                                              directions_in_which_periodic_b_cond_are_to_be_imposed);
  }
}

template <typename T>
Bitmap_cubical_complex_periodic_boundary_conditions_base<T>::Bitmap_cubical_complex_periodic_boundary_conditions_base(
    const std::vector<unsigned>& dimensions, const std::vector<T>& cells,
    const std::vector<bool>& directions_in_which_periodic_b_cond_are_to_be_imposed,
    bool input_top_cells) {
  if(input_top_cells) {
    this->construct_complex_based_on_top_dimensional_cells(dimensions, cells,
                                                           directions_in_which_periodic_b_cond_are_to_be_imposed);
  } else {
    this->construct_complex_based_on_vertices(dimensions, cells,
                                              directions_in_which_periodic_b_cond_are_to_be_imposed);
  }
}

// ***********************Methods************************ //

template <typename T>
std::vector<std::size_t> Bitmap_cubical_complex_periodic_boundary_conditions_base<T>::get_boundary_of_a_cell(
    std::size_t cell) const {
#ifdef DEBUG_TRACES
  std::clog << "Computations of boundary of a cell : " << cell << std::endl;
#endif

  std::vector<std::size_t> boundary_elements;
  boundary_elements.reserve(this->dimension() * 2);
  std::size_t cell1 = cell;
  std::size_t sum_of_dimensions = 0;

  for (std::size_t i = this->multipliers.size(); i != 0; --i) {
    unsigned position = cell1 / this->multipliers[i - 1];
    cell1 = cell1 % this->multipliers[i - 1];
    // this cell have a nonzero length in this direction, therefore we can compute its boundary in this direction.
    if (position % 2 == 1) {
      // if there are no periodic boundary conditions in this direction, we do not have to do anything.
      if (!directions_in_which_periodic_b_cond_are_to_be_imposed[i - 1]) {
        if (sum_of_dimensions % 2) {
          boundary_elements.push_back(cell - this->multipliers[i - 1]);
          boundary_elements.push_back(cell + this->multipliers[i - 1]);
        } else {
          boundary_elements.push_back(cell + this->multipliers[i - 1]);
          boundary_elements.push_back(cell - this->multipliers[i - 1]);
        }
#ifdef DEBUG_TRACES
        std::clog << cell - this->multipliers[i - 1] << " " << cell + this->multipliers[i - 1] << " ";
#endif
      } else {
        // in this direction we have to do boundary conditions. Therefore, we need to check if we are not at the end.
        if (position != 2 * this->sizes[i - 1] - 1) {
          if (sum_of_dimensions % 2) {
            boundary_elements.push_back(cell - this->multipliers[i - 1]);
            boundary_elements.push_back(cell + this->multipliers[i - 1]);
          } else {
            boundary_elements.push_back(cell + this->multipliers[i - 1]);
            boundary_elements.push_back(cell - this->multipliers[i - 1]);
          }
#ifdef DEBUG_TRACES
          std::clog << cell - this->multipliers[i - 1] << " " << cell + this->multipliers[i - 1] << " ";
#endif
        } else {
          if (sum_of_dimensions % 2) {
            boundary_elements.push_back(cell - this->multipliers[i - 1]);
            boundary_elements.push_back(cell - (2 * this->sizes[i - 1] - 1) * this->multipliers[i - 1]);
          } else {
            boundary_elements.push_back(cell - (2 * this->sizes[i - 1] - 1) * this->multipliers[i - 1]);
            boundary_elements.push_back(cell - this->multipliers[i - 1]);
          }
#ifdef DEBUG_TRACES
          std::clog << cell - this->multipliers[i - 1] << " "
                    << cell - (2 * this->sizes[i - 1] - 1) * this->multipliers[i - 1] << " ";
#endif
        }
      }
      ++sum_of_dimensions;
    }
  }
  return boundary_elements;
}

template <typename T>
std::vector<std::size_t> Bitmap_cubical_complex_periodic_boundary_conditions_base<T>::get_coboundary_of_a_cell(
    std::size_t cell) const {
  std::vector<unsigned> counter = this->compute_counter_for_given_cell(cell);
  std::vector<std::size_t> coboundary_elements;
  std::size_t cell1 = cell;
  for (std::size_t i = this->multipliers.size(); i != 0; --i) {
    unsigned position = cell1 / this->multipliers[i - 1];
    cell1 = cell1 % this->multipliers[i - 1];
    // if the cell has zero length in this direction, then it will have cbd in this direction.
    if (position % 2 == 0) {
      if (!this->directions_in_which_periodic_b_cond_are_to_be_imposed[i - 1]) {
        // no periodic boundary conditions in this direction
        if ((counter[i - 1] != 0) && (cell > this->multipliers[i - 1])) {
          coboundary_elements.push_back(cell - this->multipliers[i - 1]);
        }
        if ((counter[i - 1] != 2 * this->sizes[i - 1]) && (cell + this->multipliers[i - 1] < this->data.size())) {
          coboundary_elements.push_back(cell + this->multipliers[i - 1]);
        }
      } else {
        // we want to have periodic boundary conditions in this direction
        if (counter[i - 1] != 0) {
          coboundary_elements.push_back(cell - this->multipliers[i - 1]);
          coboundary_elements.push_back(cell + this->multipliers[i - 1]);
        } else {
          // in this case counter[i-1] == 0.
          coboundary_elements.push_back(cell + this->multipliers[i - 1]);
          coboundary_elements.push_back(cell + (2 * this->sizes[i - 1] - 1) * this->multipliers[i - 1]);
        }
      }
    }
  }
  return coboundary_elements;
}

// Exact same code as the non-periodic version, duplicated for now to ensure it calls the right version of vertices_iterator_begin.
template <typename T>
void Bitmap_cubical_complex_periodic_boundary_conditions_base<T>::impose_lower_star_filtration_from_vertices() {
  // this vector will be used to check which elements have already been taken care of in imposing lower star filtration
  std::vector<bool> is_this_cell_considered(this->data.size(), false);

  std::vector<std::size_t> indices_to_consider;
  // we assume here that we already have a filtration on the vertices and
  // we have to extend it to higher ones.
  for (auto it = this->vertices_iterator_begin();
       it != this->vertices_iterator_end(); ++it) {
    indices_to_consider.push_back(it.compute_index_in_bitmap());
  }

  while (indices_to_consider.size()) { // Iteration on the dimension
#ifdef DEBUG_TRACES
    std::clog << "indices_to_consider in this iteration \n";
    for (auto index : indices_to_consider) {
      std::clog << index << "  ";
    }
#endif
    std::vector<std::size_t> new_indices_to_consider;
    for (auto index : indices_to_consider) {
      std::vector<std::size_t> cbd = this->get_coboundary_of_a_cell(index);
      for (auto coboundary : cbd) {
#ifdef DEBUG_TRACES
        std::clog << "filtration of a cell : " << coboundary << " is : " << this->data[coboundary]
                  << " while of a cell: " << index << " is: " << this->data[index]
                  << std::endl;
#endif
        if (this->data[coboundary] < this->data[index]) {
          this->data[coboundary] = this->data[index];
#ifdef DEBUG_TRACES
          std::clog << "Setting the value of a cell : " << coboundary
                    << " to : " << this->data[index] << std::endl;
#endif
        }
        if (is_this_cell_considered[coboundary] == false) {
          new_indices_to_consider.push_back(coboundary);
          is_this_cell_considered[coboundary] = true;
        }
      }
    }
    indices_to_consider.swap(new_indices_to_consider);
  }
}
}  // namespace cubical_complex

namespace Cubical_complex = cubical_complex;

}  // namespace Gudhi

#endif  // BITMAP_CUBICAL_COMPLEX_PERIODIC_BOUNDARY_CONDITIONS_BASE_H_
