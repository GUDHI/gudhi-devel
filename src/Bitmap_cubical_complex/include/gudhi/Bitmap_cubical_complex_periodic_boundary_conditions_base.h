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
  Bitmap_cubical_complex_periodic_boundary_conditions_base(const char* perseusStyleFile);
  /**
   * A constructor of Bitmap_cubical_complex_periodic_boundary_conditions_base class that takes the following
   * parameters: (1) vector with numbers of top dimensional cells in all dimensions and (2) vector of top dimensional
   * cells (ordered lexicographically) and (3) vector of booleans. If at i-th position of this vector there is true
   * value, that means that periodic boundary conditions are to be imposed in this direction. In case of false, the
   * periodic boundary conditions will not be imposed in the direction i.
   */
  Bitmap_cubical_complex_periodic_boundary_conditions_base(
      const std::vector<unsigned>& dimensions, const std::vector<T>& topDimensionalCells,
      const std::vector<bool>& directions_in_which_periodic_b_cond_are_to_be_imposed);

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

 protected:
  std::vector<bool> directions_in_which_periodic_b_cond_are_to_be_imposed;

  void set_up_containers(const std::vector<unsigned>& sizes) {
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
    this->data = std::vector<T>(multiplier, std::numeric_limits<T>::infinity());
    this->total_number_of_cells = multiplier;
  }
  Bitmap_cubical_complex_periodic_boundary_conditions_base(const std::vector<unsigned>& sizes);
  Bitmap_cubical_complex_periodic_boundary_conditions_base(const std::vector<unsigned>& dimensions,
                                                           const std::vector<T>& topDimensionalCells);

  /**
   * A procedure used to construct the data structures in the class.
  **/
  void construct_complex_based_on_top_dimensional_cells(
      const std::vector<unsigned>& dimensions, const std::vector<T>& topDimensionalCells,
      const std::vector<bool>& directions_in_which_periodic_b_cond_are_to_be_imposed);
};

template <typename T>
void Bitmap_cubical_complex_periodic_boundary_conditions_base<T>::construct_complex_based_on_top_dimensional_cells(
    const std::vector<unsigned>& dimensions, const std::vector<T>& topDimensionalCells,
    const std::vector<bool>& directions_in_which_periodic_b_cond_are_to_be_imposed) {
  this->directions_in_which_periodic_b_cond_are_to_be_imposed = directions_in_which_periodic_b_cond_are_to_be_imposed;
  this->set_up_containers(dimensions);

  std::size_t i = 0;
  for (auto it = this->top_dimensional_cells_iterator_begin(); it != this->top_dimensional_cells_iterator_end(); ++it) {
    this->get_cell_data(*it) = topDimensionalCells[i];
    ++i;
  }
  this->impose_lower_star_filtration();
}

template <typename T>
Bitmap_cubical_complex_periodic_boundary_conditions_base<T>::Bitmap_cubical_complex_periodic_boundary_conditions_base(
    const std::vector<unsigned>& sizes,
    const std::vector<bool>& directions_in_which_periodic_b_cond_are_to_be_imposed) {
  this->directions_in_which_periodic_b_cond_are_to_be_imposed(directions_in_which_periodic_b_cond_are_to_be_imposed);
  this->set_up_containers(sizes);
}

template <typename T>
Bitmap_cubical_complex_periodic_boundary_conditions_base<T>::Bitmap_cubical_complex_periodic_boundary_conditions_base(
    const char* perseus_style_file) {
  // for Perseus style files:
  bool dbg = false;

  std::ifstream inFiltration;
  inFiltration.open(perseus_style_file);
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
  this->set_up_containers(sizes);

  typename Bitmap_cubical_complex_periodic_boundary_conditions_base<T>::Top_dimensional_cells_iterator it(*this);
  it = this->top_dimensional_cells_iterator_begin();

  while (!inFiltration.eof()) {
    double filtrationLevel;
    inFiltration >> filtrationLevel;
    if (inFiltration.eof()) break;

    if (dbg) {
      std::clog << "Cell of an index : " << it.compute_index_in_bitmap()
                << " and dimension: " << this->get_dimension_of_a_cell(it.compute_index_in_bitmap())
                << " get the value : " << filtrationLevel << std::endl;
    }
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
  this->set_up_containers(sizes);
}

template <typename T>
Bitmap_cubical_complex_periodic_boundary_conditions_base<T>::Bitmap_cubical_complex_periodic_boundary_conditions_base(
    const std::vector<unsigned>& dimensions, const std::vector<T>& topDimensionalCells) {
  std::vector<bool> directions_in_which_periodic_b_cond_are_to_be_imposed = std::vector<bool>(dimensions.size(), false);
  this->construct_complex_based_on_top_dimensional_cells(dimensions, topDimensionalCells,
                                                         directions_in_which_periodic_b_cond_are_to_be_imposed);
}

template <typename T>
Bitmap_cubical_complex_periodic_boundary_conditions_base<T>::Bitmap_cubical_complex_periodic_boundary_conditions_base(
    const std::vector<unsigned>& dimensions, const std::vector<T>& topDimensionalCells,
    const std::vector<bool>& directions_in_which_periodic_b_cond_are_to_be_imposed) {
  this->construct_complex_based_on_top_dimensional_cells(dimensions, topDimensionalCells,
                                                         directions_in_which_periodic_b_cond_are_to_be_imposed);
}

// ***********************Methods************************ //

template <typename T>
std::vector<std::size_t> Bitmap_cubical_complex_periodic_boundary_conditions_base<T>::get_boundary_of_a_cell(
    std::size_t cell) const {
  bool dbg = false;
  if (dbg) {
    std::clog << "Computations of boundary of a cell : " << cell << std::endl;
  }

  std::vector<std::size_t> boundary_elements;
  boundary_elements.reserve(this->dimension() * 2);
  std::size_t cell1 = cell;
  std::size_t sum_of_dimensions = 0;

  for (std::size_t i = this->multipliers.size(); i != 0; --i) {
    unsigned position = cell1 / this->multipliers[i - 1];
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
        if (dbg) {
          std::clog << cell - this->multipliers[i - 1] << " " << cell + this->multipliers[i - 1] << " ";
        }
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
          if (dbg) {
            std::clog << cell - this->multipliers[i - 1] << " " << cell + this->multipliers[i - 1] << " ";
          }
        } else {
          if (sum_of_dimensions % 2) {
            boundary_elements.push_back(cell - this->multipliers[i - 1]);
            boundary_elements.push_back(cell - (2 * this->sizes[i - 1] - 1) * this->multipliers[i - 1]);
          } else {
            boundary_elements.push_back(cell - (2 * this->sizes[i - 1] - 1) * this->multipliers[i - 1]);
            boundary_elements.push_back(cell - this->multipliers[i - 1]);
          }
          if (dbg) {
            std::clog << cell - this->multipliers[i - 1] << " "
                      << cell - (2 * this->sizes[i - 1] - 1) * this->multipliers[i - 1] << " ";
          }
        }
      }
      ++sum_of_dimensions;
    }
    cell1 = cell1 % this->multipliers[i - 1];
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

    cell1 = cell1 % this->multipliers[i - 1];
  }
  return coboundary_elements;
}

}  // namespace cubical_complex

namespace Cubical_complex = cubical_complex;

}  // namespace Gudhi

#endif  // BITMAP_CUBICAL_COMPLEX_PERIODIC_BOUNDARY_CONDITIONS_BASE_H_
