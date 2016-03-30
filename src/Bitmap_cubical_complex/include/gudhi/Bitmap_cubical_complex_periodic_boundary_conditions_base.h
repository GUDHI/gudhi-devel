/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2015  INRIA Sophia-Saclay (France)
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

#ifndef BITMAP_CUBICAL_COMPLEX_PERIODIC_BOUNDARY_CONDITIONS_BASE_H_
#define BITMAP_CUBICAL_COMPLEX_PERIODIC_BOUNDARY_CONDITIONS_BASE_H_

#include <gudhi/Bitmap_cubical_complex_base.h>

#include <cmath>
#include <limits>  // for numeric_limits<>
#include <vector>

namespace Gudhi {

namespace Cubical_complex {

// in this class, we are storing all the elements which are in normal bitmap (i.e. the bitmap without the periodic
// boundary conditions). But, we set up the iterators and the procedures to compute boundary and coboundary in the way
// that it is all right. We assume here that all the cells that are on the left / bottom and so on remains, while all
// the cells on the right / top are not in the Bitmap_cubical_complex_periodic_boundary_conditions_base

/**
 * @class Bitmap_cubical_complex_periodic_boundary_conditions_base
 * @brief Cubical complex with periodic boundary conditions represented as a bitmap.
 * @ingroup cubical_complex
 */
/**
 * This is a class implementing a bitmap data structure with periodic boundary conditions. Most of the functions are
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
  Bitmap_cubical_complex_periodic_boundary_conditions_base() { }
  /**
   * A constructor of Bitmap_cubical_complex_periodic_boundary_conditions_base class that takes the following
   * parameters: (1) vector with numbers of top dimensional cells in all dimensions and (2) vector of booleans. If
   * at i-th position of this vector there is true value, that means that periodic boundary conditions are to be
   * imposed in this direction. In case of false, the periodic boundary conditions will not be imposed in the direction
   * i.
   */
  Bitmap_cubical_complex_periodic_boundary_conditions_base(const std::vector<unsigned>& sizes,
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
  Bitmap_cubical_complex_periodic_boundary_conditions_base(const std::vector<unsigned>& dimensions,
                                                           const std::vector<T>& topDimensionalCells,
                                                           const std::vector< bool >& directions_in_which_periodic_b_cond_are_to_be_imposed);

  /**
   * Destructor of the Bitmap_cubical_complex_periodic_boundary_conditions_base class.
   **/
  virtual ~Bitmap_cubical_complex_periodic_boundary_conditions_base() {}

  // overwritten methods co compute boundary and coboundary
  /**
   * A version of a function that return boundary of a given cell for an object of
   * Bitmap_cubical_complex_periodic_boundary_conditions_base class.
   */
  virtual std::vector< size_t > get_boundary_of_a_cell(size_t cell) const;

  /**
   * A version of a function that return coboundary of a given cell for an object of
   * Bitmap_cubical_complex_periodic_boundary_conditions_base class.
   */
  virtual std::vector< size_t > get_coboundary_of_a_cell(size_t cell) const;

 protected:
  std::vector< bool > directions_in_which_periodic_b_cond_are_to_be_imposed;

  void set_up_containers(const std::vector<unsigned>& sizes) {
    unsigned multiplier = 1;
    for (size_t i = 0; i != sizes.size(); ++i) {
      this->sizes.push_back(sizes[i]);
      this->multipliers.push_back(multiplier);

      if (directions_in_which_periodic_b_cond_are_to_be_imposed[i]) {
        multiplier *= 2 * sizes[i];
      } else {
        multiplier *= 2 * sizes[i] + 1;
      }
    }
    // std::reverse( this->sizes.begin() , this->sizes.end() );
    this->data = std::vector<T>(multiplier, std::numeric_limits<T>::max());
    this->total_number_of_cells = multiplier;
  }
  Bitmap_cubical_complex_periodic_boundary_conditions_base(const std::vector<unsigned>& sizes);
  Bitmap_cubical_complex_periodic_boundary_conditions_base(const std::vector<unsigned>& dimensions,
                                                           const std::vector<T>& topDimensionalCells);
  void construct_complex_based_on_top_dimensional_cells(const std::vector<unsigned>& dimensions,
                                                        const std::vector<T>& topDimensionalCells,
                                                        const std::vector<bool>& directions_in_which_periodic_b_cond_are_to_be_imposed);
};

template <typename T>
void Bitmap_cubical_complex_periodic_boundary_conditions_base<T>::construct_complex_based_on_top_dimensional_cells(const std::vector<unsigned>& dimensions,
                                                                                                                   const std::vector<T>& topDimensionalCells,
                                                                                                                   const std::vector<bool>& directions_in_which_periodic_b_cond_are_to_be_imposed) {
  this->directions_in_which_periodic_b_cond_are_to_be_imposed = directions_in_which_periodic_b_cond_are_to_be_imposed;
  this->set_up_containers(dimensions);

  size_t i = 0;
  for (auto it = this->top_dimensional_cells_iterator_begin(); it != this->top_dimensional_cells_iterator_end(); ++it) {
    this->get_cell_data(*it) = topDimensionalCells[i];
    ++i;
  }
  this->impose_lower_star_filtration();
}

template <typename T>
Bitmap_cubical_complex_periodic_boundary_conditions_base<T>::Bitmap_cubical_complex_periodic_boundary_conditions_base(const std::vector<unsigned>& sizes,
                                                                                                                      const std::vector<bool>& directions_in_which_periodic_b_cond_are_to_be_imposed) {
  this->directions_in_which_periodic_b_cond_are_to_be_imposed(directions_in_which_periodic_b_cond_are_to_be_imposed);
  this->set_up_containers(sizes);
}

template <typename T>
Bitmap_cubical_complex_periodic_boundary_conditions_base<T>::Bitmap_cubical_complex_periodic_boundary_conditions_base(const char* perseus_style_file) {
  // for Perseus style files:
  bool dbg = false;

  std::ifstream inFiltration;
  inFiltration.open(perseus_style_file);
  unsigned dimensionOfData;
  inFiltration >> dimensionOfData;

  this->directions_in_which_periodic_b_cond_are_to_be_imposed = std::vector<bool>(dimensionOfData, false);

  std::vector<unsigned> sizes;
  sizes.reserve(dimensionOfData);
  for (size_t i = 0; i != dimensionOfData; ++i) {
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
    if (inFiltration.eof())break;

    if (dbg) {
      std::cerr << "Cell of an index : "
          << it.compute_index_in_bitmap()
          << " and dimension: "
          << this->get_dimension_of_a_cell(it.compute_index_in_bitmap())
          << " get the value : " << filtrationLevel << std::endl;
    }
    this->get_cell_data(*it) = filtrationLevel;
    ++it;
  }
  inFiltration.close();
  this->impose_lower_star_filtration();
}

template <typename T>
Bitmap_cubical_complex_periodic_boundary_conditions_base<T>::Bitmap_cubical_complex_periodic_boundary_conditions_base(const std::vector<unsigned>& sizes) {
  this->directions_in_which_periodic_b_cond_are_to_be_imposed = std::vector<bool>(sizes.size(), false);
  this->set_up_containers(sizes);
}

template <typename T>
Bitmap_cubical_complex_periodic_boundary_conditions_base<T>::Bitmap_cubical_complex_periodic_boundary_conditions_base(const std::vector<unsigned>& dimensions,
                                                                                                                      const std::vector<T>& topDimensionalCells) {
  std::vector<bool> directions_in_which_periodic_b_cond_are_to_be_imposed = std::vector<bool>(dimensions.size(), false);
  this->construct_complex_based_on_top_dimensional_cells(dimensions, topDimensionalCells,
                                                         directions_in_which_periodic_b_cond_are_to_be_imposed);
}

template <typename T>
Bitmap_cubical_complex_periodic_boundary_conditions_base<T>::
Bitmap_cubical_complex_periodic_boundary_conditions_base(const std::vector<unsigned>& dimensions,
                                                         const std::vector<T>& topDimensionalCells,
                                                         const std::vector<bool>& directions_in_which_periodic_b_cond_are_to_be_imposed) {
  this->construct_complex_based_on_top_dimensional_cells(dimensions, topDimensionalCells,
                                                         directions_in_which_periodic_b_cond_are_to_be_imposed);
}

// ***********************Methods************************ //

template <typename T>
std::vector< size_t > Bitmap_cubical_complex_periodic_boundary_conditions_base<T>::get_boundary_of_a_cell(size_t cell) const {
  bool dbg = false;
  if (dbg) {
    std::cerr << "Computations of boundary of a cell : " << cell << std::endl;
  }

  std::vector< size_t > boundary_elements;
  size_t cell1 = cell;
  for (size_t i = this->multipliers.size(); i != 0; --i) {
    unsigned position = cell1 / this->multipliers[i - 1];
    // this cell have a nonzero length in this direction, therefore we can compute its boundary in this direction.

    if (position % 2 == 1) {
      // if there are no periodic boundary conditions in this direction, we do not have to do anything.
      if (!directions_in_which_periodic_b_cond_are_to_be_imposed[i - 1]) {
        // std::cerr << "A\n";
        boundary_elements.push_back(cell - this->multipliers[ i - 1 ]);
        boundary_elements.push_back(cell + this->multipliers[ i - 1 ]);
        if (dbg) {
          std::cerr << cell - this->multipliers[ i - 1 ] << " " << cell + this->multipliers[ i - 1 ] << " ";
        }
      } else {
        // in this direction we have to do boundary conditions. Therefore, we need to check if we are not at the end.
        if (position != 2 * this->sizes[ i - 1 ] - 1) {
          // std::cerr << "B\n";
          boundary_elements.push_back(cell - this->multipliers[ i - 1 ]);
          boundary_elements.push_back(cell + this->multipliers[ i - 1 ]);
          if (dbg) {
            std::cerr << cell - this->multipliers[ i - 1 ] << " " << cell + this->multipliers[ i - 1 ] << " ";
          }
        } else {
          // std::cerr << "C\n";
          boundary_elements.push_back(cell - this->multipliers[ i - 1 ]);
          boundary_elements.push_back(cell - (2 * this->sizes[ i - 1 ] - 1) * this->multipliers[ i - 1 ]);
          if (dbg) {
            std::cerr << cell - this->multipliers[ i - 1 ] << " " <<
                cell - (2 * this->sizes[ i - 1 ] - 1) * this->multipliers[ i - 1 ] << " ";
          }
        }
      }
    }
    cell1 = cell1 % this->multipliers[i - 1];
  }
  return boundary_elements;
}

template <typename T>
std::vector< size_t > Bitmap_cubical_complex_periodic_boundary_conditions_base<T>::get_coboundary_of_a_cell(size_t cell) const {
  std::vector<unsigned> counter = this->compute_counter_for_given_cell(cell);
  std::vector< size_t > coboundary_elements;
  size_t cell1 = cell;
  for (size_t i = this->multipliers.size(); i != 0; --i) {
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
          coboundary_elements.push_back(cell + (2 * this->sizes[ i - 1 ] - 1) * this->multipliers[i - 1]);
        }
      }
    }

    cell1 = cell1 % this->multipliers[i - 1];
  }
  return coboundary_elements;
}

}  // namespace Cubical_complex

}  // namespace Gudhi

#endif  // BITMAP_CUBICAL_COMPLEX_PERIODIC_BOUNDARY_CONDITIONS_BASE_H_
