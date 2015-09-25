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

#ifndef BITMAP_CUBICAL_COMPLEX_BASE_H_
#define BITMAP_CUBICAL_COMPLEX_BASE_H_

#include <iostream>
#include <string>
#include <vector>
#include <string>
#include <cstdlib>
#include <climits>
#include <fstream>
#include <algorithm>
#include <iterator>

#include <gudhi/counter.h>

using namespace std;

/**
 * This is a class implementing a basic bitmap data structure to store cubical complexes. It implements only the most basic subroutines.
 * The idea of the bitmap is the following. Our aim is to have a memory efficient data structure to store d-dimensional cubical complex C being a cubical decomposition
 * of a rectangular region of a space. This is achieved by storing C as a vector of bits (this is where the name 'bitmap' came from). Each cell is represented by a single
 * bit (in case of black and white bitmaps, or by a single element of a type T (here T is a filtration type of a bitmap, typically a double). All the informations needed for homology and
 * persistent homology computations (like dimension of a cell, boundary and coboundary elements of a cell, are then obtained from the position of the element in C.
 */
template <typename T>
class Bitmap_cubical_complex_base {
 public:
  /**
   * There are a few constructors of a Bitmap_cubical_complex_base class. First one, that takes vector<unsigned>, creates an empty bitmap of a dimension equal to the number of elements in the
   * input vector and size in the i-th dimension equal the number in the position i-of the input vector.
   */
  Bitmap_cubical_complex_base(std::vector<unsigned> sizes_);
  /**
   * The second constructor takes as a input a Perseus style file. For more details, please consult the documentations of Perseus software as well as examples attached to this
   * implementation.
   **/
  Bitmap_cubical_complex_base(char* perseusStyleFile_);
  /**
   * The last constructor of a Bitmap_cubical_complex_base class accepts vector of dimensions (as the first one) together with vector of filtration values of top dimensional cells.
   **/
  Bitmap_cubical_complex_base(std::vector<unsigned> dimensions_, std::vector<T> topDimensionalCells_);

  /**
   * The functions get_boundary_of_a_cell, get_coboundary_of_a_cell and get_cell_data are the basic functions that compute boundary / coboundary / dimension and the filtration
   * value form a position of a cell in the structure of a bitmap. The input parameter of all of those function is a non-negative integer, indicating a position of a cube in the data structure.
   * In the case of functions that compute (co)boundary, the output is a vector if non-negative integers pointing to the positions of (co)boundary element of the input cell.
   */
  inline std::vector< size_t > get_boundary_of_a_cell(size_t cell_);
  /**
   * The functions get_boundary_of_a_cell, get_coboundary_of_a_cell, get_dimension_of_a_cell and get_cell_data are the basic functions that compute boundary / coboundary / dimension and the filtration
   * value form a position of a cell in the structure of a bitmap. The input parameter of all of those function is a non-negative integer, indicating a position of a cube in the data structure.
   * In the case of functions that compute (co)boundary, the output is a vector if non-negative integers pointing to the positions of (co)boundary element of the input cell.
   **/
  inline std::vector< size_t > get_coboundary_of_a_cell(size_t cell_);
  /**
   * In the case of get_dimension_of_a_cell function, the output is a non-negative integer indicating the dimension of a cell.
   **/
  inline unsigned get_dimension_of_a_cell(size_t cell_);
  /**
   * In the case of get_cell_data, the output parameter is a reference to the value of a cube in a given position.
   **/
  inline T& get_cell_data(size_t cell_);


  /**
   * Typical input used to construct a baseBitmap class is a filtration given at the top dimensional cells. Then, there are a few ways one can pick the filtration of lower dimensional
   * cells. The most typical one is by so called lower star filtration. This function is always called by any constructor which takes the top dimensional cells. If you use such a constructor,
   * then there is no need to call this function. Call it only if you are putting the filtration of the cells by your own (for instance by using topDimensionalCellsIterator).
   **/
  void impose_lower_star_filtration(); //assume that top dimensional cells are already set.

  /**
   * Returns dimension of a complex.
   **/
  inline unsigned dimension() {
    return sizes.size();
  }

  /**
   * Returns number of all cubes in the data structure.
   **/
  inline unsigned size_of_bitmap() {
    return this->data.size();
  }

  /**
   * Writing to stream operator.
   **/
  template <typename K>
  friend ostream& operator<<(ostream & os_, const Bitmap_cubical_complex_base<K>& b_);

  //ITERATORS

  /**
   * Iterator through all cells in the complex (in order they appear in the structure -- i.e. in lexicographical order).
   **/
  typedef typename std::vector< T >::iterator all_cells_iterator;

  all_cells_iterator all_cells_begin()const {
    return this->data.begin();
  }

  all_cells_iterator all_cells_end()const {
    return this->data.end();
  }


  typedef typename std::vector< T >::const_iterator all_cells_const_iterator;

  all_cells_const_iterator all_cells_const_begin()const {
    return this->data.begin();
  }

  all_cells_const_iterator all_cells_const_end()const {
    return this->data.end();
  }

  /**
   * Iterator through top dimensional cells of the complex. The cells appear in order they are stored in the structure (i.e. in lexicographical order)
   **/
  class Top_dimensional_cells_iterator : std::iterator< std::input_iterator_tag, double > {
   public:

    Top_dimensional_cells_iterator(Bitmap_cubical_complex_base& b_) : b(b_) {
      for (size_t i = 0; i != b_.dimension(); ++i) {
        this->counter.push_back(0);
      }
    }

    Top_dimensional_cells_iterator operator++() {
      //first find first element of the counter that can be increased:
      size_t dim = 0;
      while ((dim != this->b.dimension()) && (this->counter[dim] == this->b.sizes[dim] - 1))++dim;

      if (dim != this->b.dimension()) {
        ++this->counter[dim];
        for (size_t i = 0; i != dim; ++i) {
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

    Top_dimensional_cells_iterator operator=(const Top_dimensional_cells_iterator& rhs_) {
      this->counter = rhs_.counter;
      this->b = rhs_.b;
      return *this;
    }

    bool operator==(const Top_dimensional_cells_iterator& rhs_) {
      if (&this->b != &rhs_.b)return false;
      if (this->counter.size() != rhs_.counter.size())return false;
      for (size_t i = 0; i != this->counter.size(); ++i) {
        if (this->counter[i] != rhs_.counter[i])return false;
      }
      return true;
    }

    bool operator!=(const Top_dimensional_cells_iterator& rhs_) {
      return !(*this == rhs_);
    }

    T& operator*() {
      //given the counter, compute the index in the array and return this element.
      unsigned index = 0;
      for (size_t i = 0; i != this->counter.size(); ++i) {
        index += (2 * this->counter[i] + 1) * this->b.multipliers[i];
      }
      return this->b.data[index];
    }

    size_t computeIndexInBitmap() {
      size_t index = 0;
      for (size_t i = 0; i != this->counter.size(); ++i) {
        index += (2 * this->counter[i] + 1) * this->b.multipliers[i];
      }
      return index;
    }

    void printCounter() {
      for (size_t i = 0; i != this->counter.size(); ++i) {
        cout << this->counter[i] << " ";
      }
    }
    friend class Bitmap_cubical_complex_base;
   protected:
    std::vector< unsigned > counter;
    Bitmap_cubical_complex_base& b;
  };

  Top_dimensional_cells_iterator top_dimensional_cells_begin() {
    Top_dimensional_cells_iterator a(*this);
    return a;
  }

  Top_dimensional_cells_iterator top_dimensional_cells_end() {
    Top_dimensional_cells_iterator a(*this);
    for (size_t i = 0; i != this->dimension(); ++i) {
      a.counter[i] = this->sizes[i] - 1;
    }
    a.counter[0]++;
    return a;
  }


  //****************************************************************************************************************//
  //****************************************************************************************************************//
  //****************************************************************************************************************//
  //****************************************************************************************************************//


  //****************************************************************************************************************//
  //****************************************************************************************************************//
  //****************************************************************************************************************//
  //****************************************************************************************************************//

 protected:
  std::vector<unsigned> sizes;
  std::vector<unsigned> multipliers;
  std::vector<T> data;
  size_t totalNumberOfCells;

  void set_up_containers(std::vector<unsigned> sizes_) {
    unsigned multiplier = 1;
    for (size_t i = 0; i != sizes_.size(); ++i) {
      this->sizes.push_back(sizes_[i]);
      this->multipliers.push_back(multiplier);
      //multiplier *= 2*(sizes[i]+1)+1;
      multiplier *= 2 * sizes_[i] + 1;
    }
    //std::reverse( this->sizes.begin() , this->sizes.end() );
    std::vector<T> data(multiplier);
    std::fill(data.begin(), data.end(), INT_MAX);
    this->totalNumberOfCells = multiplier;
    this->data = data;
  }

  size_t compute_position_in_bitmap(std::vector< int > counter_) {
    size_t position = 0;
    for (size_t i = 0; i != this->multipliers.size(); ++i) {
      position += this->multipliers[i] * counter_[i];
    }
    return position;
  }

  std::vector<unsigned> compute_counter_for_given_cell(size_t cell_) {
    std::vector<unsigned> counter;
    for (size_t dim = this->sizes.size(); dim != 0; --dim) {
      counter.push_back(cell_ / this->multipliers[dim - 1]);
      cell_ = cell_ % this->multipliers[dim - 1];
    }
    std::reverse(counter.begin(), counter.end());
    return counter;
  }

  std::vector< size_t > generate_vector_of_shifts_for_bitmaps_with_periodic_boundary_conditions(std::vector< bool > directionsForPeriodicBCond_);
};

template <typename K>
ostream& operator<<(ostream & out_, const Bitmap_cubical_complex_base<K>& b_) {
  //for ( typename bitmap<K>::all_cells_const_iterator it = b.all_cells_const_begin() ; it != b.all_cells_const_end() ; ++it )
  for (typename Bitmap_cubical_complex_base<K>::all_cells_const_iterator it = b_.all_cells_const_begin(); it != b_.all_cells_const_end(); ++it) {
    out_ << *it << " ";
  }
  return out_;
}

template <typename T>
Bitmap_cubical_complex_base<T>::Bitmap_cubical_complex_base(std::vector<unsigned> sizes_) {
  this->set_up_containers(sizes_);
}

template <typename T>
Bitmap_cubical_complex_base<T>::Bitmap_cubical_complex_base(std::vector<unsigned> sizesInFollowingDirections_, std::vector<T> topDimensionalCells_) {
  this->set_up_containers(sizesInFollowingDirections_);

  size_t numberOfTopDimensionalElements = 1;
  for (size_t i = 0; i != sizesInFollowingDirections_.size(); ++i) {
    numberOfTopDimensionalElements *= sizesInFollowingDirections_[i];
  }
  if (numberOfTopDimensionalElements != topDimensionalCells_.size()) {
    cerr << "Error in constructor Bitmap_cubical_complex_base( std::vector<size_t> sizesInFollowingDirections_ , std::vector<float> topDimensionalCells_ ). Number of top dimensional elements that follow from sizesInFollowingDirections vector is different than the size of topDimensionalCells vector." << endl;
    throw ("Error in constructor Bitmap_cubical_complex_base( std::vector<size_t> sizesInFollowingDirections_ , std::vector<float> topDimensionalCells_ ). Number of top dimensional elements that follow from sizesInFollowingDirections vector is different than the size of topDimensionalCells vector.");
  }

  Bitmap_cubical_complex_base<T>::Top_dimensional_cells_iterator it(*this);
  size_t index = 0;
  for (it = this->top_dimensional_cells_begin(); it != this->top_dimensional_cells_end(); ++it) {
    (*it) = topDimensionalCells_[index];
    ++index;
  }
  this->impose_lower_star_filtration();
}

template <typename T>
Bitmap_cubical_complex_base<T>::Bitmap_cubical_complex_base(char* perseusStyleFile_) {
  bool dbg = false;
  ifstream inFiltration, inIds;
  inFiltration.open(perseusStyleFile_);
  unsigned dimensionOfData;
  inFiltration >> dimensionOfData;

  if (dbg) {
    cerr << "dimensionOfData : " << dimensionOfData << endl;
  }

  std::vector<unsigned> sizes;
  for (size_t i = 0; i != dimensionOfData; ++i) {
    int sizeInThisDimension;
    inFiltration >> sizeInThisDimension;
    sizeInThisDimension = abs(sizeInThisDimension);
    sizes.push_back(sizeInThisDimension);
    if (dbg) {
      cerr << "sizeInThisDimension : " << sizeInThisDimension << endl;
    }
  }
  this->set_up_containers(sizes);

  Bitmap_cubical_complex_base<T>::Top_dimensional_cells_iterator it(*this);
  it = this->top_dimensional_cells_begin();

  //TODO -- over here we also need to read id's of cell and put them to bitmapElement structure!
  while (!inFiltration.eof()) {
    double filtrationLevel;
    inFiltration >> filtrationLevel;
    if (dbg) {
      cerr << "Cell of an index : " << it.computeIndexInBitmap() << " and dimension: " << this->get_dimension_of_a_cell(it.computeIndexInBitmap()) << " get the value : " << filtrationLevel << endl;
    }
    *it = filtrationLevel;
    ++it;
  }
  inFiltration.close();
  this->impose_lower_star_filtration();
}

template <typename T>
std::vector< size_t > Bitmap_cubical_complex_base<T>::get_boundary_of_a_cell(size_t cell_) {
  bool bdg = false;
  //first of all, we need to take the list of coordinates in which the cell has nonzero length. We do it by using modified version to compute dimension of a cell:
  std::vector< unsigned > dimensionsInWhichCellHasNonzeroLength;
  unsigned dimension = 0;
  size_t cell1 = cell_;
  for (size_t i = this->multipliers.size(); i != 0; --i) {
    unsigned position = cell1 / multipliers[i - 1];
    if (position % 2 == 1) {
      dimensionsInWhichCellHasNonzeroLength.push_back(i - 1);
      dimension++;
    }
    cell1 = cell1 % multipliers[i - 1];
  }

  if (bdg) {
    cerr << "dimensionsInWhichCellHasNonzeroLength : \n";
    for (size_t i = 0; i != dimensionsInWhichCellHasNonzeroLength.size(); ++i) {
      cerr << dimensionsInWhichCellHasNonzeroLength[i] << endl;
    }
    getchar();
  }

  std::vector< size_t > boundaryElements;
  if (dimensionsInWhichCellHasNonzeroLength.size() == 0)return boundaryElements;
  for (size_t i = 0; i != dimensionsInWhichCellHasNonzeroLength.size(); ++i) {
    boundaryElements.push_back(cell_ - multipliers[ dimensionsInWhichCellHasNonzeroLength[i] ]);
    boundaryElements.push_back(cell_ + multipliers[ dimensionsInWhichCellHasNonzeroLength[i] ]);

    if (bdg) cerr << "multipliers[dimensionsInWhichCellHasNonzeroLength[i]] : " << multipliers[dimensionsInWhichCellHasNonzeroLength[i]] << endl;
    if (bdg) cerr << "cell_ - multipliers[dimensionsInWhichCellHasNonzeroLength[i]] : " << cell_ - multipliers[dimensionsInWhichCellHasNonzeroLength[i]] << endl;
    if (bdg) cerr << "cell_ + multipliers[dimensionsInWhichCellHasNonzeroLength[i]] : " << cell_ + multipliers[dimensionsInWhichCellHasNonzeroLength[i]] << endl;
  }
  return boundaryElements;
}

template <typename T>
std::vector< size_t > Bitmap_cubical_complex_base<T>::get_coboundary_of_a_cell(size_t cell_) {
  bool bdg = false;
  //first of all, we need to take the list of coordinates in which the cell has nonzero length. We do it by using modified version to compute dimension of a cell:
  std::vector< unsigned > dimensionsInWhichCellHasZeroLength;
  unsigned dimension = 0;
  size_t cell1 = cell_;
  for (size_t i = this->multipliers.size(); i != 0; --i) {
    unsigned position = cell1 / multipliers[i - 1];
    if (position % 2 == 0) {
      dimensionsInWhichCellHasZeroLength.push_back(i - 1);
      dimension++;
    }
    cell1 = cell1 % multipliers[i - 1];
  }

  std::vector<unsigned> counter = this->compute_counter_for_given_cell(cell_);
  //reverse(counter.begin() , counter.end());

  if (bdg) {
    cerr << "dimensionsInWhichCellHasZeroLength : \n";
    for (size_t i = 0; i != dimensionsInWhichCellHasZeroLength.size(); ++i) {
      cerr << dimensionsInWhichCellHasZeroLength[i] << endl;
    }
    cerr << "\n counter : " << endl;
    for (size_t i = 0; i != counter.size(); ++i) {
      cerr << counter[i] << endl;
    }
    getchar();
  }

  std::vector< size_t > coboundaryElements;
  if (dimensionsInWhichCellHasZeroLength.size() == 0)return coboundaryElements;
  for (size_t i = 0; i != dimensionsInWhichCellHasZeroLength.size(); ++i) {
    if (bdg) {
      cerr << "Dimension : " << i << endl;
      if (counter[dimensionsInWhichCellHasZeroLength[i]] == 0) {
        cerr << "In dimension : " << i << " we cannot substract, since we will jump out of a Bitmap_cubical_complex_base \n";
      }
      if (counter[dimensionsInWhichCellHasZeroLength[i]] == 2 * this->sizes[dimensionsInWhichCellHasZeroLength[i]]) {
        cerr << "In dimension : " << i << " we cannot substract, since we will jump out of a Bitmap_cubical_complex_base \n";
      }
    }


    if ((cell_ > multipliers[dimensionsInWhichCellHasZeroLength[i]]) && (counter[dimensionsInWhichCellHasZeroLength[i]] != 0))
      //if ( counter[dimensionsInWhichCellHasZeroLength[i]] != 0 )
    {
      if (bdg)cerr << "Subtracting : " << cell_ - multipliers[dimensionsInWhichCellHasZeroLength[i]] << endl;
      coboundaryElements.push_back(cell_ - multipliers[dimensionsInWhichCellHasZeroLength[i]]);
    }
    if ((cell_ + multipliers[dimensionsInWhichCellHasZeroLength[i]] < this->data.size()) && (counter[dimensionsInWhichCellHasZeroLength[i]] != 2 * this->sizes[dimensionsInWhichCellHasZeroLength[i]]))
      //if ( counter[dimensionsInWhichCellHasZeroLength[i]] != 2*this->sizes[dimensionsInWhichCellHasZeroLength[i]] )
    {
      coboundaryElements.push_back(cell_ + multipliers[dimensionsInWhichCellHasZeroLength[i]]);
      if (bdg)cerr << "Adding : " << cell_ + multipliers[dimensionsInWhichCellHasZeroLength[i]] << endl;
    }
  }
  return coboundaryElements;
}

template <typename T>
unsigned Bitmap_cubical_complex_base<T>::get_dimension_of_a_cell(size_t cell_) {
  bool dbg = false;
  if (dbg)cerr << "\n\n\n Computing position o a cell of an index : " << cell_ << endl;
  unsigned dimension = 0;
  for (size_t i = this->multipliers.size(); i != 0; --i) {
    unsigned position = cell_ / multipliers[i - 1];

    if (dbg)cerr << "i-1 :" << i - 1 << endl;
    if (dbg)cerr << "cell_ : " << cell_ << endl;
    if (dbg)cerr << "position : " << position << endl;
    if (dbg)cerr << "multipliers[" << i - 1 << "] = " << multipliers[i - 1] << endl;
    if (dbg)getchar();

    if (position % 2 == 1) {
      if (dbg)cerr << "Nonzero length in this direction \n";
      dimension++;
    }
    cell_ = cell_ % multipliers[i - 1];
  }
  return dimension;
}

template <typename T>
T& Bitmap_cubical_complex_base<T>::get_cell_data(size_t cell_) {
  return this->data[cell_];
}

template <typename T>
void Bitmap_cubical_complex_base<T>::impose_lower_star_filtration() {
  bool dbg = false;

  //this vector will be used to check which elements have already been taken care of in imposing lower star filtration:
  std::vector<bool> isThisCellConsidered(this->data.size(), false);

  std::vector<size_t> indicesToConsider;
  //we assume here that we already have a filtration on the top dimensional cells and we have to extend it to lower ones.
  typename Bitmap_cubical_complex_base<T>::Top_dimensional_cells_iterator it(*this);
  for (it = this->top_dimensional_cells_begin(); it != this->top_dimensional_cells_end(); ++it) {
    indicesToConsider.push_back(it.computeIndexInBitmap());
  }

  while (indicesToConsider.size()) {
    if (dbg) {
      cerr << "indicesToConsider in this iteration \n";
      for (size_t i = 0; i != indicesToConsider.size(); ++i) {
        cout << indicesToConsider[i] << "  ";
      }
      getchar();
    }
    std::vector<size_t> newIndicesToConsider;
    for (size_t i = 0; i != indicesToConsider.size(); ++i) {
      std::vector<size_t> bd = this->get_boundary_of_a_cell(indicesToConsider[i]);
      for (size_t boundaryIt = 0; boundaryIt != bd.size(); ++boundaryIt) {
        if (this->data[ bd[boundaryIt] ] > this->data[ indicesToConsider[i] ]) {
          this->data[ bd[boundaryIt] ] = this->data[ indicesToConsider[i] ];
        }
        if (isThisCellConsidered[ bd[boundaryIt] ] == false) {
          newIndicesToConsider.push_back(bd[boundaryIt]);
          isThisCellConsidered[ bd[boundaryIt] ] = true;
        }
      }
    }
    indicesToConsider.swap(newIndicesToConsider);
  }
}

template <typename T>
std::vector< size_t > Bitmap_cubical_complex_base<T>::generate_vector_of_shifts_for_bitmaps_with_periodic_boundary_conditions(std::vector< bool > directionsForPeriodicBCond_) {
  bool dbg = false;
  if (this->sizes.size() != directionsForPeriodicBCond_.size())throw "directionsForPeriodicBCond_ vector size is different from the size of the bitmap. The program will now terminate \n";

  std::vector<int> sizes(this->sizes.size());
  for (size_t i = 0; i != this->sizes.size(); ++i)sizes[i] = 2 * this->sizes[i];

  counter c(sizes);

  std::vector< size_t > result;

  for (size_t i = 0; i != this->data.size(); ++i) {
    size_t position;
    if (!c.isFinal()) {
      position = i;
      //result.push_back( i );
    } else {
      std::vector< bool > finals = c.directionsOfFinals();
      bool jumpInPosition = false;
      for (size_t dir = 0; dir != finals.size(); ++dir) {
        if (finals[dir] == false)continue;
        if (directionsForPeriodicBCond_[dir]) {
          jumpInPosition = true;
        }
      }
      if (jumpInPosition == true) {
        //in this case this guy is final, so we need to find 'the opposite one'
        position = compute_position_in_bitmap(c.findOpposite(directionsForPeriodicBCond_));
      } else {
        position = i;
      }
    }
    result.push_back(position);
    if (dbg) {
      cerr << " position : " << position << endl;
      cerr << c << endl;
      getchar();
    }

    c.increment();
  }

  return result;
}

#endif  // BITMAP_CUBICAL_COMPLEX_BASE_H_
