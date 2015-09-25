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

#ifndef BITMAP_CUBICAL_COMPLEX_H_
#define BITMAP_CUBICAL_COMPLEX_H_

#include <gudhi/Bitmap_cubical_complex_base.h>

#include <limits>
#include <utility>  // for pair
#include <algorithm>  // for sort
#include <vector>  // for vector

// global variable, was used just for debugging.
bool globalDbg = false;

template <typename T = double>
class Bitmap_cubical_complex : public Bitmap_cubical_complex_base<T> {
 public:
  //******************************************************************************************************************//
  // Typedefs and typenames
  //******************************************************************************************************************//
  friend class Simplex_handle;
  typedef size_t Simplex_key;
  typedef T Filtration_value;


  //******************************************************************************************************************//
  // Simplex handle class
  //******************************************************************************************************************//

  /**
   * Handle of a cell, required for compatibility with the function to compute persistence in Gudhi. Elements of this
   * class are: the pointer to the bitmap B in which the considered cell is together with a position of this cell in B.
   * Given this data, one can get all the information about the considered cell.
   **/
  class Simplex_handle {
   public:
    Simplex_handle() {
      if (globalDbg) {
        std::cerr << "Simplex_handle()\n";
      }
      this->b = 0;
      this->position = 0;
    }

    Simplex_handle(Bitmap_cubical_complex<T>* b) {
      if (globalDbg) {
        std::cerr << "Simplex_handle(Bitmap_cubical_complex<T>* b)\n";
      }
      this->b = b;
      this->position = 0;
    }

    Simplex_handle(const Simplex_handle& org) : b(org.b) {
      if (globalDbg) {
        std::cerr << "Simplex_handle( const Simplex_handle& org )\n";
      }
      this->position = org.position;
    }

    Simplex_handle& operator=(const Simplex_handle& rhs) {
      if (globalDbg) {
        std::cerr << "Simplex_handle operator =  \n";
      }
      this->position = rhs.position;
      this->b = rhs.b;
      return *this;
    }

    Simplex_handle(Bitmap_cubical_complex<T>* b, Simplex_key position) {
      if (globalDbg) {
        std::cerr << "Simplex_handle(Bitmap_cubical_complex<T>* b , Simplex_key position)\n";
        std::cerr << "Position : " << position << std::endl;
      }
      this->b = b;
      this->position = position;
    }
    friend class Bitmap_cubical_complex<T>;
   private:
    Bitmap_cubical_complex<T>* b;
    Simplex_key position;
    // Assumption -- this field always keep the REAL position of simplex in the bitmap, no matter what keys have been.
    // to deal with the keys, the class Bitmap_cubical_complex have extra vectors: keyAssociatedToSimplex and
    // simplexAssociatedToKey that allow to move between actual cell and the key assigned to it.
  };


  //******************************************************************************************************************//
  // Constructors
  //******************************************************************************************************************//
  // Over here we need to definie various input types. I am proposing the following ones:
  // Perseus style
  // TODO(Pawel Dlotko): H5 files?
  // TODO(Pawel Dlotko): binary files with little endiangs / big endians?
  // TODO(Pawel Dlotko): constructor from a vector of elements of a type T.

  /**
   * Constructor form a Perseus-style file.
   **/
  Bitmap_cubical_complex(char* perseusStyleFile) : Bitmap_cubical_complex_base<T>(perseusStyleFile) {
    if (globalDbg) {
      std::cerr << "Bitmap_cubical_complex( char* perseusStyleFile )\n";
    }
    std::vector< size_t > keyAssociatedToSimplex(this->totalNumberOfCells + 1);
    std::vector< size_t > simplexAssociatedToKey(this->totalNumberOfCells + 1);

    for (size_t i = 0; i != this->totalNumberOfCells; ++i) {
      keyAssociatedToSimplex[i] = simplexAssociatedToKey[i] = i;
    }
    this->keyAssociatedToSimplex = keyAssociatedToSimplex;
    this->simplexAssociatedToKey = simplexAssociatedToKey;
    // we initialize this only once, in each constructor, when the bitmap is constructed. If the user decide to change
    // some elements of the bitmap, then this procedure need to be called again.
    this->initializeElementsOrderedAccordingToFiltration();
  }

  /**
   * Constructor that requires vector of elements of type unsigned, which gives number of top dimensional cells in the
   * following directions and vector of element of a type T with filtration on top dimensional cells.
   **/
  Bitmap_cubical_complex(std::vector<unsigned> dimensions, std::vector<T> topDimensionalCells)
  : Bitmap_cubical_complex_base<T>(dimensions, topDimensionalCells) {
    std::vector< size_t > keyAssociatedToSimplex(this->totalNumberOfCells + 1);
    std::vector< size_t > simplexAssociatedToKey(this->totalNumberOfCells + 1);

    for (size_t i = 0; i != this->totalNumberOfCells; ++i) {
      keyAssociatedToSimplex[i] = simplexAssociatedToKey[i] = i;
    }
    this->keyAssociatedToSimplex = keyAssociatedToSimplex;
    this->simplexAssociatedToKey = simplexAssociatedToKey;
    // we initialize this only once, in each constructor, when the bitmap is constructed. If the user decide to change
    // some elements of the bitmap, then this procedure need to be called again.
    this->initializeElementsOrderedAccordingToFiltration();
  }

  //******************************************************************************************************************//
  // Other 'easy' functions
  //******************************************************************************************************************//

  /**
   * Returns number of all cubes in the complex.
   **/
  size_t num_simplices()const {
    return this->totalNumberOfCells;
  }

  /**
   * Returns a Simplex_handle to a cube that do not exist in this complex.
   **/
  Simplex_handle null_simplex() {
    return Simplex_handle(this, this->data.size());
  }

  /**
   * Returns dimension of the complex.
   **/
  size_t dimension() {
    return this->sizes.size();
  }

  /**
   * Return dimension of a cell pointed by the Simplex_handle.
   **/
  size_t dimension(const Simplex_handle& sh) {
    if (globalDbg) {
      std::cerr << "int dimension(const Simplex_handle& sh)\n";
    }
    if (sh.position != this->data.size()) return sh.b->get_dimension_of_a_cell(sh.position);
    return std::numeric_limits<size_t>::max();
  }

  /**
   * Return the filtration of a cell pointed by the Simplex_handle.
   **/
  T filtration(const Simplex_handle& sh) {
    if (globalDbg) {
      std::cerr << "T filtration(const Simplex_handle& sh)\n";
    }
    // Returns the filtration value of a simplex.
    if (sh.position != this->data.size()) return sh.b->data[ sh.position ];
    return INT_MAX;
  }

  /**
   * Return a key which is not a key of any cube in the considered data structure.
   **/
  Simplex_key null_key() {
    if (globalDbg) {
      std::cerr << "Simplex_key null_key()\n";
    }
    return this->data.size();
  }

  /**
   * Return the key of a cube pointed by the Simplex_handle.
   **/
  Simplex_key key(const Simplex_handle& sh) {
    if (globalDbg) {
      std::cerr << "Simplex_key key(const Simplex_handle& sh)\n";
    }
    return sh.b->keyAssociatedToSimplex[ sh.position ];
  }

  /**
   * Return the Simplex_handle given the key of the cube.
   **/
  Simplex_handle simplex(Simplex_key key) {
    if (globalDbg) {
      std::cerr << "Simplex_handle simplex(Simplex_key key)\n";
    }
    return Simplex_handle(this, this->simplexAssociatedToKey[ key ]);
  }

  /**
   * Assign key to a cube pointed by the Simplex_handle
   **/
  void assign_key(Simplex_handle& sh, Simplex_key key) {
    if (globalDbg) {
      std::cerr << "void assign_key(Simplex_handle& sh, Simplex_key key)\n";
    }
    this->keyAssociatedToSimplex[sh.position] = key;
    this->simplexAssociatedToKey[key] = sh.position;
  }

  /**
   * Function called from a constructor. It is needed for Filtration_simplex_iterator to work.
   **/
  void initializeElementsOrderedAccordingToFiltration();



  //******************************************************************************************************************//
  // Iterators
  //******************************************************************************************************************//

  /**
   * Boundary_simplex_iterator class allows iteration on boundary of each cube.
   **/
  class Boundary_simplex_range;

  class Boundary_simplex_iterator : std::iterator< std::input_iterator_tag, Simplex_handle > {
  // Iterator on the simplices belonging to the boundary of a simplex.
  // value_type must be 'Simplex_handle'.
   public:
    Boundary_simplex_iterator(Simplex_handle& sh) : sh(sh) {
      if (globalDbg) {
        std::cerr << "Boundary_simplex_iterator( Simplex_handle& sh )\n";
      }
      this->position = 0;
      this->boundaryElements = this->sh.b->get_boundary_of_a_cell(this->sh.position);
    }

    Boundary_simplex_iterator operator++() {
      if (globalDbg) {
        std::cerr << "Boundary_simplex_iterator operator++()\n";
      }
      ++this->position;
      return *this;
    }

    Boundary_simplex_iterator operator++(int) {
      Boundary_simplex_iterator result = *this;
      ++(*this);
      return result;
    }

    Boundary_simplex_iterator operator=(const Boundary_simplex_iterator& rhs) {
      if (globalDbg) {
        std::cerr << "Boundary_simplex_iterator operator =\n";
      }
      this->sh = rhs.sh;
      this->boundaryElements.clear();
      this->boundaryElementsinsert(this->boundaryElements.end(),
                                   rhs.boundaryElements.begin(), rhs.boundaryElements.end());
    }

    bool operator==(const Boundary_simplex_iterator& rhs) {
      if (globalDbg) {
        std::cerr << "bool operator ==\n";
      }
      if (this->position == rhs.position) {
        if (this->boundaryElements.size() != rhs.boundaryElements.size())return false;
        for (size_t i = 0; i != this->boundaryElements.size(); ++i) {
          if (this->boundaryElements[i] != rhs.boundaryElements[i])return false;
        }
        return true;
      }
      return false;
    }

    bool operator!=(const Boundary_simplex_iterator& rhs) {
      if (globalDbg) {
        std::cerr << "bool operator != \n";
      }
      return !(*this == rhs);
    }

    Simplex_handle operator*() {
      if (globalDbg) {
        std::cerr << "Simplex_handle operator*\n";
      }
      return Simplex_handle(this->sh.b, this->boundaryElements[this->position]);
    }

    friend class Boundary_simplex_range;
   private:
    Simplex_handle sh;
    std::vector< size_t > boundaryElements;
    size_t position;
  };

  /**
   * Boundary_simplex_range class provides ranges for boundary iterators.
   **/
  class Boundary_simplex_range {
  // Range giving access to the simplices in the boundary of a simplex.
  // .begin() and .end() return type Boundary_simplex_iterator.
   public:
    Boundary_simplex_range(const Simplex_handle& sh) : sh(sh) { }

    Boundary_simplex_iterator begin() {
      if (globalDbg) {
        std::cerr << "Boundary_simplex_iterator begin\n";
      }
      Boundary_simplex_iterator it(this->sh);
      return it;
    }

    Boundary_simplex_iterator end() {
      if (globalDbg) {
        std::cerr << "Boundary_simplex_iterator end()\n";
      }
      Boundary_simplex_iterator it(this->sh);
      it.position = it.boundaryElements.size();
      return it;
    }

   private:
    Simplex_handle sh;
  };


  /**
   * Filtration_simplex_iterator class provides an iterator though the whole structure in the order of filtration. Secondary criteria for filtration are:
   * (1) Dimension of a cube (lower dimensional comes first).
   * (2) Position in the data structure (the ones that are earlies in the data structure comes first).
   **/
  class Filtration_simplex_range;

  class Filtration_simplex_iterator : std::iterator< std::input_iterator_tag, Simplex_handle > {
  // Iterator over all simplices of the complex in the order of the indexing scheme.
  // 'value_type' must be 'Simplex_handle'.
   public:
    Filtration_simplex_iterator(Bitmap_cubical_complex* b) : b(b), position(0) { }

    Filtration_simplex_iterator() : b(NULL) { }

    Filtration_simplex_iterator operator++() {
      if (globalDbg) {
        std::cerr << "Filtration_simplex_iterator operator++\n";
      }
      ++this->position;
      return (*this);
    }

    Filtration_simplex_iterator operator++(int) {
      Filtration_simplex_iterator result = *this;
      ++(*this);
      return result;
    }

    Filtration_simplex_iterator operator=(const Filtration_simplex_iterator& rhs) {
      if (globalDbg) {
        std::cerr << "Filtration_simplex_iterator operator =\n";
      }
      this->b = rhs.b;
      this->position = rhs.position;
    }

    bool operator==(const Filtration_simplex_iterator& rhs) {
      if (globalDbg) {
        std::cerr << "bool operator == ( const Filtration_simplex_iterator& rhs )\n";
      }
      if (this->position == rhs.position) {
        return true;
      }
      return false;
    }

    bool operator!=(const Filtration_simplex_iterator& rhs) {
      if (globalDbg) {
        std::cerr << "bool operator != ( const Filtration_simplex_iterator& rhs )\n";
      }
      return !(*this == rhs);
    }

    Simplex_handle operator*() {
      if (globalDbg) {
        std::cerr << "Simplex_handle operator*()\n";
      }
      return Simplex_handle(this->b, this->b->elementsOrderedAccordingToFiltration[ this->position ]);
    }

    friend class Filtration_simplex_range;
   private:
    Bitmap_cubical_complex<T>* b;
    size_t position;
  };

  /**
   * Filtration_simplex_range provides the ranges for Filtration_simplex_iterator.
   **/
  class Filtration_simplex_range {
  // Range over the simplices of the complex in the order of the filtration.
  // .begin() and .end() return type Filtration_simplex_iterator.
   public:
    Filtration_simplex_range(Bitmap_cubical_complex<T>* b) : b(b) { }

    Filtration_simplex_iterator begin() {
      if (globalDbg) {
        std::cerr << "Filtration_simplex_iterator begin() \n";
      }
      return Filtration_simplex_iterator(this->b);
    }

    Filtration_simplex_iterator end() {
      if (globalDbg) {
        std::cerr << "Filtration_simplex_iterator end()\n";
      }
      Filtration_simplex_iterator it(this->b);
      it.position = this->b->elementsOrderedAccordingToFiltration.size();
      return it;
    }
   private:
    Bitmap_cubical_complex<T>* b;
  };



  //******************************************************************************************************************//
  // Methods to access iterators from the container:

  /**
   * boundary_simplex_range creates an object of a Boundary_simplex_range class that provides ranges for the Boundary_simplex_iterator.
   **/
  Boundary_simplex_range boundary_simplex_range(Simplex_handle& sh) {
    if (globalDbg) {
      std::cerr << "Boundary_simplex_range boundary_simplex_range(Simplex_handle& sh)\n";
    }
    // Returns a range giving access to all simplices of the boundary of a simplex, i.e. the set of
    // codimension 1 subsimplices of the Simplex.
    return Boundary_simplex_range(sh);
  }

  /**
   * filtration_simplex_range creates an object of a Filtration_simplex_range class that provides ranges for the
   * Filtration_simplex_iterator.
   **/
  Filtration_simplex_range filtration_simplex_range() {
    if (globalDbg) {
      std::cerr << "Filtration_simplex_range filtration_simplex_range()\n";
    }
    // Returns a range over the simplices of the complex in the order of the filtration
    return Filtration_simplex_range(this);
  }
  //******************************************************************************************************************//



  //******************************************************************************************************************//
  // Elements which are in Gudhi now, but I (and in all the cases I asked also Marc) do not understand why they are
  // there.
  // TODO(Pawel Dlotko): The file IndexingTag.h in the Gudhi library contains an empty structure, so I understand that
  // this is something that was planned (for simplicial maps?) but was never finished. The only idea I have here is
  // to use the same empty structure from IndexingTag.h file, but only if the compiler needs it. If the compiler
  // do not need it, then I would rather not add here elements which I do not understand.
  // typedef Indexing_tag

  /**
   * Function needed for compatibility with Gudhi. Not useful for other purposes.
   **/
  std::pair<Simplex_handle, Simplex_handle> endpoints(Simplex_handle sh) {
    std::vector< size_t > bdry = this->get_boundary_of_a_cell(sh.position);
    if (globalDbg) {
      std::cerr << "std::pair<Simplex_handle, Simplex_handle> endpoints( Simplex_handle sh )\n";
      std::cerr << "bdry.size() : " << bdry.size() << std::endl;
    }
    // this method returns two first elements from the boundary of sh.
    if (bdry.size() < 2)
      throw("Error in endpoints in Bitmap_cubical_complex class. "
             "The cell for which this method was called have less than two elements in the boundary.");
    return std::make_pair(Simplex_handle(this, bdry[0]), Simplex_handle(this, bdry[1]));
  }


  /**
   * Class needed for compatibility with Gudhi. Not useful for other purposes.
   **/
  class Skeleton_simplex_range;

  class Skeleton_simplex_iterator : std::iterator< std::input_iterator_tag, Simplex_handle > {
    // Iterator over all simplices of the complex in the order of the indexing scheme.
    // 'value_type' must be 'Simplex_handle'.
   public:
    Skeleton_simplex_iterator(Bitmap_cubical_complex* b, size_t d) : b(b), dimension(d) {
      if (globalDbg) {
        std::cerr << "Skeleton_simplex_iterator ( Bitmap_cubical_complex* b , size_t d )\n";
      }
      // find the position of the first simplex of a dimension d
      this->position = 0;
      while ((this->position != b->data.size()) &&
             (this->b->get_dimension_of_a_cell(this->position) != this->dimension)) {
        ++this->position;
      }
    }

    Skeleton_simplex_iterator() : b(NULL), dimension(0) { }

    Skeleton_simplex_iterator operator++() {
      if (globalDbg) {
        std::cerr << "Skeleton_simplex_iterator operator++()\n";
      }
      // increment the position as long as you did not get to the next element of the dimension dimension.
      ++this->position;
      while ((this->position != this->b->data.size()) &&
             (this->b->get_dimension_of_a_cell(this->position) != this->dimension)) {
        ++this->position;
      }
      return (*this);
    }

    Skeleton_simplex_iterator operator++(int) {
      Skeleton_simplex_iterator result = *this;
      ++(*this);
      return result;
    }

    Skeleton_simplex_iterator operator=(const Skeleton_simplex_iterator& rhs) {
      if (globalDbg) {
        std::cerr << "Skeleton_simplex_iterator operator =\n";
      }
      this->b = rhs.b;
      this->position = rhs.position;
    }

    bool operator==(const Skeleton_simplex_iterator& rhs) {
      if (globalDbg) {
        std::cerr << "bool operator ==\n";
      }
      if (this->position == rhs.position) {
        return true;
      }
      return false;
    }

    bool operator!=(const Skeleton_simplex_iterator& rhs) {
      if (globalDbg) {
        std::cerr << "bool operator != ( const Skeleton_simplex_iterator& rhs )\n";
      }
      return !(*this == rhs);
    }

    Simplex_handle operator*() {
      if (globalDbg) {
        std::cerr << "Simplex_handle operator*() \n";
      }
      return Simplex_handle(this->b, this->position);
    }

    friend class Skeleton_simplex_range;
   private:
    Bitmap_cubical_complex<T>* b;
    size_t position;
    int dimension;
  };

  /**
   * Class needed for compatibility with Gudhi. Not useful for other purposes.
   **/
  class Skeleton_simplex_range {
    // Range over the simplices of the complex in the order of the filtration.
    // .begin() and .end() return type Filtration_simplex_iterator.
   public:
    Skeleton_simplex_range(Bitmap_cubical_complex<T>* b, int dimension) : b(b), dimension(dimension) { }

    Skeleton_simplex_iterator begin() {
      if (globalDbg) {
        std::cerr << "Skeleton_simplex_iterator begin()\n";
      }
      return Skeleton_simplex_iterator(this->b, this->dimension);
    }

    Skeleton_simplex_iterator end() {
      if (globalDbg) {
        std::cerr << "Skeleton_simplex_iterator end()\n";
      }
      Skeleton_simplex_iterator it(this->b, this->dimension);
      it.position = this->b->data.size();
      return it;
    }

   private:
    Bitmap_cubical_complex<T>* b;
    int dimension;
  };

  /**
   * Function needed for compatibility with Gudhi. Not useful for other purposes.
   **/
  Skeleton_simplex_range skeleton_simplex_range(int dimension) {
    if (globalDbg) {
      std::cerr << "Skeleton_simplex_range skeleton_simplex_range( int dimension )\n";
    }
    return Skeleton_simplex_range(this, dimension);
  }



  //******************************************************************************************************************//
  // functions used for debugging:

  /**
   * Function used for debugging purposes.
   **/
  void printKeyAssociatedToSimplex() {
    for (size_t i = 0; i != this->data.size(); ++i) {
      std::cerr << i << " -> " << this->simplexAssociatedToKey[i] << std::endl;
    }
  }

  /**
   * Function used for debugging purposes.
   **/
  size_t printRealPosition(const Simplex_handle& sh) {
    return sh.position;
  }

 private:
  std::vector< size_t > keyAssociatedToSimplex;
  std::vector< size_t > simplexAssociatedToKey;
  // needed by Filtration_simplex_iterator. If this iterator is not used, this field is not initialized.
  std::vector< size_t > elementsOrderedAccordingToFiltration;
};

template <typename T>
bool compareElementsForElementsOrderedAccordingToFiltration(const std::pair< size_t,
                                                            std::pair< T, char > >& f,
                                                            const std::pair< size_t,
                                                            std::pair< T, char > >& s) {
  if (globalDbg) {
    std::cerr << "ompareElementsForElementsOrderedAccordingToFiltration\n";
  }
  if (f.second.first < s.second.first) {
    return true;
  } else {
    if (f.second.first > s.second.first) {
      return false;
    } else {
      // in this case f.second.first == s.second.first, and we use dimension to compare:
      if (f.second.second < s.second.second) {
        return true;
      } else {
        if (f.second.second > s.second.second) {
          return false;
        } else {
          // in this case, both the filtration value and the dimensions for those cells are the same.
          // Since it may be nice to have a stable sorting procedure, in this case, we compare positions in the bitmap:
          return ( f.first < s.first);
        }
      }
    }
  }
}

template <typename T>
void Bitmap_cubical_complex<T>::initializeElementsOrderedAccordingToFiltration() {
  if (globalDbg) {
    std::cerr << "void Bitmap_cubical_complex<T>::initializeElementsOrderedAccordingToFiltration() \n";
  }
  // ( position , (filtration , dimension) )
  std::vector< std::pair< size_t, std::pair< T, char > > > dataOfElementsFromBitmap(this->data.size());
  for (size_t i = 0; i != this->data.size(); ++i) {
    // TODO(Pawel Dlotko): This can be optimized by having a counter here. We do not need to re-compute the dimension
    // for every cell from scratch
    dataOfElementsFromBitmap[i] = std::make_pair(i, std::make_pair(this->data[i], this->get_dimension_of_a_cell(i)));
  }
  std::sort(dataOfElementsFromBitmap.begin(), dataOfElementsFromBitmap.end(),
            compareElementsForElementsOrderedAccordingToFiltration<T>);

  // Elements of bitmap ordered according to filtration then according to dimension then according to position in bitmap
  std::vector< size_t > elements_of_bitmap_ordered(this->data.size());
  for (size_t i = 0; i != dataOfElementsFromBitmap.size(); ++i) {
    elements_of_bitmap_ordered[i] = dataOfElementsFromBitmap[i].first;
  }
  this->elementsOrderedAccordingToFiltration = elements_of_bitmap_ordered;
}


//****************************************************************************************************************//
//****************************************************************************************************************//
//****************************************************************************************************************//
//****************************************************************************************************************//

#endif  // BITMAP_CUBICAL_COMPLEX_H_
