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
#include <gudhi/Bitmap_cubical_complex_periodic_boundary_conditions_base.h>

#include <limits>
#include <utility>  // for pair<>
#include <algorithm>  // for sort
#include <vector>

namespace Gudhi {

namespace Cubical_complex {

// global variable, was used just for debugging.
const bool globalDbg = false;

template <typename T> class is_before_in_filtration;

/**
* This is a Bitmap_cubical_complex class. It joints a functionalities of Bitmap_cubical_complex_base and Bitmap_cubical_complex_periodic_boundary_conditions_base classes into
* Gudhi persistent homology engine. It is a template class that inherit from its template parameter. The template parameter is supposed to be either Bitmap_cubical_complex_base or Bitmap_cubical_complex_periodic_boundary_conditions_base class.
**/

/**
 *@class Bitmap_cubical_complex
 *@brief Cubical complex represented as a bitmap.
 *@ingroup cubical_complex
 */
template <typename T>
class Bitmap_cubical_complex : public T {
 public:
  //*********************************************//
  // Typedefs and typenames
  //*********************************************//
  typedef size_t Simplex_key;
  typedef typename T::filtration_type Filtration_value;
  typedef Simplex_key Simplex_handle;


  //*********************************************//
  // Constructors
  //*********************************************//
  // Over here we need to define various input types. I am proposing the following ones:
  // Perseus style
  // TODO(PD) H5 files?
  // TODO(PD) binary files with little endiangs / big endians ?
  // TODO(PD) constructor from a vector of elements of a type T. ?

  /**
   * Constructor form a Perseus-style file.
   **/
  Bitmap_cubical_complex(const char* perseus_style_file) :
      T(perseus_style_file), key_associated_to_simplex(this->total_number_of_cells + 1) {
    if (globalDbg) {
      std::cerr << "Bitmap_cubical_complex( const char* perseus_style_file )\n";
    }
    for (size_t i = 0; i != this->total_number_of_cells; ++i) {
      this->key_associated_to_simplex[i] = i;
    }
    // we initialize this only once, in each constructor, when the bitmap is constructed.
    // If the user decide to change some elements of the bitmap, then this procedure need
    // to be called again.
    this->initialize_simplex_associated_to_key();
  }

  /**
   * Constructor that requires vector of elements of type unsigned, which gives number of top dimensional cells
   * in the following directions and vector of element of a type T
   * with filtration on top dimensional cells.
   **/
  Bitmap_cubical_complex(const std::vector<unsigned>& dimensions,
                         const std::vector<typename T::filtration_type>& top_dimensional_cells) :
      T(dimensions, top_dimensional_cells),
      key_associated_to_simplex(this->total_number_of_cells + 1) {
    for (size_t i = 0; i != this->total_number_of_cells; ++i) {
      this->key_associated_to_simplex[i] = i;
    }
    // we initialize this only once, in each constructor, when the bitmap is constructed.
    // If the user decide to change some elements of the bitmap, then this procedure need
    // to be called again.
    this->initialize_simplex_associated_to_key();
  }

  /**
   * Constructor that requires vector of elements of type unsigned, which gives number of top dimensional cells
   * in the following directions and vector of element of a type T::filtration_type
   * with filtration on top dimensional cells. The last parameter of the constructor is a vector of bools of a length equal to the dimension of cubical complex.
   * If the position i on this vector is true, then we impose periodic boundary conditions in this direction.
   **/
  Bitmap_cubical_complex(const std::vector<unsigned>& dimensions,
                         const std::vector<typename T::filtration_type>& top_dimensional_cells,
                         std::vector< bool > directions_in_which_periodic_b_cond_are_to_be_imposed) :
      T(dimensions, top_dimensional_cells, directions_in_which_periodic_b_cond_are_to_be_imposed),
      key_associated_to_simplex(this->total_number_of_cells + 1) {
    for (size_t i = 0; i != this->total_number_of_cells; ++i) {
      this->key_associated_to_simplex[i] = i;
    }
    // we initialize this only once, in each constructor, when the bitmap is constructed.
    // If the user decide to change some elements of the bitmap, then this procedure need
    // to be called again.
    this->initialize_simplex_associated_to_key();
  }

  /**
   * Destructor of the Bitmap_cubical_complex class.
   **/
  virtual ~Bitmap_cubical_complex(){}

  //*********************************************//
  // Other 'easy' functions
  //*********************************************//

  /**
   * Returns number of all cubes in the complex.
   **/
  size_t num_simplices()const {
    return this->total_number_of_cells;
  }

  /**
   * Returns a Simplex_handle to a cube that do not exist in this complex.
   **/
  static Simplex_handle null_simplex() {
    if (globalDbg) {
      std::cerr << "Simplex_handle null_simplex()\n";
    }
    return std::numeric_limits<Simplex_handle>::max();
  }

  /**
   * Returns dimension of the complex.
   **/
  inline size_t dimension()const {
    return this->sizes.size();
  }

  /**
   * Return dimension of a cell pointed by the Simplex_handle.
   **/
  inline unsigned dimension(Simplex_handle sh)const {
    if (globalDbg) {
      std::cerr << "unsigned dimension(const Simplex_handle& sh)\n";
    }
    if (sh != std::numeric_limits<Simplex_handle>::max()) return this->get_dimension_of_a_cell(sh);
    return -1;
  }

  /**
   * Return the filtration of a cell pointed by the Simplex_handle.
   **/
  typename T::filtration_type filtration(Simplex_handle sh) {
    if (globalDbg) {
      std::cerr << "T::filtration_type filtration(const Simplex_handle& sh)\n";
    }
    // Returns the filtration value of a simplex.
    if (sh != std::numeric_limits<Simplex_handle>::max()) return this->data[sh];
    return std::numeric_limits<Simplex_handle>::max();
  }

  /**
   * Return a key which is not a key of any cube in the considered data structure.
   **/
  static Simplex_key null_key() {
    if (globalDbg) {
      std::cerr << "Simplex_key null_key()\n";
    }
    return std::numeric_limits<Simplex_handle>::max();
  }

  /**
   * Return the key of a cube pointed by the Simplex_handle.
   **/
  Simplex_key key(Simplex_handle sh)const {
    if (globalDbg) {
      std::cerr << "Simplex_key key(const Simplex_handle& sh)\n";
    }
    if (sh != std::numeric_limits<Simplex_handle>::max()) {
      return this->key_associated_to_simplex[sh];
    }
    return this->null_key();
  }

  /**
   * Return the Simplex_handle given the key of the cube.
   **/
  Simplex_handle simplex(Simplex_key key) {
    if (globalDbg) {
      std::cerr << "Simplex_handle simplex(Simplex_key key)\n";
    }
    if (key != std::numeric_limits<Simplex_handle>::max()) {
      return this->simplex_associated_to_key[ key ];
    }
    return null_simplex();
  }

  /**
   * Assign key to a cube pointed by the Simplex_handle
   **/
  void assign_key(Simplex_handle sh, Simplex_key key) {
    if (globalDbg) {
      std::cerr << "void assign_key(Simplex_handle& sh, Simplex_key key)\n";
    }
    if (key == std::numeric_limits<Simplex_handle>::max()) return;
    this->key_associated_to_simplex[sh] = key;
    this->simplex_associated_to_key[key] = sh;
  }

  /**
   * Function called from a constructor. It is needed for Filtration_simplex_iterator to work.
   **/
  void initialize_simplex_associated_to_key();

  //*********************************************//
  // Iterators
  //*********************************************//

  /**
   * Boundary_simplex_range class provides ranges for boundary iterators.
   **/
  typedef typename std::vector< Simplex_handle >::iterator Boundary_simplex_iterator;
  typedef typename std::vector< Simplex_handle > Boundary_simplex_range;

  /**
   * Filtration_simplex_iterator class provides an iterator though the whole structure in the order of filtration.
   * Secondary criteria for filtration are:
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

    Filtration_simplex_iterator& operator=(const Filtration_simplex_iterator& rhs) {
      if (globalDbg) {
        std::cerr << "Filtration_simplex_iterator operator =\n";
      }
      this->b = rhs.b;
      this->position = rhs.position;
    }

    bool operator==(const Filtration_simplex_iterator& rhs)const {
      if (globalDbg) {
        std::cerr << "bool operator == ( const Filtration_simplex_iterator& rhs )\n";
      }
      return ( this->position == rhs.position);
    }

    bool operator!=(const Filtration_simplex_iterator& rhs)const {
      if (globalDbg) {
        std::cerr << "bool operator != ( const Filtration_simplex_iterator& rhs )\n";
      }
      return !(*this == rhs);
    }

    Simplex_handle operator*() {
      if (globalDbg) {
        std::cerr << "Simplex_handle operator*()\n";
      }
      return this->b->simplex_associated_to_key[ this->position ];
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
    typedef Filtration_simplex_iterator const_iterator;
    typedef Filtration_simplex_iterator iterator;

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
      it.position = this->b->simplex_associated_to_key.size();
      return it;
    }

   private:
    Bitmap_cubical_complex<T>* b;
  };



  //*********************************************//
  // Methods to access iterators from the container:

  /**
   * boundary_simplex_range creates an object of a Boundary_simplex_range class
   * that provides ranges for the Boundary_simplex_iterator.
   **/
  Boundary_simplex_range boundary_simplex_range(Simplex_handle sh) {
    return this->get_boundary_of_a_cell(sh);
  }

  /**
   * filtration_simplex_range creates an object of a Filtration_simplex_range class
   * that provides ranges for the Filtration_simplex_iterator.
   **/
  Filtration_simplex_range filtration_simplex_range() {
    if (globalDbg) {
      std::cerr << "Filtration_simplex_range filtration_simplex_range()\n";
    }
    // Returns a range over the simplices of the complex in the order of the filtration
    return Filtration_simplex_range(this);
  }
  //*********************************************//



  //*********************************************//
  // Elements which are in Gudhi now, but I (and in all the cases I asked also Marc) do not understand why they are
  // there.
  // TODO(PD) the file IndexingTag.h in the Gudhi library contains an empty structure, so
  // I understand that this is something that was planned (for simplicial maps?)
  // but was never finished. The only idea I have here is to use the same empty structure from
  // IndexingTag.h file, but only if the compiler needs it. If the compiler
  // do not need it, then I would rather not add here elements which I do not understand.
  // typedef Indexing_tag

  /**
   * Function needed for compatibility with Gudhi. Not useful for other purposes.
   **/
  std::pair<Simplex_handle, Simplex_handle> endpoints(Simplex_handle sh) {
    std::vector< size_t > bdry = this->get_boundary_of_a_cell(sh);
    if (globalDbg) {
      std::cerr << "std::pair<Simplex_handle, Simplex_handle> endpoints( Simplex_handle sh )\n";
      std::cerr << "bdry.size() : " << bdry.size() << std::endl;
    }
    // this method returns two first elements from the boundary of sh.
    if (bdry.size() < 2)
      throw("Error in endpoints in Bitmap_cubical_complex class. The cell have less than two elements in the boundary.");
    return std::make_pair(bdry[0], bdry[1]);
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
      while (
             (this->position != b->data.size()) &&
             (this->b->get_dimension_of_a_cell(this->position) != this->dimension)
             ) {
        ++this->position;
      }
    }

    Skeleton_simplex_iterator() : b(NULL), position(0), dimension(0) { }

    Skeleton_simplex_iterator operator++() {
      if (globalDbg) {
        std::cerr << "Skeleton_simplex_iterator operator++()\n";
      }
      // increment the position as long as you did not get to the next element of the dimension dimension.
      ++this->position;
      while (
             (this->position != this->b->data.size()) &&
             (this->b->get_dimension_of_a_cell(this->position) != this->dimension)
             ) {
        ++this->position;
      }
      return (*this);
    }

    Skeleton_simplex_iterator operator++(int) {
      Skeleton_simplex_iterator result = *this;
      ++(*this);
      return result;
    }

    Skeleton_simplex_iterator& operator=(const Skeleton_simplex_iterator& rhs) {
      if (globalDbg) {
        std::cerr << "Skeleton_simplex_iterator operator =\n";
      }
      this->b = rhs.b;
      this->position = rhs.position;
      this->dimension = rhs.dimension;
    }

    bool operator==(const Skeleton_simplex_iterator& rhs)const {
      if (globalDbg) {
        std::cerr << "bool operator ==\n";
      }
      return ( this->position == rhs.position);
    }

    bool operator!=(const Skeleton_simplex_iterator& rhs)const {
      if (globalDbg) {
        std::cerr << "bool operator != ( const Skeleton_simplex_iterator& rhs )\n";
      }
      return !(*this == rhs);
    }

    Simplex_handle operator*() {
      if (globalDbg) {
        std::cerr << "Simplex_handle operator*() \n";
      }
      return this->position;
    }

    friend class Skeleton_simplex_range;
   private:
    Bitmap_cubical_complex<T>* b;
    size_t position;
    unsigned dimension;
  };

  /**
   * Class needed for compatibility with Gudhi. Not useful for other purposes.
   **/
  class Skeleton_simplex_range {
    // Range over the simplices of the complex in the order of the filtration.
    // .begin() and .end() return type Filtration_simplex_iterator.
   public:
    typedef Skeleton_simplex_iterator const_iterator;
    typedef Skeleton_simplex_iterator iterator;

    Skeleton_simplex_range(Bitmap_cubical_complex<T>* b, unsigned dimension) : b(b), dimension(dimension) { }

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
    unsigned dimension;
  };

  /**
   * Function needed for compatibility with Gudhi. Not useful for other purposes.
   **/
  Skeleton_simplex_range skeleton_simplex_range(unsigned dimension) {
    if (globalDbg) {
      std::cerr << "Skeleton_simplex_range skeleton_simplex_range( unsigned dimension )\n";
    }
    return Skeleton_simplex_range(this, dimension);
  }

  friend class is_before_in_filtration<T>;

 protected:
  std::vector< size_t > key_associated_to_simplex;
  std::vector< size_t > simplex_associated_to_key;
};  // Bitmap_cubical_complex

template <typename T>
void Bitmap_cubical_complex<T>::initialize_simplex_associated_to_key() {
  if (globalDbg) {
    std::cerr << "void Bitmap_cubical_complex<T>::initialize_elements_ordered_according_to_filtration() \n";
  }
  this->simplex_associated_to_key = std::vector<size_t>(this->data.size());
  std::iota(std::begin(simplex_associated_to_key), std::end(simplex_associated_to_key), 0);
  std::sort(simplex_associated_to_key.begin(),
            simplex_associated_to_key.end(),
            is_before_in_filtration<T>(this));

  // we still need to deal here with a key_associated_to_simplex:
  for ( size_t i = 0  ; i != simplex_associated_to_key.size() ; ++i )
  {
    this->key_associated_to_simplex[ simplex_associated_to_key[i] ] = i;
  }
}

template <typename T>
class is_before_in_filtration {
 public:
  explicit is_before_in_filtration(Bitmap_cubical_complex<T> * CC)
      : CC_(CC) { }

  bool operator()(const typename Bitmap_cubical_complex<T>::Simplex_handle& sh1,
                  const typename Bitmap_cubical_complex<T>::Simplex_handle& sh2) const {
    // Not using st_->filtration(sh1) because it uselessly tests for null_simplex.
    typename T::filtration_type fil1 = CC_->data[sh1];
    typename T::filtration_type fil2 = CC_->data[sh2];
    if (fil1 != fil2) {
      return fil1 < fil2;
    }
    // in this case they are on the same filtration level, so the dimension decide.
    size_t dim1 = CC_->get_dimension_of_a_cell(sh1);
    size_t dim2 = CC_->get_dimension_of_a_cell(sh2);
    if (dim1 != dim2) {
      return dim1 < dim2;
    }
    // in this case both filtration and dimensions of the considered cubes are the same. To have stable sort, we simply
    // compare their positions in the bitmap:
    return sh1 < sh2;
  }

 protected:
  Bitmap_cubical_complex<T>* CC_;
};

}  // namespace Cubical_complex

}  // namespace Gudhi

#endif  // BITMAP_CUBICAL_COMPLEX_H_
