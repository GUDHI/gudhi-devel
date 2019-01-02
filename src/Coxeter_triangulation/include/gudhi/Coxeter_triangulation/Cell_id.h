#ifndef CELL_ID_H_
#define CELL_ID_H_

#include <vector>
#include <iostream>

namespace Gudhi {

  /** \private
   * \class Cell_id
   * \brief A data structure to represent the cell coordinates of a given cell in a Coxeter triangulation.
   * \details Each coordinate in a tuple of cell coordinates has an integer and a boolean value.
   The integer component indicates the position, while the boolean value indicates whether the cell lies in a slab (if the value is false) or on a hyperplane (if the value is true).
   The coordinates are arranged in a random-access stack.
   * \ingroup coxeter_triangulation
   */
  struct Cell_id {

    typedef long int value_type; // TODO: Maybe put it as a template

    /** \brief A const iterator over the integer values of the tuple of cell coordinates 
     */
    typedef typename std::vector<value_type>::const_iterator const_iterator;
    typedef typename std::vector<bool>::const_iterator const_mask_iterator;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /* @name Constructor
   */
    
    //@{

    /** \brief Constructs an empty tuple of cell coordinates.
	\details The level and the dimension are assigned to the default values (0).
     */
    Cell_id()
      : level_(0), dimension_(0) {}

    /** \brief Constructs an empty tuple of cell coordinates of a given level.
	\details The dimension are assigned to the default value (0).
     */
    Cell_id(double level)
      : level_(level), dimension_(0) {}

    /** \brief Constructs an empty tuple of cell coordinates of given level and dimension.
     */
    Cell_id(double level, unsigned dimension)
      : level_(level), dimension_(dimension) {}

    //@}

    
    /** \brief Returns the integer value of a coordinate.
	\details The function returns the integer value of a coordinate by its index.
	@param[out] The integer value of a coordinate.
	@param[in] index The index of the coordinate.
     */
    value_type value(std::size_t index) const {
      return value_[index];
    }

    /** \brief Returns the boolean value of a coordinate.
	\details The function returns the boolean value of a coordinate by its index.
	@param[out] The boolean value of a coordinate.
	@param[in] index The index of the coordinate.
     */
    bool mask(std::size_t index) const {
      return mask_[index];
    }

    /** \brief Returns the level of the tuple of cell coordinates.
	@param[out] The level of the tuple of cell coordinates.
     */
    double level() const {
      return level_;
    }
  
    /** \brief Returns the dimension of the tuple of cell coordinates.
	@param[out] The dimension of the tuple of cell coordinates.
     */
    double dimension() const {
      return dimension_;
    }

    /** \brief Returns the size of the tuple of cell coordinates.
	@param[out] The size of the tuple.
    */
    std::size_t size() const {
      return value_.size();
    }

    /** \brief Checks if the tuple of cell coordinates is empty.
	@param[out] A boolean value, which is true if and only if the tuple is empty.
    */
    bool empty() const {
      return value_.size() == 0;
    }

    
    /** \brief Pushes a coordinate with the given integer and boolean values to the back of the tuple
	\details The function returns the integer value of a coordinate by its index.
	@param[in] value The integer value of the input coordinate.
	@param[in] mask  The boolean value of the input coordinate.
     */    
    void push_back(value_type value, bool mask = false) {
      value_.push_back(value);
      mask_.push_back(mask);
    }

    /** \brief Pops a coordinate from the back of the tuple
     */    
    void pop_back() {
      value_.pop_back();
      mask_.pop_back();
    }

    /** \brief Sets the level to a given value.
	@param[in] new_level The new level of the tuple of cell coordinates.
     */
    void set_level(double new_level) {
      level_ = new_level;
    }
    
    /** \brief Sets the dimension to a given value.
	@param[in] new_dim The new dimension of the tuple of cell coordinates.
     */
    void set_dimension(std::size_t new_dim) {
      dimension_ = new_dim;
    }
    
    /** \brief Reserves storage.
	\detail Allocates additional storage if the input parameter new_cap is greater than the current capacity. Otherwise, does nothing.
	Invalidates all iterators if the parameter new_cap is greater than the current capacity.
	@param[in] new_cap New capacity of the tuple of cell coordinates.
     */  
    void reserve(std::size_t new_cap) {
      value_.reserve(new_cap);
      mask_.reserve(new_cap);
    }
  
    /** \brief Resizes the tuple to a given size.
	\details If the parameter new_size is smaller than the current size, the tuple is reduced to its first new_size coordinates.
	If the parameter new_size is greater than the current size, the additional coordinates are initialized with values (0, false).
	@param[in] new_size The new size of the tuple of cell coordinates.
     */  
    void resize(std::size_t new_size) {
      value_.resize(new_size);
      mask_.resize(new_size);
    }

    /** \brief Deletes all elements in the tuple of cell cordinates.
     */  
    void clear() {
      resize(0);
    }    

    
    /** \brief Returns the constant iterator to the integer value of the first coordinate.
	@param[out] The constant iterator to the integer value of the first coordinate.
     */
    const_iterator begin() const {
      return value_.begin();
    }

    /** \brief Returns the constant integer value iterator to the position after the last coordinate in the tuple.
	@param[out] The constant integer value iterator to the position after the last coordinate in the tuple.
     */
    const_iterator end() const {
      return value_.end();
    }

    /** \brief Returns the constant iterator to the boolean value of the first coordinate.
	@param[out] The constant iterator to the boolean value of the first coordinate.
     */
    const_mask_iterator mask_begin() const {
      return mask_.begin();
    }

    /** \brief Returns the constant boolean value iterator to the position after the last coordinate in the tuple.
	@param[out] The constant boolean value iterator to the position after the last coordinate in the tuple.
     */
    const_mask_iterator mask_end() const {
      return mask_.end();
    }

    /** \brief Returns true if all coordinates of the tuple other are equal to the coordinates of the current tuple.
     * @param[in] other A tuple of cell coordinates to compare.
     * @param[out] The result of the comparison. The result is true if the two sizes are equal, and all coordinates are equal in both integer and boolean values. Otherwise, the result is false.
     */
    bool operator==(const Cell_id& other) const {
      if (this->size() != other.size())
	return false;
      for (std::size_t k = 0; k < this->size(); ++k)
	if (this->value(k) != other.value(k) || this->mask(k) != other.mask(k))
	  return false;
      return true;
    }

    bool operator!=(const Cell_id& other) const {
      return !(*this == other);
    }

    /** \brief Check if a simplex is a face of another simplex.
     * \detail Returns true if c1 is a face of c2. The two tuples c1 and c2 are required to be valid.
     */
    static bool is_face(const Gudhi::Cell_id& c1, const Gudhi::Cell_id& c2) {
      std::size_t k = 0; 
      for (; k < c1.size(); ++k) {
	if (c2.mask(k)) {
	  if (!c1.mask(k) || c1.value(k) != c2.value(k))
	    return false;
	}
	else { // !c2.mask(k)
	  if (!c1.mask(k) && c1.value(k) != c2.value(k))
	    return false;
	  if (c1.mask(k))
	    if (c1.value(k) < c2.value(k) || c1.value(k) > c2.value(k)+1)
	      return false;
	} 
      }
      return true;
    }
    
  protected:
    double level_;
    unsigned dimension_;
    std::vector<value_type> value_;
    std::vector<bool> mask_;
  };
  
}

/** \brief The default comparison operator between two tuples of cell coordinates.
    \detail If the sizes of the two tuples are different, the smallest tuple is returned as the smallest. If the sizes are the same, then order for the comparison is the lexicographical order, with a true coordinate smaller than a false coordinate, in the case of difference.
*/
bool operator< (const Gudhi::Cell_id& lhs, const Gudhi::Cell_id& rhs) {
  if (lhs.size() < rhs.size())
    return true;
  else if (lhs.size() > rhs.size())
    return false;
  else
    for (std::size_t k = 0; k < lhs.size(); ++k)
      if (lhs.value(k) < rhs.value(k))
	return true;
      else if (lhs.value(k) < rhs.value(k))
	return false;
      else if (lhs.mask(k) && !rhs.mask(k))
	return true;
      else if (!lhs.mask(k) && rhs.mask(k))
	return false;
  return false;
}

/** \brief Inserts a tuple of cell coordinates to a stream.
 *  \details The output format of the cell coordinates is the integer values in the order between square brackets, followed by underscore and the dimension. The coordinates with true boolean values are highlighted in bold red. 
 */
std::ostream& operator<<(std::ostream& os, const Gudhi::Cell_id& c_id) {
  os << "[";
  if (c_id.empty()) {
    std::cout << "]_" << c_id.dimension();
    return os;
  }
  if (c_id.mask(0))
    os << "\033[1;31m" << c_id.value(0) << "\033[0m";
  else
    os << c_id.value(0);
  for (std::size_t i = 1; i < c_id.size(); ++i)
    if (c_id.mask(i))
      os << ", \033[1;31m" << c_id.value(i) << "\033[0m";
    else
      os << ", " << c_id.value(i);
  std::cout << "]_" << c_id.dimension();
  return os;
}


#endif
