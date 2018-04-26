/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2015 Inria
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

#ifndef BITMAP_CUBICAL_COMPLEX_COUNTER_H_
#define BITMAP_CUBICAL_COMPLEX_COUNTER_H_

#include <iostream>
#include <vector>
#include <cstddef>

namespace Gudhi {

namespace cubical_complex {

/**
 * @brief This is an implementation of a counter being a vector of integers. 
 * @details The constructor of the class takes as an input two vectors W and V. 
 * It assumes that W < V coordinatewise.
 * If the initial counter W is not specified, it is assumed to be vector of zeros. 
 * The class allows to iterate between W and V by using increment() function.
 * The increment() function returns a bool value. 
 * The current counter reach the end counter V if the value returned by the increment function is FALSE.
 * This class is needed for the implementation of a bitmapCubicalComplex.
 **/
class counter {
 public:
  /**
   * Constructor of a counter class. It takes only the parameter which is the end value of the counter. 
   * The default beginning value is a vector of the same length as the endd, filled-in with zeros.
   **/
  counter(const std::vector<unsigned>& endd) : begin(endd.size(), 0), end(endd), current(endd.size(), 0) { }

  /**
   * Constructor of a counter class. It takes as the input beginn and end vector. 
   * It assumes that begin vector is lexicographically below the end vector.
   **/
  counter(const std::vector< unsigned >& beginn, const std::vector< unsigned >& endd) : begin(beginn), end(endd), current(endd.size(), 0) {
    if (beginn.size() != endd.size())
      throw "In constructor of a counter, begin and end vectors do not have the same size. Program terminate";
  }

  /**
   * Function to increment the counter. If the value returned by the function is true, 
   * then the incrementation process was successful.
   * If the value of the function is false, that means, that the counter have reached its end-value.
   **/
  bool increment() {
    std::size_t i = 0;
    while ((i != this->end.size()) && (this->current[i] == this->end[i])) {
      ++i;
    }

    if (i == this->end.size())return false;
    ++this->current[i];
    for (std::size_t j = 0; j != i; ++j) {
      this->current[j] = this->begin[j];
    }
    return true;
  }

  /**
   * Function to check if we are at the end of counter.
   **/
  bool isFinal() {
    for (std::size_t i = 0; i != this->current.size(); ++i) {
      if (this->current[i] == this->end[i])return true;
    }
    return false;
  }

  /**
   * Function required in the implementation of bitmapCubicalComplexWPeriodicBoundaryCondition. 
   * Its aim is to find an counter corresponding to the element the following
   * boundary element is identified with when periodic boundary conditions are imposed.
   **/
  std::vector< unsigned > find_opposite(const std::vector< bool >& directionsForPeriodicBCond) {
    std::vector< unsigned > result;
    for (std::size_t i = 0; i != this->current.size(); ++i) {
      if ((this->current[i] == this->end[i]) && (directionsForPeriodicBCond[i] == true)) {
        result.push_back(this->begin[i]);
      } else {
        result.push_back(this->current[i]);
      }
    }
    return result;
  }

  /**
   * Function checking at which positions the current value of a counter is the final value of the counter.
   **/
  std::vector< bool > directions_of_finals() {
    std::vector< bool > result;
    for (std::size_t i = 0; i != this->current.size(); ++i) {
      if (this->current[i] == this->end[i]) {
        result.push_back(true);
      } else {
        result.push_back(false);
      }
    }
    return result;
  }

  /**
   * Function to write counter to the stream.
   **/
  friend std::ostream& operator<<(std::ostream& out, const counter& c) {
    // std::cerr << "c.current.size() : " << c.current.size() << endl;
    for (std::size_t i = 0; i != c.current.size(); ++i) {
      out << c.current[i] << " ";
    }
    return out;
  }

 private:
  std::vector< unsigned > begin;
  std::vector< unsigned > end;
  std::vector< unsigned > current;
};

}  // namespace cubical_complex

namespace Cubical_complex = cubical_complex;

}  // namespace Gudhi

#endif  // BITMAP_CUBICAL_COMPLEX_COUNTER_H_
