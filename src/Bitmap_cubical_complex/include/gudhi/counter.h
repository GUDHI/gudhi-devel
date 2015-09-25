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

#ifndef COUNTER_H_
#define COUNTER_H_

#include <iostream>
#include <vector>

using namespace std;

/**
 * This is an implementation of a simple counter. It is needed for the implementation of a bitmapCubicalComplex.
 **/

class counter {
 public:

  /**
   * Constructor of a counter class. It takes only the parameter which is the end value of the counter. The default beginning value is a vector of the same length as the endd, filled-in with zeros.
   **/
  counter(std::vector< int > endd) {
    for (size_t i = 0; i != endd.size(); ++i) {
      this->current.push_back(0);
      this->begin.push_back(0);
      this->end.push_back(endd[i]);
    }
  }

  /**
   * Constructor of a counter class. It takes as the input beginn and end vector. It assumes that begin vector is lexicographically below the end vector.
   **/
  counter(std::vector< int > beginn, std::vector< int > endd) {
    if (beginn.size() != endd.size())throw "In constructor of a counter, begin and end vectors do not have the same size. Program terminate";
    for (size_t i = 0; i != endd.size(); ++i) {
      this->current.push_back(0);
      this->begin.push_back(0);
      this->end.push_back(endd[i]);
    }
  }

  /**
   * Function to increment the counter. If the value returned by the function is true, then the incrementation process was successful.
   * If the value of the function is false, that means, that the counter have reached its end-value.
   **/
  bool increment() {
    size_t i = 0;
    while ((i != this->end.size()) && (this->current[i] == this->end[i])) {
      ++i;
    }

    if (i == this->end.size())return false;
    ++this->current[i];
    for (size_t j = 0; j != i; ++j) {
      this->current[j] = this->begin[j];
    }
    return true;
  }

  /**
   * Function to check if we are at the end of counter.
   **/
  bool isFinal() {
    for (size_t i = 0; i != this->current.size(); ++i) {
      if (this->current[i] == this->end[i])return true;
    }
    return false;
  }

  /**
   * Function required in the implementation of bitmapCubicalComplexWPeriodicBoundaryCondition. Its aim is to find an counter corresponding to the element the following
   * boundary element is identified with when periodic boundary conditions are imposed.
   **/
  std::vector< int > findOpposite(std::vector< bool > directionsForPeriodicBCond) {
    std::vector< int > result;
    for (size_t i = 0; i != this->current.size(); ++i) {
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
  std::vector< bool > directionsOfFinals() {
    std::vector< bool > result;
    for (size_t i = 0; i != this->current.size(); ++i) {
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
    //cerr << "c.current.size() : " << c.current.size() << endl;
    for (size_t i = 0; i != c.current.size(); ++i) {
      out << c.current[i] << " ";
    }
    return out;
  }
 private:
  std::vector< int > begin;
  std::vector< int > end;
  std::vector< int > current;
};

#endif  // COUNTER_H_
