/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2017  Swansea University, UK
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

#ifndef WRITING_PERSISTENCE_TO_FILE_H
#define WRITING_PERSISTENCE_TO_FILE_H

#include <iostream>
#include <string>
#include <limits>

namespace Gudhi {

/**
* This is a class to store persistence intervals. Its main purpose is to
* exchange data in between different packages and provide unified way
* of writing a collection of persistence intervals to file.
**/
template <typename Filtration_type, typename Coefficient_field>
class Persistence_interval_common {
 public:
  /**
   * Constructor taking as an input birth and death of the pair.
  **/
  Persistence_interval_common(Filtration_type birth, Filtration_type death)
      : birth_(birth),
        death_(death),
        dimension_(std::numeric_limits<unsigned>::max),
        arith_element_(std::numeric_limits<Coefficient_field>::max()) {}

  /**
   * Constructor taking as an input birth, death and dimension of the pair.
  **/
  Persistence_interval_common(Filtration_type birth, Filtration_type death, unsigned dim)
      : birth_(birth), death_(death), dimension_(dim), arith_element_(std::numeric_limits<Coefficient_field>::max()) {}

  /**
* Constructor taking as an input birth, death, dimension of the pair as well
* as the number p such that this interval is present over Z_p field.
**/
  Persistence_interval_common(Filtration_type birth, Filtration_type death, unsigned dim, Coefficient_field field)
      : birth_(birth), death_(death), dimension_(dim), arith_element_(field) {}

  /**
   * Operator to compare two persistence pairs. During the comparision all the
   * fields: birth, death, dimensiona and arith_element_ are taken into account
   * and they all have to be equal for two pairs to be equal.
  **/
  inline bool operator==(const Persistence_interval_common& i2) {
    return ((this->birth_ == i2.birth_) && (this->death_ == i2.death_) && (this->dimension_ == i2.dimension_) &&
            (this->arith_element_ == i2.arith_element_));
  }

  /**
   * Check if two persistence paris are not equal.
  **/
  inline bool operator!=(const Persistence_interval_common& i2) { return (!((*this) == i2)); }

  /**
   * Operator to compare objects of a type Persistence_interval_common.
   * One intervals is smaller than the other if it has lower persistence.
   * Note that this operator do not take Arith_element into account when doing comparisions.
  **/
  inline bool operator<(const Persistence_interval_common& i2) {
    return fabs(this->death_ - this->birth_) < fabs(i2.death_ - i2.birth_);
  }

  friend std::ostream& operator<<(std::ostream& out, const Persistence_interval_common& it) {
    if (it.arith_element_ != std::numeric_limits<Coefficient_field>::max()) {
      out << it.arith_element_ << "  ";
    }
    if (it.dimension_ != std::numeric_limits<unsigned>::max()) {
      out << it.dimension_ << " ";
    }
    out << it.birth_ << " " << it.death_ << " ";
    return out;
  }

 private:
  Filtration_type birth_;
  Filtration_type death_;
  unsigned dimension_;
  Coefficient_field arith_element_;
};

/**
 * This function write a vector<Persistence_interval_common> to a stream
**/
template <typename Persistence_interval_range>
void write_persistence_intervals_to_stream(const Persistence_interval_range& intervals,
                                           std::ostream& out = std::cout) {
  for (auto interval : intervals) {
    out << interval << "\n";
  }
}
}

#endif  // WRITING_PERSISTENCE_TO_FILE_H
