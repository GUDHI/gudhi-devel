/* This file is part of the Gudhi Library. The Gudhi library 
 *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
 *    library for computational topology.
 *
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014  INRIA Sophia Antipolis-Mediterranee (France)
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
 * 
 */

#ifndef UTILS_MCLOCK_H_
#define UTILS_MCLOCK_H_

#include <sys/time.h>

class MClock {
 public:
  MClock() {
    end_called = false;
    begin();
  }

  void begin() const {
    end_called = false;
    gettimeofday(&startTime, NULL);
  }

  void end() const {
    end_called = true;
    gettimeofday(&endTime, NULL);
  }

  friend std::ostream& operator<<(std::ostream& stream, const MClock& clock) {
    // if(!clock.end_called) clock.end();
    if (!clock.end_called) {
      stream << "end not called";
    } else {
      long totalTime = (clock.endTime.tv_sec - clock.startTime.tv_sec) * 1000000L;
      totalTime += (clock.endTime.tv_usec - clock.startTime.tv_usec);
      stream << clock.num_seconds() << "s";
    }
    return stream;
  }

  double num_seconds() const {
    if (!end_called) end();
    long totalTime = (endTime.tv_sec - startTime.tv_sec) * 1000000L;
    totalTime += (endTime.tv_usec - startTime.tv_usec);
    return (totalTime / 1000L) / static_cast<double>(1000);
  }

 private:
  mutable struct timeval startTime, endTime;
  mutable bool end_called;
};

#endif  // UTILS_MCLOCK_H_
