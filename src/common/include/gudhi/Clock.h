/*    This file is part of the Gudhi Library. The Gudhi library 
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
 */

#ifndef CLOCK_H_
#define CLOCK_H_

#include <iostream>
#include <string>
#include <chrono>

namespace Gudhi {

class Clock {
 public:
  // Construct and start the timer
  Clock(const std::string& msg_ = std::string())
      : startTime(std::chrono::system_clock::now()),
      end_called(false),
      msg(msg_) { }

  // Restart the timer
  void begin() const {
    end_called = false;
    startTime = std::chrono::system_clock::now();
  }

  // Stop the timer
  void end() const {
    end_called = true;
    endTime = std::chrono::system_clock::now();
  }

  std::string message() const {
    return msg;
  }

  // Print current value to std::cout
  void print() const {
    std::cout << *this << std::endl;
  }

  friend std::ostream& operator<<(std::ostream& stream, const Clock& clock) {
    if (!clock.msg.empty())
      stream << clock.msg << ": ";

    stream << clock.num_seconds() << "s\n";
    return stream;
  }

  // Get the number of seconds between the timer start and:
  // - the last call of end() if it was called
  // - or now otherwise. In this case, the timer is not stopped.
  double num_seconds() const {
    if (!end_called) {
      auto end = std::chrono::system_clock::now();
      return std::chrono::duration_cast<std::chrono::milliseconds>(end-startTime).count() / 1000.;
    } else {
      return std::chrono::duration_cast<std::chrono::milliseconds>(endTime-startTime).count() / 1000.;
    }
  }

 private:
  mutable std::chrono::time_point<std::chrono::system_clock> startTime, endTime;
  mutable bool end_called;
  std::string msg;
};

}  // namespace Gudhi

#endif   // CLOCK_H_
