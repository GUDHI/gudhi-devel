/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
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
