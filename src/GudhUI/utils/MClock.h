/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
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
