/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2018 Inria
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


#include <tbb/task_scheduler_init.h>  // for tbb::task_scheduler_init
#include <tbb/pipeline.h>  // for tbb::pipeline

#include <random>  // for std::default_random_engine, uniform_int_distribution
#include <chrono>  // for std::chrono::seconds
#include <thread>  // for std::this_thread::sleep_for
#include <iostream>  // for std::cout
#include <numeric>  // for std::accumulate
#include <utility>  // for std::pair
#include <mutex>  // for std::mutex
#include <cassert>  // for assert
#include <vector>  // for std::vector
#include <cstddef>  // for std::size_t

class Random_parallel_filter: public tbb::filter {
public:
  Random_parallel_filter(int number) :
      filter(parallel) {
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0,9);

    for (int i = 0; i < number; i++) {
      random_values_.push_back(static_cast<int>(distribution(generator)));
      std::cout << "Insert random_values " << random_values_[i] << std::endl;
    }
    random_values_iter_ = random_values_.begin();
    accumulate_ = std::accumulate(random_values_.begin(), random_values_.end(), 0);
    std::cout << "Result will be  " << accumulate_ << std::endl;
  }

  int get_accumulate() const {
    return accumulate_;
  }

private:
  std::vector<int> random_values_;
  std::vector<int>::iterator random_values_iter_;
  std::mutex iterator_mutex;
  int accumulate_;

  void *operator()(void *) {
    if (random_values_iter_ >= random_values_.end())
      // End condition
      return NULL;

    // Save what is required before incrementing for next iteration
    // Dereference vector iterator
    iterator_mutex.lock();
    int random_time = *random_values_iter_;
    std::size_t index = random_values_iter_ - random_values_.begin();
    std::cout << "index = " << index << " - random_time = " << random_time << std::endl;
    // For next iteration - to be done as soon as possible
    ++random_values_iter_;
    iterator_mutex.unlock();

    std::this_thread::sleep_for(std::chrono::seconds(random_time));

    std::pair<std::size_t, int>* index_and_value_ptr = new std::pair<std::size_t, int>(index, random_time);
    return index_and_value_ptr;
  }
};

class Accumulate_filter: public tbb::filter {
public:
  Accumulate_filter(int number) :
      filter(serial_in_order),
      index_(0),
      accumulate_(0),
      parallel_values_(number, -1) {
    std::cout << "Accumulate_filter creation" << std::endl;
  }

  int get_accumulate() const {
    return accumulate_;
  }

private:
  void *operator()(void * item ) {
    // Cast as sender format
    std::pair<std::size_t, int>* index_and_value = static_cast<std::pair<std::size_t, int>*>(item);

    std::cout << "index " << index_and_value->first << " - value " << index_and_value->second << std::endl;
    parallel_values_[index_and_value->first] = index_and_value->second;

    if (index_ == index_and_value->first) {
      index_mutex.lock();
      while ((parallel_values_[index_] != -1) && index_ < parallel_values_.size()) {
        accumulate_ += parallel_values_[index_];
        std::cout << "accumulate index= " << index_ << " - value= " << parallel_values_[index_] << std::endl;
        ++index_;
      }
      index_mutex.unlock();
    }

    std::cout << "index " << index_and_value->first << " - value " << index_and_value->second << " - OUT" << std::endl;
    // Was new'ed by Random_parallel_filter
    delete index_and_value;
    return NULL;
  }

  std::size_t index_;
  int accumulate_;
  std::vector<int> parallel_values_;
  std::mutex index_mutex;

};

int main(int argc, char * argv[]) {

  tbb::task_scheduler_init init_parallel;

  // Create the pipeline
  tbb::pipeline pipeline;

  int number = atoi(argv[1]);
  // Create file-reading writing stage and add it to the pipeline
  Random_parallel_filter random_parallel(number);
  pipeline.add_filter( random_parallel );

  Accumulate_filter accumulate(number);
  pipeline.add_filter( accumulate );

  // Run the pipeline
  // Need more than one token in flight per thread to keep all threads
  // busy; 2-4 works
  pipeline.run( init_parallel.default_num_threads() * 4 );
  // Sequential
  //pipeline.run(1);

  assert(accumulate.get_accumulate() == random_parallel.get_accumulate());
  std::cout << "Accumulate " << accumulate.get_accumulate() << " - versus " << random_parallel.get_accumulate() << std::endl;

  return 0;
}
