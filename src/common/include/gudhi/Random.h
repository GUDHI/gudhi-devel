/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef GUDHI_RANDOM_H_
#define GUDHI_RANDOM_H_

#include <random>  // for std::random_device, std::mt19937, std::uniform_real_distribution
#include <vector>
#include <algorithm>  // for std::generate
#include <cstddef>  // for std::size_t
#ifdef DEBUG_TRACES
#include <iostream>
#endif  // DEBUG_TRACES

namespace Gudhi {

namespace random {
  
  template<typename Engine = std::mt19937_64>
  class Random {
   public:
    explicit Random(typename Engine::result_type seed = std::random_device{}())
      : gen_(seed) {
#ifdef DEBUG_TRACES
        std::clog << "Random ctor " << this << " - " << seed << "\n";
#endif  // DEBUG_TRACES
    }

    void set_seed(typename Engine::result_type seed) {
      gen_.seed(seed);
    }
    
    template <typename Type>
    Type operator() (const Type& min, const Type& max) {
      if constexpr (std::is_floating_point_v<Type>) {
        std::uniform_real_distribution<Type> dis(min, max);
        return dis(gen_);
      } else if constexpr (std::is_integral_v<Type>) {
        std::uniform_int_distribution<Type> dis(min, max);
        return dis(gen_);
      }
    }

    template <typename Type>
    Type operator() (const Type& max) {
      return operator()<Type>(static_cast<Type>(0), max);
    }

    template <typename Type>
    Type operator() () {
      return operator()<Type>(static_cast<Type>(0), static_cast<Type>(1));
    }

    Engine get_engine() const {
      return gen_;
    }

   private:
    Engine gen_;
  };
  
  inline Random<>& get_default_random()
  {
    static thread_local Random<> default_random;
#ifdef DEBUG_TRACES
    std::clog << "get_default_random() " << &default_random << "\n";
#endif  // DEBUG_TRACES
    return default_random;
  }

}  // namespace random

}  // namespace Gudhi

#endif  // GUDHI_RANDOM_H_
