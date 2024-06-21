/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef RANDOM__H_
#define RANDOM__H_

#include <random>  // for std::random_device, std::mt19937, std::uniform_real_distribution
#include <vector>
#include <algorithm>  // for std::generate
#include <cstddef>  // for std::size_t

namespace Gudhi {
  std::random_device rd;

  class Random {
   public:
    Random() : gen_(rd()) {}
    Random(std::uint_fast32_t seed) : gen_(seed) {}

    template <typename Type>
    Type get(const Type& min = 0, const Type& max = 1) {
      if constexpr (std::is_floating_point_v<Type>) {
        std::uniform_real_distribution<Type> dis(min, max);
        return dis(gen_);
      } else if constexpr (std::is_integral_v<Type>) {
        std::uniform_int_distribution<Type> dis(min, max);
        return dis(gen_);
      }
    }

    template <typename Type>
    std::vector<Type> get_range(std::size_t nbr, const Type& min = 0, const Type& max = 1) {
      std::vector<Type> result(nbr);
      if constexpr (std::is_floating_point_v<Type>) {
        std::uniform_real_distribution<Type> dis(min, max);
        std::generate(result.begin(), result.end(), [&]() { return dis(gen_); });
      } else if constexpr (std::is_integral_v<Type>) {
        std::uniform_int_distribution<Type> dis(min, max);
        std::generate(result.begin(), result.end(), [&]() { return dis(gen_); });
      }
      return result;
    }

   private:
    std::mt19937 gen_;
  };
}  // namespace Gudhi

#endif  // RANDOM__H_