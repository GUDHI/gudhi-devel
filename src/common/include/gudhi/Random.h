/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef GUDHI_RANDOM_H_
#define GUDHI_RANDOM_H_

#include <random>  // for std::random_device, std::mt19937_64, std::uniform_real_distribution
#include <vector>
#include <algorithm>  // for std::generate
#include <cstddef>  // for std::size_t
#ifdef DEBUG_TRACES
#include <iostream>
#endif  // DEBUG_TRACES

namespace Gudhi {

namespace random {
  
  /**
   * @class Random
   * @brief A generic random number generator class using C++11 random facilities.
   *
   * This class provides a flexible interface for generating random numbers of various types
   * (integral and floating-point) and ranges, using a customizable random engine.
   * If you do not need a specific one, consider using the one provided by `Gudhi::random::get_default_random()`.
   *
   * @tparam Engine The random engine type to use (default: `std::mt19937_64`).
   */
  template<typename Engine = std::mt19937_64>
  class Random {
   public:
    /**
     * @brief Constructs a `Random` object with an optional seed. If no seed is provided, a random seed is generated.
     *
     * @param seed (Optional) The seed value for the random engine.
     */
    explicit Random(typename Engine::result_type seed = std::random_device{}())
      : gen_(seed) {
#ifdef DEBUG_TRACES
        std::clog << "Random ctor " << this << " - " << seed << "\n";
#endif  // DEBUG_TRACES
    }

    /**
     * @brief Sets the seed for the random engine.
     *
     * @param seed The new seed value.
     */
    void set_seed(typename Engine::result_type seed) {
      gen_.seed(seed);
    }
    
    /**
     * @brief Generates a random number in the range [min, max].
     *
     * @tparam Type The type of the random number (must be integral or floating-point).
     * @param min The minimum value (inclusive).
     * @param max The maximum value (inclusive).
     * @return A random number of type `Type` in the range [min, max].
     */
    template <typename Type>
    Type get(const Type& min, const Type& max) {
      if constexpr (std::is_floating_point_v<Type>) {
        std::uniform_real_distribution<Type> dis(min, max);
        return dis(gen_);
      } else if constexpr (std::is_integral_v<Type>) {
        std::uniform_int_distribution<Type> dis(min, max);
        return dis(gen_);
      }
    }

    /**
     * @brief Generates a random number in the range [0, max].
     *
     * @tparam Type The type of the random number (must be integral or floating-point).
     * @param max The maximum value (inclusive).
     * @return A random number of type `Type` in the range [0, max].
     */
    template <typename Type>
    Type get(const Type& max) {
      return get<Type>(static_cast<Type>(0), max);
    }

    /**
     * @brief Generates a random number in the range [0, 1].
     *
     * @tparam Type The type of the random number (must be integral or floating-point).
     * @return A random number of type `Type` in the range [0, 1].
     */
    template <typename Type>
    Type get() {
      return get<Type>(static_cast<Type>(0), static_cast<Type>(1));
    }

    /**
     * @brief Generates a vector of N random numbers in the range [min, max].
     *
     * @tparam Type The type of the random numbers (must be integral or floating-point).
     * @param nbr The number of random numbers to generate.
     * @param min The minimum value (inclusive).
     * @param max The maximum value (inclusive).
     * @return A vector of random numbers of type `Type` in the range [min, max].
     */
    template <typename Type>
    std::vector<Type> get_range(std::size_t nbr, const Type& min, const Type& max) {
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
    
    /**
     * @brief Generates a vector of N random numbers in the range [0, max].
     *
     * @tparam Type The type of the random numbers (must be integral or floating-point).
     * @param nbr The number of random numbers to generate.
     * @param max The maximum value (inclusive).
     * @return A vector of random numbers of type `Type` in the range [0, max].
     */
    template <typename Type>
    std::vector<Type> get_range(std::size_t nbr, const Type& max) {
      return get_range<Type>(nbr, static_cast<Type>(0), max);
    }
    
    /**
     * @brief Generates a N vector of random numbers in the range [0, 1].
     *
     * @tparam Type The type of the random numbers (must be integral or floating-point).
     * @param nbr The number of random numbers to generate.
     * @return A vector of random numbers of type `Type` in the range [0, 1].
     */
    template <typename Type>
    std::vector<Type> get_range(std::size_t nbr) {
      return get_range<Type>(nbr, static_cast<Type>(0), static_cast<Type>(1));
    }
    
    /**
     * @return The underlying random engine.
     */
    Engine& get_engine() {
      return gen_;
    }

   private:
    Engine gen_;
  };
  
  /**
   * @brief Returns a thread-local default random number generator instance.
   *
   * @return A reference to the thread-local default @ref Random instance.
   */
  inline Random<>& get_default_random() {
    static thread_local Random<> default_random;
#ifdef DEBUG_TRACES
    std::clog << "get_default_random() " << &default_random << "\n";
#endif  // DEBUG_TRACES
    return default_random;
  }

}  // namespace random

}  // namespace Gudhi

#endif  // GUDHI_RANDOM_H_
