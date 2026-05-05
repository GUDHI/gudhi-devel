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
   * @defgroup reproducibility Reproducibility
   * 
   * Some of the GUDHI functionnalities are using randomness. Here is a list of classes and methods that are using
   * internally @ref Gudhi::random::Random
   * 
   * - @ref Gudhi::subsampling::choose_n_farthest_points
   * - @ref Gudhi::subsampling::pick_n_random_points
   * - @ref Gudhi::Persistence_representations::Sliced_Wasserstein
   * - @ref Gudhi::cover_complex::Cover_complex
   * 
   * In order to reproduce the results, a specific GUDHI random generator is available, and one can set its seed:
   * 
   * @code{.cpp}
   * #include <gudhi/Random.h>
   * 
   * // Set the seed for the default generator that is used internally in GUDHI
   * Gudhi::random::set_seed(42);
   * @endcode
   * 
   * Other GUDHI functionnalities, that are based on CGAL, are using also randomness, but depends on
   * <a href="https://doc.cgal.org/latest/Generator/classCGAL_1_1Random.html">CGAL::Random</a>.
   * Here is the list:
   * 
   * - Gudhi::random::generate_points_on_plane 
   * - Gudhi::random::generate_points_on_moment_curve
   * - Gudhi::random::generate_points_on_torus_3D
   * - Gudhi::random::generate_grid_points_on_torus_d
   * - Gudhi::random::generate_points_on_torus_d
   * - Gudhi::random::generate_points_on_sphere_d
   * - Gudhi::random::generate_points_in_ball_d
   * - Gudhi::random::generate_points_in_cube_d
   * - Gudhi::random::generate_points_on_two_spheres_d
   * - Gudhi::random::generate_points_on_3sphere_and_circle
   * - Gudhi::random::generate_points_on_klein_bottle_3D
   * - Gudhi::random::generate_points_on_klein_bottle_4D
   * - Gudhi::random::generate_points_on_klein_bottle_variant_5D
   * - Gudhi::coxeter_triangulation::random_orthogonal_matrix
   * 
   * In order to reproduce the results, one can set the seed of the default CGAL random generator, that is used
   * internally in GUDHI, by doing:
   * 
   * @code{.cpp}
   * #include <CGAL/Random.h>
   * 
   * // Get the default static thread local CGAL random generator and set its seed
   * // It is the one that is used internally
   * CGAL::get_default_random() = CGAL::Random(42);
   * @endcode
   *
   * @{
   */

   /**
    * The default random number generator type returned by @ref get_default_random
    */
  using Random_generator = std::mt19937_64;
  
  /**
   * @brief Returns a thread-local default random number generator instance.
   * 
   * @return A reference to the thread-local default random instance.
   */
  inline Random_generator& get_default_random() {
    static thread_local Random_generator default_random{std::random_device{}()};
#ifdef DEBUG_TRACES
    std::clog << "get_default_random() " << &default_random << "\n";
#endif  // DEBUG_TRACES
    return default_random;
  }

  /**
   * @brief Sets the seed for the default random engine.
   *
   * @param seed The new seed value.
   */
  void set_seed(Random_generator::result_type seed) {
    get_default_random().seed(seed);
  }

  /**
   * @brief Generates a random number in the range [min, max].
   *
   * @tparam Type The type of the random number (must be integral or floating-point).
   * @tparam CustomRandomGenerator The type of the random generator. Default is @ref Random_generator.
   * @param min The minimum value (inclusive).
   * @param max The maximum value (inclusive).
   * @param rng A random generator. Use the default one if not set.
   * @return A random number of type `Type` in the range [min, max].
   */
   template <typename Type, typename CustomRandomGenerator = Random_generator&>
   Type get_uniform(const Type& min, const Type& max, CustomRandomGenerator&& rng = get_default_random()) {
       if constexpr (std::is_floating_point_v<Type>) {
           std::uniform_real_distribution<Type> dis(min, max);
           return dis(rng);
       } else if constexpr (std::is_integral_v<Type>) {
           std::uniform_int_distribution<Type> dis(min, max);
           return dis(rng);
       }
   }

  /**
   * @brief Generates a vector of N random numbers in the range [min, max].
   *
   * @tparam Type The type of the random numbers (must be integral or floating-point).
   * @tparam CustomRandomGenerator The type of the random generator. Default is @ref Random_generator.
   * @param nbr The number of random numbers to generate.
   * @param min The minimum value (inclusive).
   * @param max The maximum value (inclusive).
   * @param rng A random generator. Use the default one if not set.
   * @return A vector of random numbers of type `Type` in the range [min, max].
   */
  template <typename Type, typename CustomRandomGenerator = Random_generator&>
  std::vector<Type> get_uniform_range(std::size_t nbr, const Type& min, const Type& max,
                              CustomRandomGenerator&& rng = get_default_random()) {
    // With C++26, this should be replaced with std::ranges::generate_random
    std::vector<Type> result(nbr);
    if constexpr (std::is_floating_point_v<Type>) {
      std::uniform_real_distribution<Type> dis(min, max);
      std::generate(result.begin(), result.end(), [&]() { return dis(rng); });
    } else if constexpr (std::is_integral_v<Type>) {
      std::uniform_int_distribution<Type> dis(min, max);
      std::generate(result.begin(), result.end(), [&]() { return dis(rng); });
    }
    return result;
  }
  
  
  /** @} */  // end defgroup random

}  // namespace random

}  // namespace Gudhi

#endif  // GUDHI_RANDOM_H_
