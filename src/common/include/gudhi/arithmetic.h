/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef GUDHI_ARITHMETIC_H_
#define GUDHI_ARITHMETIC_H_

namespace Gudhi {

// This method is templated but its aim is to test small primes (not the largest prime that fits in the number type)
template<typename Integer>
static constexpr bool is_prime(const Integer p) {
  if (p <= 1) return false;
  if (p <= 3) return true;
  if (p % 2 == 0 || p % 3 == 0) return false;
  
  for (long i = 5; i * i <= p; i = i + 6)
    if (p % i == 0 || p % (i + 2) == 0) return false;
  
  return true;
}

}  // namespace Gudhi

#endif   // GUDHI_ARITHMETIC_H_
