/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <iostream>
#include <random>
#include <numeric>  // for std::iota
#include <algorithm>  // for std::shuffle
#include <vector>

std::random_device rd;

template <class Vertex_handle>
std::vector<Vertex_handle> rand_int_range(int subset_min_size,
                                          int subset_max_size,
                                          Vertex_handle range_max_value)
{
#ifdef DEBUG_TRACES
  std::clog << "rand_int_range - subset_min_size = " << subset_min_size << " - subset_max_size = " << subset_max_size
            << " - range_max_value = " << range_max_value << std::endl;
#endif  // DEBUG_TRACES

  std::vector<Vertex_handle> range(range_max_value);
  std::iota(range.begin(), range.end(), 0); // range is { 0, 1, 2, ..., 99 } when range_max_value is 100

  std::shuffle(range.begin(), range.end(), std::mt19937 { rd() });

  std::uniform_int_distribution<int> dist(subset_min_size, subset_max_size);
  // Return a subset, which size is in between [subset_min_size; subset_max_size], of shuffled range
  // {1, 9, 13, 19, 86, 36}, for example
  return std::vector<Vertex_handle> {range.begin(), range.begin() + dist(rd)};
}
