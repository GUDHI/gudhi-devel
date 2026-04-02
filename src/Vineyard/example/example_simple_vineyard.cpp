/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <iostream>
#include <limits>
#include <vector>

#include <gudhi/Rips_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/distance_functions.h>
#include <gudhi/vineyard_builder.h>
#include <gudhi/vineyard_helper.h>

using Filtration_value = double;
// short for Gudhi::vineyard::Vineyard_builder<double, Gudhi::vineyard::Default_vineyard_options, false>
using Vineyard_builder = Gudhi::vineyard::Vineyard_builder<Filtration_value>;
using Vine = Vineyard_builder::Vine_t;
using Cycle = Vineyard_builder::Cycle;
using Dimension = Vineyard_builder::Dimension;
using Index = Vineyard_builder::Index;
using Simplex_tree = Gudhi::Simplex_tree<>;  // uses double as Filtration_value by default
using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;
using Point = std::vector<Filtration_value>;

bool is_trivial(const Vine& vine, Filtration_value threshold)
{
  for (const auto& pair : vine.get_pairs()) {
    if (pair[1] - pair[0] > threshold) return false;
  }
  return true;
}

void print_vineyard(const std::vector<Vine>& vineyard, Filtration_value threshold, Dimension maxDim)
{
  std::cout << "Vines above threshold " << threshold << ":\n";
  unsigned int i = 0;
  for (const Vine& vine : vineyard) {
    auto dim = vine.get_dimension();
    if (dim <= maxDim && !is_trivial(vine, threshold)) {
      std::cout << "vine n° " << i << ": ";
      std::cout << "[" << dim << "] ";
      for (const auto& pair : vine.get_pairs()) {
        std::cout << "(" << pair[0] << ", " << pair[1] << ") ";
      }
      std::cout << "\n";
    }
    ++i;
  }
  std::cout << "\n";
}

void print_cycles(const std::vector<std::tuple<Cycle, Dimension, Index>>& cycles, int step)
{
  std::cout << "Non-trivial representative 1-cycles at step " << step << ":\n";
  for (const auto& cycle : cycles) {
    std::cout << "[" << get<1>(cycle) << "] ";
    for (const auto& idx : get<0>(cycle)) {
      std::cout << idx << " ";
    }
    std::cout << "(rep. vine n° " << get<2>(cycle) << " at step n° " << step << ")";
    std::cout << "\n";
  }
  std::cout << "\n";
}

void build_rips_complex(Simplex_tree& st, const std::vector<Point>& points)
{
  Rips_complex rips_complex_from_points(
      points, std::numeric_limits<Filtration_value>::infinity(), Gudhi::Euclidean_distance());
  rips_complex_from_points.create_complex(st, 2);
}

std::vector<Point> point_cloud_1()
{
  return {{-5.923219, 3.448312, 4.423653},  {-4.537442, 3.979718, 4.258528},   {-2.184740, 3.722509, 2.830755},
          {-1.536782, 3.765269, 1.400056},  {-1.468656, 3.360516, -0.794839},  {-2.160547, 3.318545, -2.527798},
          {-3.541939, 3.575503, -4.488310}, {-4.929805, 3.444428, -4.907573},  {-7.000496, 3.291212, -4.986328},
          {-9.106589, 3.596875, -4.436673}, {-10.469606, 3.413014, -2.375957}, {-10.736621, 3.255241, -0.744842},
          {-10.792444, 3.580073, 1.362501}, {-9.738022, 3.651807, 2.965904},   {-8.152821, 3.352247, 3.992848},
          {6.434311, 3.709249, -0.243290},  {5.702681, 5.899886, -0.546323},   {4.228264, 7.930254, -0.338698},
          {2.138856, 9.396997, -0.156139},  {-0.415187, 9.585111, -0.205344},  {-2.470828, 8.977603, -0.447381},
          {-4.477954, 6.951411, -0.307180}, {-5.555579, 4.885035, -0.152296},  {-5.690052, 2.472443, -0.491344},
          {-4.755301, 0.083429, -0.553459}, {-2.756040, -1.558657, -0.254191}, {-0.042665, -2.612909, -0.292575},
          {2.029499, -2.032840, -0.163618}, {4.460381, -0.778737, -0.461743},  {6.180998, 1.164147, -0.155191}};
}

std::vector<Point> point_cloud_2()
{
  return {{-8.584966, 3.473100, 4.471038},   {-7.199189, 4.004506, 4.305912},   {-4.846487, 3.747298, 2.878140},
          {-4.198530, 3.790057, 1.447441},   {-4.130403, 3.385305, -0.747455},  {-4.822294, 3.343333, -2.480413},
          {-6.203686, 3.600291, -4.440925},  {-7.591552, 3.469216, -4.860188},  {-9.662243, 3.316001, -4.938943},
          {-11.768336, 3.621663, -4.389288}, {-13.131353, 3.437802, -2.328572}, {-13.398368, 3.280029, -0.697457},
          {-13.454191, 3.604861, 1.409886},  {-12.399769, 3.676595, 3.013289},  {-10.814568, 3.377035, 4.040232},
          {9.053499, 3.623095, -0.168799},   {8.433733, 6.054553, -0.536405},   {7.290290, 8.035612, -0.458752},
          {5.062159, 9.549855, -0.184824},   {2.370435, 9.638249, -0.198651},   {0.124010, 9.032944, -0.404372},
          {-2.062278, 6.897104, -0.244886},  {-2.993217, 4.935229, -0.494346},  {-2.735058, 2.510149, -0.356950},
          {-2.054853, -0.015501, -0.282280}, {-0.093336, -1.671380, -0.435406}, {2.234651, -2.616742, -0.238271},
          {4.735187, -2.336024, -0.521192},  {6.991530, -1.088274, -0.489849},  {8.580836, 1.151733, -0.464560}};
}

std::vector<Point> point_cloud_3()
{
  return {{-11.059176, 3.512883, 4.418220},  {-9.673399, 4.044289, 4.253094},   {-7.320697, 3.787080, 2.825321},
          {-6.672740, 3.829840, 1.394622},   {-6.604613, 3.425087, -0.800273},  {-7.296504, 3.383116, -2.533231},
          {-8.677896, 3.640074, -4.493743},  {-10.065762, 3.508999, -4.913006}, {-12.136453, 3.355783, -4.991761},
          {-14.242546, 3.661445, -4.442107}, {-15.605563, 3.477585, -2.381390}, {-15.872579, 3.319812, -0.750275},
          {-15.928401, 3.644644, 1.357068},  {-14.873979, 3.716378, 2.960471},  {-13.288778, 3.416818, 3.987414},
          {11.778702, 3.401393, -0.451256},  {11.056639, 6.109902, -0.592193},  {9.767080, 8.204863, -0.473468},
          {7.140577, 9.114285, -0.233917},   {4.912468, 9.755665, -0.575986},   {2.341775, 8.735551, -0.301358},
          {0.328925, 7.253019, -0.149653},   {-0.701924, 4.727598, -0.165476},  {-0.416099, 2.290828, -0.121612},
          {0.539196, -0.068614, -0.333255},  {2.353096, -1.808696, -0.276212},  {4.591209, -2.286783, -0.269431},
          {7.527898, -2.437508, -0.344064},  {9.643542, -1.089942, -0.250708},  {11.187654, 1.202300, -0.148680}};
}

int main()
{
  // `true` to save non-trivial representative cycles and `1` to only store 1-cycles
  // if all non-trivial cycles are desired, just don't specify the last argument
  Vineyard_builder vyb(true, 1);

  // structures to fill
  Simplex_tree st;
  std::vector<std::vector<int>> boundaries;
  std::vector<int> dimensions;
  std::vector<Filtration_value> filtrationValues;

  // what is important here, is that the underlying complex remains the same at each step (if it is necessary to make
  // cells appear and/or disappear from one step to another, set their filtration value to infinity) and that
  // the boundaries remain in the same order!

  // builds a full rips complex from the first point cloud
  build_rips_complex(st, point_cloud_1());
  // retrieves the boundaries from the rips complex and their corresponding dimensions and filtration values
  Gudhi::vineyard::build_boundary_matrix_from_complex(st, boundaries, dimensions, filtrationValues);
  // initializes the vineyard with a first barcode
  vyb.initialize(boundaries, dimensions, filtrationValues);

  // prints non-trivial representative 1-cycles of last computed barcode
  print_cycles(vyb.get_latest_representative_cycles(), 0);

  st.clear();
  // builds a new rips complex from the second point cloud. As the point cloud represents exactly the same points
  // than before just with new filtration values and in the exact same order, with no threshold, the new complex will
  // be exactly the same than before just with new filtration values.
  build_rips_complex(st, point_cloud_2());
  filtrationValues.clear();
  // retrieves the new filtration values such that the order corresponds to the order of `boundaries`.
  // (possible because the rips complex method preserves the order of vertices and so the order is also preserved
  // in the final simplex tree)
  Gudhi::vineyard::build_boundary_matrix_from_complex(st, filtrationValues);
  // adds the second layer to the current vineyard
  vyb.update(filtrationValues);

  // prints representative 1-cycles of last computed barcode
  print_cycles(vyb.get_latest_representative_cycles(), 1);

  // same again with 3rd point cloud
  st.clear();
  build_rips_complex(st, point_cloud_3());
  filtrationValues.clear();
  Gudhi::vineyard::build_boundary_matrix_from_complex(st, filtrationValues);
  // adds the third layer to the current vineyard
  vyb.update(filtrationValues);

  // prints representative 1-cycles of last computed barcode
  print_cycles(vyb.get_latest_representative_cycles(), 2);
  // prints vines in vineyard with pairs of length > 3.5 and of maximal dimension 1
  print_vineyard(vyb.get_current_vineyard(), 3.5, 1);
}
