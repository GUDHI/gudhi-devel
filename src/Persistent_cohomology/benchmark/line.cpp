/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s): Marc Glisse
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Persistence_on_a_line.h>
#include <gudhi/Clock.h>

#include <vector>
#include <cstdlib>
#include <random>
#include <algorithm>

std::random_device rd;
std::mt19937 gen(rd());

float get_random_f()
{
  // values essentially unique
  std::uniform_real_distribution<float> dist(0., 1.);
  return dist(gen);
}
float get_random_i()
{
  // lots of repetitions in the values
  std::uniform_int_distribution<int> dist(0, 7);
  return dist(gen);
}

// $ prog [n_points [n_repeat]]
int main(int argc, char* argv[]) {

  using Cubical = Gudhi::cubical_complex::Bitmap_cubical_complex<Gudhi::cubical_complex::Bitmap_cubical_complex_base<float>>;
  using ST = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_fast_persistence>;
  using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;

  int N = 1000000;
  int repeat = 1;
  if (argc >= 2) N = atoi(argv[1]);
  if (argc >= 3) repeat = atoi(argv[2]);

  for (int rep = 0; rep < repeat; ++rep) {
    for (bool b : { false, true }) {
      std::vector<float> data(N);
      if (b) {
        std::clog << "Many repeated values\n";
        std::generate(data.begin(), data.end(), get_random_i);
      } else {
        std::clog << "Mostly unique values\n";
        std::generate(data.begin(), data.end(), get_random_f);
      }
      auto convert_result = [](auto& cplx, auto& pers){
        std::vector<std::pair<float,float>> res;
        for (auto p: pers.get_persistent_pairs()){
          res.emplace_back(cplx.filtration(std::get<0>(p)), cplx.filtration(std::get<1>(p)));
        }
        return res;
      };

      Gudhi::Clock clock_cub1("Cubical complex from top cells");
      Cubical cub1({{ (unsigned)data.size() }}, data, false);
      Gudhi::persistent_cohomology::Persistent_cohomology<Cubical, Field_Zp> pcub1(cub1);
      pcub1.init_coefficients(2);
      pcub1.compute_persistent_cohomology();
      clock_cub1.print();
      auto res_cub1 = convert_result(cub1, pcub1);
      std::sort(res_cub1.begin(), res_cub1.end());

      Gudhi::Clock clock_cub2("Cubical complex from vertices");
      Cubical cub2({{ (unsigned)data.size() }}, data, true);
      Gudhi::persistent_cohomology::Persistent_cohomology<Cubical, Field_Zp> pcub2(cub2);
      pcub2.init_coefficients(2);
      pcub2.compute_persistent_cohomology();
      clock_cub2.print();
      auto res_cub2 = convert_result(cub2, pcub2);
      std::sort(res_cub2.begin(), res_cub2.end());

      Gudhi::Clock clock_st("Simplex tree");
      ST st;
      for (int i=0; i < (int)data.size()-1; ++i) {
        st.insert_simplex({i}, data[i]);
        st.insert_simplex({i, i+1}, -std::numeric_limits<float>::infinity());
      }
      st.insert_simplex({(int)data.size()-1}, data.back());
      st.make_filtration_non_decreasing();
      Gudhi::persistent_cohomology::Persistent_cohomology<ST, Field_Zp> pst(st);
      pst.init_coefficients(2);
      pst.compute_persistent_cohomology();
      clock_st.print();
      auto res_st = convert_result(st, pst);
      std::sort(res_st.begin(), res_st.end());

      Gudhi::Clock clock_line("Special case for function on a line");
      std::vector<std::pair<float,float>> res_line;
      Gudhi::persistent_cohomology::compute_persistence_of_function_on_line(data,
          [&res_line](float a, float b){ res_line.emplace_back(a, b); });
      clock_line.print();
      std::sort(res_line.begin(), res_line.end());

      if (res_cub2 != res_cub1) std::exit(-1);
      if (res_st   != res_cub1) std::exit(-2);
      if (res_line != res_cub1) std::exit(-3);
    }
  }
}
