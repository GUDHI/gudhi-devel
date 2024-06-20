/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "zigzag_persistence"
#include <boost/test/unit_test.hpp>

#include <gudhi/Zigzag_persistence.h>

using ZP = Gudhi::zigzag_persistence::Zigzag_persistence<>;
// using ZP = Gudhi::zigzag_persistence::Zigzag_persistence<Gudhi::zigzag_persistence::Default_zigzag_options, false>;

struct Interval {
  Interval() {}
  Interval(int dim, ZP::index b, ZP::index d) : dim_(dim), b_(b), d_(d) {}

  int dim() const { return dim_; }
  int birth() const { return b_; }
  int death() const { return d_; }

private:
  int dim_;
  ZP::index b_;
  ZP::index d_;
};

BOOST_AUTO_TEST_CASE(constructor) {
  std::vector<Interval> pairs;
  auto stream = [&](int dim, ZP::index birth, ZP::index death){ pairs.emplace_back(dim, birth, death); };
  BOOST_CHECK_NO_THROW(ZP zp(stream));
  BOOST_CHECK_NO_THROW(ZP zp(stream, 28));

  ZP zp(stream);
  BOOST_CHECK(pairs.empty());
}

void test_indices(std::vector<Interval>& zp_indices, std::vector<Interval>& witness_indices) {
  auto it = witness_indices.begin();
  for (const Interval& interval : zp_indices) {
    BOOST_CHECK_EQUAL(interval.dim(), it->dim());
    BOOST_CHECK_EQUAL(interval.birth(), it->birth());
    BOOST_CHECK_EQUAL(interval.death(), it->death());
    ++it;
  }
  BOOST_CHECK(it == witness_indices.end());
}

std::vector<std::vector<int> > get_boundaries() {
  return {{},
          {},
          {},
          {0, 1},
          {0, 2},
          {},
          {1, 2},
          {},
          {5, 7},
          {},
          {3, 4, 6},
          {7, 9},
          {5, 9},
          {8, 11, 12},
          {10},                         // remove
          {13},                         // remove
          {1, 7},
          {3, 4, 6},
          {2, 7},
          {8, 11, 12},
          {0, 7},
          {4, 18, 20},
          {6, 16, 18},
          {3, 16, 20},
          {19},                         // remove
          {8},                          // remove
          {12},                         // remove
          {17, 21, 22, 23},
          {27}};                        // remove
}

BOOST_AUTO_TEST_CASE(zigzag_persistence_single) {
  std::vector<Interval> pairs;
  auto stream = [&](int dim, ZP::index birth, ZP::index death) { pairs.emplace_back(dim, birth, death); };
  auto stream_inf = [&](int dim, ZP::index birth) { pairs.emplace_back(dim, birth, -1); };
  ZP zp(stream, 28);
  std::vector<Interval> realIndices;
  realIndices.reserve(13);

  std::vector<std::vector<int> > simplices = get_boundaries();

  for (unsigned int i = 0; i < 14; ++i) {
    zp.insert_face(simplices[i], simplices[i].size() == 0 ? 0 : simplices[i].size() - 1);
  }

  realIndices.emplace_back(0, 1, 3);
  realIndices.emplace_back(0, 2, 4);
  realIndices.emplace_back(0, 7, 8);
  realIndices.emplace_back(1, 6, 10);
  realIndices.emplace_back(0, 9, 11);
  realIndices.emplace_back(1, 12, 13);

  for (unsigned int i = 14; i < 16; ++i) {
    auto id = simplices[i][0];
    zp.remove_face(id, simplices[id].size() == 0 ? 0 : simplices[id].size() - 1);
  }

  for (unsigned int i = 16; i < 24; ++i) {
    zp.insert_face(simplices[i], simplices[i].size() == 0 ? 0 : simplices[i].size() - 1);
  }

  realIndices.emplace_back(0, 5, 16);
  realIndices.emplace_back(1, 14, 17);
  realIndices.emplace_back(1, 15, 19);
  realIndices.emplace_back(1, 20, 21);
  realIndices.emplace_back(1, 18, 22);

  for (unsigned int i = 24; i < 27; ++i) {
    auto id = simplices[i][0];
    zp.remove_face(id, simplices[id].size() == 0 ? 0 : simplices[id].size() - 1);
  }

  realIndices.emplace_back(1, 24, 25);

  zp.insert_face(simplices[27], simplices[27].size() == 0 ? 0 : simplices[27].size() - 1);

  realIndices.emplace_back(2, 23, 27);

  auto id = simplices[28][0];
  zp.remove_face(id, simplices[id].size() == 0 ? 0 : simplices[id].size() - 1);

  realIndices.emplace_back(0, 0, -1);
  realIndices.emplace_back(0, 26, -1);
  realIndices.emplace_back(2, 28, -1);

  auto start = pairs.size();
  zp.get_current_infinit_intervals(stream_inf);
  std::sort(pairs.begin() + start, pairs.end(), [](const Interval& i1, const Interval& i2){
    if (i1.dim() != i2.dim()) return i1.dim() < i2.dim();
    return i1.birth() < i2.birth();
  });

  test_indices(pairs, realIndices);
}

BOOST_AUTO_TEST_CASE(zigzag_persistence_single_max1) {
  std::vector<Interval> pairs;
  auto stream = [&](int dim, ZP::index birth, ZP::index death) {
    if (dim < 1) pairs.emplace_back(dim, birth, death);
  };
  auto stream_inf = [&](int dim, ZP::index birth) {
    if (dim < 1) pairs.emplace_back(dim, birth, -1);
  };
  ZP zp(stream, 28);
  std::vector<Interval> realIndices;
  realIndices.reserve(5);

  std::vector<std::vector<int> > simplices = get_boundaries();

  for (unsigned int i = 0; i < 14; ++i) {
    zp.insert_face(simplices[i], simplices[i].size() == 0 ? 0 : simplices[i].size() - 1);
  }

  realIndices.emplace_back(0, 1, 3);
  realIndices.emplace_back(0, 2, 4);
  realIndices.emplace_back(0, 7, 8);
  realIndices.emplace_back(0, 9, 11);

  for (unsigned int i = 14; i < 16; ++i) {
    auto id = simplices[i][0];
    zp.remove_face(id, simplices[id].size() == 0 ? 0 : simplices[id].size() - 1);
  }

  for (unsigned int i = 16; i < 24; ++i) {
    zp.insert_face(simplices[i], simplices[i].size() == 0 ? 0 : simplices[i].size() - 1);
  }

  realIndices.emplace_back(0, 5, 16);

  for (unsigned int i = 24; i < 27; ++i) {
    auto id = simplices[i][0];
    zp.remove_face(id, simplices[id].size() == 0 ? 0 : simplices[id].size() - 1);
  }

  zp.insert_face(simplices[27], simplices[27].size() == 0 ? 0 : simplices[27].size() - 1);
  auto id = simplices[28][0];
  zp.remove_face(id, simplices[id].size() == 0 ? 0 : simplices[id].size() - 1);

  realIndices.emplace_back(0, 0, -1);
  realIndices.emplace_back(0, 26, -1);

  auto start = pairs.size();
  zp.get_current_infinit_intervals(stream_inf);
  std::sort(pairs.begin() + start, pairs.end(), [](const Interval& i1, const Interval& i2){
    if (i1.dim() != i2.dim()) return i1.dim() < i2.dim();
    return i1.birth() < i2.birth();
  });

  test_indices(pairs, realIndices);
}
