/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <tuple>
#include <vector>
#include <limits>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "zigzag_persistence"
#include <boost/test/unit_test.hpp>

#include <gudhi/filtered_zigzag_persistence.h>

struct Interval_comparator {
  Interval_comparator() {}
  template<class Interval_filtration>
  bool operator()(const Interval_filtration& p, const Interval_filtration& q) {
    if (p.dim != q.dim) {
      return p.dim < q.dim;
    }
    if (p.birth != q.birth) {
      return p.birth < q.birth;
    }
    return p.death < q.death;
  }
};

BOOST_AUTO_TEST_CASE(constructor) {
  using ZP1 = Gudhi::zigzag_persistence::Filtered_zigzag_persistence_with_storage<>;
  BOOST_CHECK_NO_THROW(ZP1 zp);
  BOOST_CHECK_NO_THROW(ZP1 zp(28));
  BOOST_CHECK_NO_THROW(ZP1 zp(28, 2));
  ZP1 zp1;
  BOOST_CHECK(zp1.get_persistence_diagram(0).empty());

  using ZP2 = Gudhi::zigzag_persistence::Filtered_zigzag_persistence<>;
  std::vector<std::tuple<int,double,double> > pairs;
  auto stream = [&](int dim, double birth, double death){ pairs.emplace_back(dim, birth, death); };
  BOOST_CHECK_NO_THROW(ZP2 zp(stream));
  BOOST_CHECK_NO_THROW(ZP2 zp(stream, 28));

  ZP2 zp2(stream);
  BOOST_CHECK(pairs.empty());
}

template<class ZP>
void test_barcode(ZP& zp, std::vector<typename ZP::Filtration_value_interval>& barcode) {
  auto bars = zp.get_persistence_diagram(0, true);
  std::stable_sort(bars.begin(), bars.end(), Interval_comparator());
  std::stable_sort(barcode.begin(), barcode.end(), Interval_comparator());
  auto it = barcode.begin();
  for (const auto& interval : bars) {
    BOOST_CHECK_EQUAL(interval.dim, it->dim);
    BOOST_CHECK_EQUAL(interval.birth, it->birth);
    BOOST_CHECK_EQUAL(interval.death, it->death);
    ++it;
  }
  BOOST_CHECK(it == barcode.end());
}

template <class ZP>
void test_indices(ZP& zp, std::vector<typename ZP::Index_interval>& indices,
                  std::vector<typename ZP::filtration_value>& indexToFil) {
  auto it = indices.begin();
  for (const auto& interval : zp.get_index_persistence_diagram()) {
    BOOST_CHECK_EQUAL(interval.dim, it->dim);
    BOOST_CHECK_EQUAL(interval.birth, it->birth);
    BOOST_CHECK_EQUAL(interval.death, it->death);
    auto p = zp.map_index_to_filtration_value(interval.birth, interval.death);
    BOOST_CHECK_EQUAL(p.first, indexToFil[interval.birth]);
    BOOST_CHECK_EQUAL(p.second, indexToFil[interval.death]);
    ++it;
  }
  BOOST_CHECK(it == indices.end());
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

std::vector<double> get_filtration_values() {
  return {0, 0, 0, 
          1, 1, 1, 
          2, 2, 2, 
          3, 3, 3, 3, 
          4, 
          5, 
          6, 6, 6, 
          7, 7, 7, 7, 7, 7, 
          8, 
          9, 9, 9, 
          10};
}

template<class ZP>
void test_filtered_zigzag_with_storage() {
  using face_handle = typename ZP::face_key;
  using filtration_value = typename ZP::filtration_value;
  using Interval_index = typename ZP::Index_interval;
  using Interval_filtration = typename ZP::Filtration_value_interval;

  ZP zp(28);
  std::vector<Interval_index> realIndices;
  std::vector<Interval_filtration> realBarcode;
  realIndices.reserve(13);
  realBarcode.reserve(9);

  std::vector<std::vector<face_handle> > simplices = get_boundaries();
  std::vector<filtration_value> filValues = get_filtration_values();

  for (unsigned int i = 0; i < 14; ++i) {
    zp.insert_face(i, simplices[i], simplices[i].size() == 0 ? 0 : simplices[i].size() - 1, filValues[i]);
  }

  realIndices.emplace_back(0, 1, 3);
  realIndices.emplace_back(0, 2, 4);
  realIndices.emplace_back(0, 7, 8);
  realIndices.emplace_back(1, 6, 10);
  realIndices.emplace_back(0, 9, 11);
  realIndices.emplace_back(1, 12, 13);

  realBarcode.emplace_back(0, 0, 1);
  realBarcode.emplace_back(0, 0, 1);
  realBarcode.emplace_back(1, 2, 3);
  realBarcode.emplace_back(1, 3, 4);

  for (unsigned int i = 14; i < 16; ++i) {
    auto id = simplices[i][0];
    zp.remove_face(id, simplices[id].size() == 0 ? 0 : simplices[id].size() - 1, filValues[i]);
  }

  for (unsigned int i = 16; i < 24; ++i) {
    zp.insert_face(i, simplices[i], simplices[i].size() == 0 ? 0 : simplices[i].size() - 1, filValues[i]);
  }

  realIndices.emplace_back(0, 5, 16);
  realIndices.emplace_back(1, 14, 17);
  realIndices.emplace_back(1, 15, 19);
  realIndices.emplace_back(1, 20, 21);
  realIndices.emplace_back(1, 18, 22);

  realBarcode.emplace_back(0, 1, 6);
  realBarcode.emplace_back(1, 5, 6);
  realBarcode.emplace_back(1, 6, 7);

  for (unsigned int i = 24; i < 27; ++i) {
    auto id = simplices[i][0];
    zp.remove_face(id, simplices[id].size() == 0 ? 0 : simplices[id].size() - 1, filValues[i]);
  }

  realIndices.emplace_back(1, 24, 25);
  realBarcode.emplace_back(1, 8, 9);

  zp.insert_face(27, simplices[27], simplices[27].size() == 0 ? 0 : simplices[27].size() - 1, filValues[27]);

  realIndices.emplace_back(2, 23, 27);
  realBarcode.emplace_back(2, 7, 9);

  auto id = simplices[28][0];
  zp.remove_face(id, simplices[id].size() == 0 ? 0 : simplices[id].size() - 1, filValues[28]);

  realBarcode.emplace_back(0, 0);
  realBarcode.emplace_back(0, 9);
  realBarcode.emplace_back(2, 10);

  test_indices(zp, realIndices, filValues);
  test_barcode(zp, realBarcode);
}

template<class ZP>
void test_filtered_zigzag_with_storage_max1() {
  using face_handle = typename ZP::face_key;
  using filtration_value = typename ZP::filtration_value;
  using Interval_index = typename ZP::Index_interval;
  using Interval_filtration = typename ZP::Filtration_value_interval;

  ZP zp(28, 1);
  std::vector<Interval_index> realIndices;
  std::vector<Interval_filtration> realBarcode;
  realIndices.reserve(5);
  realBarcode.reserve(3);

  std::vector<std::vector<face_handle> > simplices = get_boundaries();
  std::vector<filtration_value> filValues = get_filtration_values();

  for (unsigned int i = 0; i < 14; ++i) {
    zp.insert_face(i, simplices[i], simplices[i].size() == 0 ? 0 : simplices[i].size() - 1, filValues[i]);
  }

  realIndices.emplace_back(0, 1, 3);
  realIndices.emplace_back(0, 2, 4);
  realIndices.emplace_back(0, 7, 8);
  realIndices.emplace_back(0, 9, 11);

  realBarcode.emplace_back(0, 0, 1);
  realBarcode.emplace_back(0, 0, 1);

  for (unsigned int i = 14; i < 16; ++i) {
    auto id = simplices[i][0];
    zp.remove_face(id, simplices[id].size() == 0 ? 0 : simplices[id].size() - 1, filValues[i]);
  }

  for (unsigned int i = 16; i < 24; ++i) {
    zp.insert_face(i, simplices[i], simplices[i].size() == 0 ? 0 : simplices[i].size() - 1, filValues[i]);
  }

  realIndices.emplace_back(0, 5, 16);
  realBarcode.emplace_back(0, 1, 6);

  for (unsigned int i = 24; i < 27; ++i) {
    auto id = simplices[i][0];
    zp.remove_face(id, simplices[id].size() == 0 ? 0 : simplices[id].size() - 1, filValues[i]);
  }

  zp.insert_face(27, simplices[27], simplices[27].size() == 0 ? 0 : simplices[27].size() - 1, filValues[27]);
  auto id = simplices[28][0];
  zp.remove_face(id, simplices[id].size() == 0 ? 0 : simplices[id].size() - 1, filValues[28]);

  realBarcode.emplace_back(0, 0);
  realBarcode.emplace_back(0, 9);

  test_indices(zp, realIndices, filValues);
  test_barcode(zp, realBarcode);
}

BOOST_AUTO_TEST_CASE(filtered_zigzag_persistence_with_storage) {
  test_filtered_zigzag_with_storage<Gudhi::zigzag_persistence::Filtered_zigzag_persistence_with_storage<> >();
  test_filtered_zigzag_with_storage_max1<Gudhi::zigzag_persistence::Filtered_zigzag_persistence_with_storage<> >();
}

template<class ZP>
void test_filtered_zigzag() {
  using face_handle = typename ZP::face_key;
  using filtration_value = typename ZP::filtration_value;
  using dimension_type = typename ZP::dimension_type;
  using Interval = std::tuple<dimension_type, filtration_value, filtration_value>;

  Interval interval;
  ZP zp([&](dimension_type dim, filtration_value birth, filtration_value death){
    BOOST_CHECK_EQUAL(std::get<0>(interval), dim);
    BOOST_CHECK_EQUAL(std::get<1>(interval), birth);
    BOOST_CHECK_EQUAL(std::get<2>(interval), death);
  },28);

  std::vector<Interval> realBarcode;
  realBarcode.reserve(28);
  realBarcode.emplace_back(3, 0, 0);    //dummy
  realBarcode.emplace_back(3, 0, 1);    //dummy
  realBarcode.emplace_back(3, 0, 2);    //dummy
  realBarcode.emplace_back(0, 0, 1);    //1-3
  realBarcode.emplace_back(0, 0, 1);    //2-4
  realBarcode.emplace_back(3, 0, 5);    //dummy
  realBarcode.emplace_back(3, 0, 6);    //dummy
  realBarcode.emplace_back(3, 0, 7);    //dummy
  realBarcode.emplace_back(0, 7, 8);    //dummy
  realBarcode.emplace_back(3, 0, 9);    //dummy
  realBarcode.emplace_back(1, 2, 3);    //6-10
  realBarcode.emplace_back(0, 9, 11);   //dummy
  realBarcode.emplace_back(3, 0, 12);   //dummy
  realBarcode.emplace_back(1, 3, 4);    //12-13
  realBarcode.emplace_back(3, 0, 14);   //dummy
  realBarcode.emplace_back(3, 0, 15);   //dummy
  realBarcode.emplace_back(0, 1, 6);    //5-16
  realBarcode.emplace_back(1, 5, 6);    //14-17
  realBarcode.emplace_back(3, 0, 18);   //dummy
  realBarcode.emplace_back(1, 6, 7);    //15-19
  realBarcode.emplace_back(3, 0, 20);   //dummy
  realBarcode.emplace_back(1, 20, 21);  //dummy
  realBarcode.emplace_back(1, 18, 22);  //dummy
  realBarcode.emplace_back(3, 0, 23);   //dummy
  realBarcode.emplace_back(3, 0, 24);   //dummy
  realBarcode.emplace_back(1, 8, 9);    //24-25
  realBarcode.emplace_back(3, 0, 26);   //dummy
  realBarcode.emplace_back(2, 7, 9);    //23-27
  realBarcode.emplace_back(3, 0, 28);   //dummy

  std::vector<std::vector<face_handle> > simplices = get_boundaries();
  std::vector<filtration_value> filValues = get_filtration_values();

  for (unsigned int i = 0; i < 14; ++i) {
    interval = realBarcode[i];
    zp.insert_face(i, simplices[i], simplices[i].size() == 0 ? 0 : simplices[i].size() - 1, filValues[i]);
  }

  for (unsigned int i = 14; i < 16; ++i) {
    interval = realBarcode[i];
    auto id = simplices[i][0];
    zp.remove_face(id, simplices[id].size() == 0 ? 0 : simplices[id].size() - 1, filValues[i]);
  }

  for (unsigned int i = 16; i < 24; ++i) {
    interval = realBarcode[i];
    zp.insert_face(i, simplices[i], simplices[i].size() == 0 ? 0 : simplices[i].size() - 1, filValues[i]);
  }

  for (unsigned int i = 24; i < 27; ++i) {
    interval = realBarcode[i];
    auto id = simplices[i][0];
    zp.remove_face(id, simplices[id].size() == 0 ? 0 : simplices[id].size() - 1, filValues[i]);
  }

  interval = realBarcode[27];
  zp.insert_face(27, simplices[27], simplices[27].size() == 0 ? 0 : simplices[27].size() - 1, filValues[27]);

  interval = realBarcode[28];
  auto id = simplices[28][0];
  zp.remove_face(id, simplices[id].size() == 0 ? 0 : simplices[id].size() - 1, filValues[28]);

  //there is no real garantee on the order of the infinite bars
  std::vector<Interval> infiniteBars;
  zp.get_current_infinite_intervals([&](dimension_type dim, filtration_value birth) {
    infiniteBars.emplace_back(dim, birth, std::numeric_limits<filtration_value>::infinity());
  });

  realBarcode.clear();
  realBarcode.emplace_back(0, 0, std::numeric_limits<filtration_value>::infinity());
  realBarcode.emplace_back(0, 9, std::numeric_limits<filtration_value>::infinity());
  realBarcode.emplace_back(2, 10, std::numeric_limits<filtration_value>::infinity());

  std::sort(infiniteBars.begin(), infiniteBars.end());
  std::sort(realBarcode.begin(), realBarcode.end());
  auto it = realBarcode.begin();
  for (const auto& interval : infiniteBars) {
    BOOST_CHECK_EQUAL(std::get<0>(interval), std::get<0>(*it));
    BOOST_CHECK_EQUAL(std::get<1>(interval), std::get<1>(*it));
    ++it;
  }
  BOOST_CHECK(it == realBarcode.end());
}

template<class ZP>
void test_filtered_zigzag_max1() {
  using face_handle = typename ZP::face_key;
  using filtration_value = typename ZP::filtration_value;
  using dimension_type = typename ZP::dimension_type;
  using Interval = std::tuple<dimension_type, filtration_value, filtration_value>;

  Interval interval;
  ZP zp([&](dimension_type dim, filtration_value birth, filtration_value death){
    if (dim < 1){
      BOOST_CHECK_EQUAL(std::get<0>(interval), dim);
      BOOST_CHECK_EQUAL(std::get<1>(interval), birth);
      BOOST_CHECK_EQUAL(std::get<2>(interval), death);
    } else {
      BOOST_CHECK_NE(std::get<0>(interval), 0);
    }
  },28);

  std::vector<Interval> realBarcode;
  realBarcode.reserve(28);
  realBarcode.emplace_back(1, 0, 0);    //dummy
  realBarcode.emplace_back(1, 0, 1);    //dummy
  realBarcode.emplace_back(1, 0, 2);    //dummy
  realBarcode.emplace_back(0, 0, 1);    //1-3
  realBarcode.emplace_back(0, 0, 1);    //2-4
  realBarcode.emplace_back(1, 0, 5);    //dummy
  realBarcode.emplace_back(1, 0, 6);    //dummy
  realBarcode.emplace_back(1, 0, 7);    //dummy
  realBarcode.emplace_back(1, 7, 8);    //dummy
  realBarcode.emplace_back(1, 0, 9);    //dummy
  realBarcode.emplace_back(1, 2, 3);    //6-10
  realBarcode.emplace_back(1, 9, 11);   //dummy
  realBarcode.emplace_back(1, 0, 12);   //dummy
  realBarcode.emplace_back(1, 3, 4);    //12-13
  realBarcode.emplace_back(1, 0, 14);   //dummy
  realBarcode.emplace_back(1, 0, 15);   //dummy
  realBarcode.emplace_back(0, 1, 6);    //5-16
  realBarcode.emplace_back(1, 5, 6);    //14-17
  realBarcode.emplace_back(1, 0, 18);   //dummy
  realBarcode.emplace_back(1, 6, 7);    //15-19
  realBarcode.emplace_back(1, 0, 20);   //dummy
  realBarcode.emplace_back(1, 20, 21);  //dummy
  realBarcode.emplace_back(1, 18, 22);  //dummy
  realBarcode.emplace_back(1, 0, 23);   //dummy
  realBarcode.emplace_back(1, 0, 24);   //dummy
  realBarcode.emplace_back(1, 8, 9);    //24-25
  realBarcode.emplace_back(1, 0, 26);   //dummy
  realBarcode.emplace_back(2, 7, 9);    //23-27
  realBarcode.emplace_back(1, 0, 28);   //dummy

  std::vector<std::vector<face_handle> > simplices = get_boundaries();
  std::vector<filtration_value> filValues = get_filtration_values();

  for (unsigned int i = 0; i < 14; ++i) {
    interval = realBarcode[i];
    zp.insert_face(i, simplices[i], simplices[i].size() == 0 ? 0 : simplices[i].size() - 1, filValues[i]);
  }

  for (unsigned int i = 14; i < 16; ++i) {
    interval = realBarcode[i];
    auto id = simplices[i][0];
    zp.remove_face(id, simplices[id].size() == 0 ? 0 : simplices[id].size() - 1, filValues[i]);
  }

  for (unsigned int i = 16; i < 24; ++i) {
    interval = realBarcode[i];
    zp.insert_face(i, simplices[i], simplices[i].size() == 0 ? 0 : simplices[i].size() - 1, filValues[i]);
  }

  for (unsigned int i = 24; i < 27; ++i) {
    interval = realBarcode[i];
    auto id = simplices[i][0];
    zp.remove_face(id, simplices[id].size() == 0 ? 0 : simplices[id].size() - 1, filValues[i]);
  }

  interval = realBarcode[27];
  zp.insert_face(27, simplices[27], simplices[27].size() == 0 ? 0 : simplices[27].size() - 1, filValues[27]);

  interval = realBarcode[28];
  auto id = simplices[28][0];
  zp.remove_face(id, simplices[id].size() == 0 ? 0 : simplices[id].size() - 1, filValues[28]);

  //there is no real garantee on the order of the infinite bars
  std::vector<Interval> infiniteBars;
  zp.get_current_infinite_intervals([&](dimension_type dim, filtration_value birth) {
    if (dim < 1){
      infiniteBars.emplace_back(dim, birth, std::numeric_limits<filtration_value>::infinity());
    }
  });

  realBarcode.clear();
  realBarcode.emplace_back(0, 0, std::numeric_limits<filtration_value>::infinity());
  realBarcode.emplace_back(0, 9, std::numeric_limits<filtration_value>::infinity());

  std::sort(infiniteBars.begin(), infiniteBars.end());
  std::sort(realBarcode.begin(), realBarcode.end());
  auto it = realBarcode.begin();
  for (const auto& interval : infiniteBars) {
    BOOST_CHECK_EQUAL(std::get<0>(interval), std::get<0>(*it));
    BOOST_CHECK_EQUAL(std::get<1>(interval), std::get<1>(*it));
    ++it;
  }
  BOOST_CHECK(it == realBarcode.end());
}

BOOST_AUTO_TEST_CASE(filtered_zigzag_persistence) {
  test_filtered_zigzag<Gudhi::zigzag_persistence::Filtered_zigzag_persistence<> >();
  test_filtered_zigzag_max1<Gudhi::zigzag_persistence::Filtered_zigzag_persistence<> >();
}
