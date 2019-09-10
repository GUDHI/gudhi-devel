/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2015 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "cubical_complex"
#include <boost/test/unit_test.hpp>

#include <gudhi/reader_utils.h>
#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Persistent_cohomology.h>

// standard stuff
#include <iostream>
#include <sstream>
#include <vector>
#include <limits>

typedef Gudhi::cubical_complex::Bitmap_cubical_complex_base<double> Bitmap_cubical_complex_base;
typedef Gudhi::cubical_complex::Bitmap_cubical_complex<Bitmap_cubical_complex_base> Bitmap_cubical_complex;

typedef Gudhi::cubical_complex::Bitmap_cubical_complex_periodic_boundary_conditions_base<double>
    Bitmap_cubical_complex_periodic_boundary_conditions_base;
typedef Gudhi::cubical_complex::Bitmap_cubical_complex<Bitmap_cubical_complex_periodic_boundary_conditions_base>
    Bitmap_cubical_complex_periodic_boundary_conditions;

BOOST_AUTO_TEST_CASE(check_dimension) {
  std::vector<double> increasingFiltrationOfTopDimensionalCells({1, 2, 3, 4, 5, 6, 7, 8, 9});

  std::vector<unsigned> dimensions({3, 3});

  Bitmap_cubical_complex increasing(dimensions, increasingFiltrationOfTopDimensionalCells);
  BOOST_CHECK(increasing.dimension() == 2);
}

BOOST_AUTO_TEST_CASE(topDimensionalCellsIterator_test) {
  std::vector<double> expectedFiltrationValues1({0, 0, 0, 0, 100, 0, 0, 0, 0});

  std::vector<double> expectedFiltrationValues2({1, 2, 3, 4, 5, 6, 7, 8, 9});

  std::vector<double> increasingFiltrationOfTopDimensionalCells({1, 2, 3, 4, 5, 6, 7, 8, 9});

  std::vector<double> oneDimensionalCycle({0, 0, 0, 0, 100, 0, 0, 0, 0});

  std::vector<unsigned> dimensions({3, 3});

  Bitmap_cubical_complex increasing(dimensions, increasingFiltrationOfTopDimensionalCells);
  Bitmap_cubical_complex hole(dimensions, oneDimensionalCycle);

  int i = 0;
  for (Bitmap_cubical_complex::Top_dimensional_cells_iterator it = increasing.top_dimensional_cells_iterator_begin();
       it != increasing.top_dimensional_cells_iterator_end(); ++it) {
    BOOST_CHECK(increasing.get_cell_data(*it) == expectedFiltrationValues2[i]);
    ++i;
  }
  i = 0;
  for (Bitmap_cubical_complex::Top_dimensional_cells_iterator it = hole.top_dimensional_cells_iterator_begin();
       it != hole.top_dimensional_cells_iterator_end(); ++it) {
    BOOST_CHECK(hole.get_cell_data(*it) == expectedFiltrationValues1[i]);
    ++i;
  }
}

BOOST_AUTO_TEST_CASE(compute_boundary_test_1) {
  std::vector<double> boundary0;
  std::vector<double> boundary1;
  boundary1.push_back(0);
  boundary1.push_back(2);
  std::vector<double> boundary2;
  std::vector<double> boundary3;
  boundary3.push_back(2);
  boundary3.push_back(4);
  std::vector<double> boundary4;
  std::vector<double> boundary5;
  boundary5.push_back(4);
  boundary5.push_back(6);
  std::vector<double> boundary6;
  std::vector<double> boundary7;
  boundary7.push_back(0);
  boundary7.push_back(14);
  std::vector<double> boundary8;
  boundary8.push_back(1);
  boundary8.push_back(15);
  boundary8.push_back(9);
  boundary8.push_back(7);
  std::vector<double> boundary9;
  boundary9.push_back(2);
  boundary9.push_back(16);
  std::vector<double> boundary10;
  boundary10.push_back(3);
  boundary10.push_back(17);
  boundary10.push_back(11);
  boundary10.push_back(9);
  std::vector<double> boundary11;
  boundary11.push_back(4);
  boundary11.push_back(18);
  std::vector<double> boundary12;
  boundary12.push_back(5);
  boundary12.push_back(19);
  boundary12.push_back(13);
  boundary12.push_back(11);
  std::vector<double> boundary13;
  boundary13.push_back(6);
  boundary13.push_back(20);
  std::vector<double> boundary14;
  std::vector<double> boundary15;
  boundary15.push_back(14);
  boundary15.push_back(16);
  std::vector<double> boundary16;
  std::vector<double> boundary17;
  boundary17.push_back(16);
  boundary17.push_back(18);
  std::vector<double> boundary18;
  std::vector<double> boundary19;
  boundary19.push_back(18);
  boundary19.push_back(20);
  std::vector<double> boundary20;
  std::vector<double> boundary21;
  boundary21.push_back(14);
  boundary21.push_back(28);
  std::vector<double> boundary22;
  boundary22.push_back(15);
  boundary22.push_back(29);
  boundary22.push_back(23);
  boundary22.push_back(21);
  std::vector<double> boundary23;
  boundary23.push_back(16);
  boundary23.push_back(30);
  std::vector<double> boundary24;
  boundary24.push_back(17);
  boundary24.push_back(31);
  boundary24.push_back(25);
  boundary24.push_back(23);
  std::vector<double> boundary25;
  boundary25.push_back(18);
  boundary25.push_back(32);
  std::vector<double> boundary26;
  boundary26.push_back(19);
  boundary26.push_back(33);
  boundary26.push_back(27);
  boundary26.push_back(25);
  std::vector<double> boundary27;
  boundary27.push_back(20);
  boundary27.push_back(34);
  std::vector<double> boundary28;
  std::vector<double> boundary29;
  boundary29.push_back(28);
  boundary29.push_back(30);
  std::vector<double> boundary30;
  std::vector<double> boundary31;
  boundary31.push_back(30);
  boundary31.push_back(32);
  std::vector<double> boundary32;
  std::vector<double> boundary33;
  boundary33.push_back(32);
  boundary33.push_back(34);
  std::vector<double> boundary34;
  std::vector<double> boundary35;
  boundary35.push_back(28);
  boundary35.push_back(42);
  std::vector<double> boundary36;
  boundary36.push_back(29);
  boundary36.push_back(43);
  boundary36.push_back(37);
  boundary36.push_back(35);
  std::vector<double> boundary37;
  boundary37.push_back(30);
  boundary37.push_back(44);
  std::vector<double> boundary38;
  boundary38.push_back(31);
  boundary38.push_back(45);
  boundary38.push_back(39);
  boundary38.push_back(37);
  std::vector<double> boundary39;
  boundary39.push_back(32);
  boundary39.push_back(46);
  std::vector<double> boundary40;
  boundary40.push_back(33);
  boundary40.push_back(47);
  boundary40.push_back(41);
  boundary40.push_back(39);
  std::vector<double> boundary41;
  boundary41.push_back(34);
  boundary41.push_back(48);
  std::vector<double> boundary42;
  std::vector<double> boundary43;
  boundary43.push_back(42);
  boundary43.push_back(44);
  std::vector<double> boundary44;
  std::vector<double> boundary45;
  boundary45.push_back(44);
  boundary45.push_back(46);
  std::vector<double> boundary46;
  std::vector<double> boundary47;
  boundary47.push_back(46);
  boundary47.push_back(48);
  std::vector<double> boundary48;
  std::vector<std::vector<double> > boundaries;
  boundaries.push_back(boundary0);
  boundaries.push_back(boundary1);
  boundaries.push_back(boundary2);
  boundaries.push_back(boundary3);
  boundaries.push_back(boundary4);
  boundaries.push_back(boundary5);
  boundaries.push_back(boundary6);
  boundaries.push_back(boundary7);
  boundaries.push_back(boundary8);
  boundaries.push_back(boundary9);
  boundaries.push_back(boundary10);
  boundaries.push_back(boundary11);
  boundaries.push_back(boundary12);
  boundaries.push_back(boundary13);
  boundaries.push_back(boundary14);
  boundaries.push_back(boundary15);
  boundaries.push_back(boundary16);
  boundaries.push_back(boundary17);
  boundaries.push_back(boundary18);
  boundaries.push_back(boundary19);
  boundaries.push_back(boundary20);
  boundaries.push_back(boundary21);
  boundaries.push_back(boundary22);
  boundaries.push_back(boundary23);
  boundaries.push_back(boundary24);
  boundaries.push_back(boundary25);
  boundaries.push_back(boundary26);
  boundaries.push_back(boundary27);
  boundaries.push_back(boundary28);
  boundaries.push_back(boundary29);
  boundaries.push_back(boundary30);
  boundaries.push_back(boundary31);
  boundaries.push_back(boundary32);
  boundaries.push_back(boundary33);
  boundaries.push_back(boundary34);
  boundaries.push_back(boundary35);
  boundaries.push_back(boundary36);
  boundaries.push_back(boundary37);
  boundaries.push_back(boundary38);
  boundaries.push_back(boundary39);
  boundaries.push_back(boundary40);
  boundaries.push_back(boundary41);
  boundaries.push_back(boundary42);
  boundaries.push_back(boundary43);
  boundaries.push_back(boundary44);
  boundaries.push_back(boundary45);
  boundaries.push_back(boundary46);
  boundaries.push_back(boundary47);
  boundaries.push_back(boundary48);

  std::vector<double> increasingFiltrationOfTopDimensionalCells({1, 2, 3, 4, 5, 6, 7, 8, 9});

  std::vector<unsigned> dimensions({3, 3});

  Bitmap_cubical_complex increasing(dimensions, increasingFiltrationOfTopDimensionalCells);
  for (size_t i = 0; i != increasing.size(); ++i) {
    std::vector<size_t> bd = increasing.get_boundary_of_a_cell(i);
    for (size_t j = 0; j != bd.size(); ++j) {
      BOOST_CHECK(boundaries[i][j] == bd[j]);
    }
  }
}

BOOST_AUTO_TEST_CASE(compute_boundary_test_2) {
  std::vector<double> increasingFiltrationOfTopDimensionalCells({1, 2, 3, 4, 5, 6, 7, 8, 9});

  std::vector<unsigned> dimensions({3, 3});

  Bitmap_cubical_complex increasing(dimensions, increasingFiltrationOfTopDimensionalCells);

  std::vector<double> coboundaryElements;
  coboundaryElements.push_back(7);
  coboundaryElements.push_back(1);
  coboundaryElements.push_back(8);
  coboundaryElements.push_back(9);
  coboundaryElements.push_back(1);
  coboundaryElements.push_back(3);
  coboundaryElements.push_back(10);
  coboundaryElements.push_back(11);
  coboundaryElements.push_back(3);
  coboundaryElements.push_back(5);
  coboundaryElements.push_back(12);
  coboundaryElements.push_back(13);
  coboundaryElements.push_back(5);
  coboundaryElements.push_back(8);
  coboundaryElements.push_back(8);
  coboundaryElements.push_back(10);
  coboundaryElements.push_back(10);
  coboundaryElements.push_back(12);
  coboundaryElements.push_back(12);
  coboundaryElements.push_back(7);
  coboundaryElements.push_back(21);
  coboundaryElements.push_back(15);
  coboundaryElements.push_back(8);
  coboundaryElements.push_back(22);
  coboundaryElements.push_back(9);
  coboundaryElements.push_back(23);
  coboundaryElements.push_back(15);
  coboundaryElements.push_back(17);
  coboundaryElements.push_back(10);
  coboundaryElements.push_back(24);
  coboundaryElements.push_back(11);
  coboundaryElements.push_back(25);
  coboundaryElements.push_back(17);
  coboundaryElements.push_back(19);
  coboundaryElements.push_back(12);
  coboundaryElements.push_back(26);
  coboundaryElements.push_back(13);
  coboundaryElements.push_back(27);
  coboundaryElements.push_back(19);
  coboundaryElements.push_back(22);
  coboundaryElements.push_back(22);
  coboundaryElements.push_back(24);
  coboundaryElements.push_back(24);
  coboundaryElements.push_back(26);
  coboundaryElements.push_back(26);
  coboundaryElements.push_back(21);
  coboundaryElements.push_back(35);
  coboundaryElements.push_back(29);
  coboundaryElements.push_back(22);
  coboundaryElements.push_back(36);
  coboundaryElements.push_back(23);
  coboundaryElements.push_back(37);
  coboundaryElements.push_back(29);
  coboundaryElements.push_back(31);
  coboundaryElements.push_back(24);
  coboundaryElements.push_back(38);
  coboundaryElements.push_back(25);
  coboundaryElements.push_back(39);
  coboundaryElements.push_back(31);
  coboundaryElements.push_back(33);
  coboundaryElements.push_back(26);
  coboundaryElements.push_back(40);
  coboundaryElements.push_back(27);
  coboundaryElements.push_back(41);
  coboundaryElements.push_back(33);
  coboundaryElements.push_back(36);
  coboundaryElements.push_back(36);
  coboundaryElements.push_back(38);
  coboundaryElements.push_back(38);
  coboundaryElements.push_back(40);
  coboundaryElements.push_back(40);
  coboundaryElements.push_back(35);
  coboundaryElements.push_back(43);
  coboundaryElements.push_back(36);
  coboundaryElements.push_back(37);
  coboundaryElements.push_back(43);
  coboundaryElements.push_back(45);
  coboundaryElements.push_back(38);
  coboundaryElements.push_back(39);
  coboundaryElements.push_back(45);
  coboundaryElements.push_back(47);
  coboundaryElements.push_back(40);
  coboundaryElements.push_back(41);
  coboundaryElements.push_back(47);

  size_t number = 0;
  for (size_t i = 0; i != increasing.size(); ++i) {
    std::vector<size_t> bd = increasing.get_coboundary_of_a_cell(i);
    for (size_t j = 0; j != bd.size(); ++j) {
      BOOST_CHECK(coboundaryElements[number] == bd[j]);
      ++number;
    }
  }
}

BOOST_AUTO_TEST_CASE(compute_boundary_test_3) {
  std::vector<double> increasingFiltrationOfTopDimensionalCells({1, 2, 3, 4, 5, 6, 7, 8, 9});

  std::vector<unsigned> dimensions({3, 3});

  Bitmap_cubical_complex increasing(dimensions, increasingFiltrationOfTopDimensionalCells);

  std::vector<unsigned> dim;
  dim.push_back(0);
  dim.push_back(1);
  dim.push_back(0);
  dim.push_back(1);
  dim.push_back(0);
  dim.push_back(1);
  dim.push_back(0);
  dim.push_back(1);
  dim.push_back(2);
  dim.push_back(1);
  dim.push_back(2);
  dim.push_back(1);
  dim.push_back(2);
  dim.push_back(1);
  dim.push_back(0);
  dim.push_back(1);
  dim.push_back(0);
  dim.push_back(1);
  dim.push_back(0);
  dim.push_back(1);
  dim.push_back(0);
  dim.push_back(1);
  dim.push_back(2);
  dim.push_back(1);
  dim.push_back(2);
  dim.push_back(1);
  dim.push_back(2);
  dim.push_back(1);
  dim.push_back(0);
  dim.push_back(1);
  dim.push_back(0);
  dim.push_back(1);
  dim.push_back(0);
  dim.push_back(1);
  dim.push_back(0);
  dim.push_back(1);
  dim.push_back(2);
  dim.push_back(1);
  dim.push_back(2);
  dim.push_back(1);
  dim.push_back(2);
  dim.push_back(1);
  dim.push_back(0);
  dim.push_back(1);
  dim.push_back(0);
  dim.push_back(1);
  dim.push_back(0);
  dim.push_back(1);
  dim.push_back(0);

  for (size_t i = 0; i != increasing.size(); ++i) {
    BOOST_CHECK(increasing.get_dimension_of_a_cell(i) == dim[i]);
  }
}

BOOST_AUTO_TEST_CASE(Filtration_simplex_iterator_test) {
  std::vector<double> increasingFiltrationOfTopDimensionalCells({1, 2, 3, 4, 5, 6, 7, 8, 9});

  std::vector<unsigned> dimensions({3, 3});

  Bitmap_cubical_complex increasing(dimensions, increasingFiltrationOfTopDimensionalCells);

  std::vector<unsigned> dim;
  dim.push_back(0);
  dim.push_back(0);
  dim.push_back(0);
  dim.push_back(0);
  dim.push_back(1);
  dim.push_back(1);
  dim.push_back(1);
  dim.push_back(1);
  dim.push_back(2);
  dim.push_back(0);
  dim.push_back(0);
  dim.push_back(1);
  dim.push_back(1);
  dim.push_back(1);
  dim.push_back(2);
  dim.push_back(0);
  dim.push_back(0);
  dim.push_back(1);
  dim.push_back(1);
  dim.push_back(1);
  dim.push_back(2);
  dim.push_back(0);
  dim.push_back(0);
  dim.push_back(1);
  dim.push_back(1);
  dim.push_back(1);
  dim.push_back(2);
  dim.push_back(0);
  dim.push_back(1);
  dim.push_back(1);
  dim.push_back(2);
  dim.push_back(0);
  dim.push_back(1);
  dim.push_back(1);
  dim.push_back(2);
  dim.push_back(0);
  dim.push_back(0);
  dim.push_back(1);
  dim.push_back(1);
  dim.push_back(1);
  dim.push_back(2);
  dim.push_back(0);
  dim.push_back(1);
  dim.push_back(1);
  dim.push_back(2);
  dim.push_back(0);
  dim.push_back(1);
  dim.push_back(1);
  dim.push_back(2);

  std::vector<double> fil;
  fil.push_back(1);
  fil.push_back(1);
  fil.push_back(1);
  fil.push_back(1);
  fil.push_back(1);
  fil.push_back(1);
  fil.push_back(1);
  fil.push_back(1);
  fil.push_back(1);
  fil.push_back(2);
  fil.push_back(2);
  fil.push_back(2);
  fil.push_back(2);
  fil.push_back(2);
  fil.push_back(2);
  fil.push_back(3);
  fil.push_back(3);
  fil.push_back(3);
  fil.push_back(3);
  fil.push_back(3);
  fil.push_back(3);
  fil.push_back(4);
  fil.push_back(4);
  fil.push_back(4);
  fil.push_back(4);
  fil.push_back(4);
  fil.push_back(4);
  fil.push_back(5);
  fil.push_back(5);
  fil.push_back(5);
  fil.push_back(5);
  fil.push_back(6);
  fil.push_back(6);
  fil.push_back(6);
  fil.push_back(6);
  fil.push_back(7);
  fil.push_back(7);
  fil.push_back(7);
  fil.push_back(7);
  fil.push_back(7);
  fil.push_back(7);
  fil.push_back(8);
  fil.push_back(8);
  fil.push_back(8);
  fil.push_back(8);
  fil.push_back(9);
  fil.push_back(9);
  fil.push_back(9);
  fil.push_back(9);

  Bitmap_cubical_complex::Filtration_simplex_range range = increasing.filtration_simplex_range();
  size_t position = 0;
  for (Bitmap_cubical_complex::Filtration_simplex_iterator it = range.begin(); it != range.end(); ++it) {
    BOOST_CHECK(increasing.dimension(*it) == dim[position]);
    BOOST_CHECK(increasing.filtration(*it) == fil[position]);
    ++position;
  }
}

BOOST_AUTO_TEST_CASE(boudary_operator_2d_bitmap_with_periodic_bcond) {
  std::vector<double> filtration({0, 0, 0, 0});

  std::vector<unsigned> dimensions({2, 2});

  std::vector<bool> periodic_directions({true, true});

  Bitmap_cubical_complex_periodic_boundary_conditions cmplx(dimensions, filtration, periodic_directions);
  BOOST_CHECK(cmplx.dimension() == 2);

  std::vector<double> boundary0;
  std::vector<double> boundary1;
  boundary1.push_back(2);
  boundary1.push_back(0);
  std::vector<double> boundary2;
  std::vector<double> boundary3;
  boundary3.push_back(0);
  boundary3.push_back(2);
  std::vector<double> boundary4;
  boundary4.push_back(8);
  boundary4.push_back(0);
  std::vector<double> boundary5;
  boundary5.push_back(9);
  boundary5.push_back(1);
  boundary5.push_back(4);
  boundary5.push_back(6);
  std::vector<double> boundary6;
  boundary6.push_back(10);
  boundary6.push_back(2);
  std::vector<double> boundary7;
  boundary7.push_back(11);
  boundary7.push_back(3);
  boundary7.push_back(6);
  boundary7.push_back(4);
  std::vector<double> boundary8;
  std::vector<double> boundary9;
  boundary9.push_back(10);
  boundary9.push_back(8);
  std::vector<double> boundary10;
  std::vector<double> boundary11;
  boundary11.push_back(8);
  boundary11.push_back(10);
  std::vector<double> boundary12;
  boundary12.push_back(0);
  boundary12.push_back(8);
  std::vector<double> boundary13;
  boundary13.push_back(1);
  boundary13.push_back(9);
  boundary13.push_back(12);
  boundary13.push_back(14);
  std::vector<double> boundary14;
  boundary14.push_back(2);
  boundary14.push_back(10);
  std::vector<double> boundary15;
  boundary15.push_back(3);
  boundary15.push_back(11);
  boundary15.push_back(14);
  boundary15.push_back(12);

  std::vector<std::vector<double> > boundaries;
  boundaries.push_back(boundary0);
  boundaries.push_back(boundary1);
  boundaries.push_back(boundary2);
  boundaries.push_back(boundary3);
  boundaries.push_back(boundary4);
  boundaries.push_back(boundary5);
  boundaries.push_back(boundary6);
  boundaries.push_back(boundary7);
  boundaries.push_back(boundary8);
  boundaries.push_back(boundary9);
  boundaries.push_back(boundary10);
  boundaries.push_back(boundary11);
  boundaries.push_back(boundary12);
  boundaries.push_back(boundary13);
  boundaries.push_back(boundary14);
  boundaries.push_back(boundary15);

  for (size_t i = 0; i != cmplx.size(); ++i) {
    std::vector<size_t> bd = cmplx.get_boundary_of_a_cell(i);
    for (size_t j = 0; j != bd.size(); ++j) {
      BOOST_CHECK(boundaries[i][j] == bd[j]);
    }
  }
}

BOOST_AUTO_TEST_CASE(coboudary_operator_2d_bitmap_with_periodic_bcond) {
  std::vector<double> filtration({0, 0, 0, 0});

  std::vector<unsigned> dimensions({2, 2});

  std::vector<bool> periodic_directions({true, true});

  Bitmap_cubical_complex_periodic_boundary_conditions cmplx(dimensions, filtration, periodic_directions);
  BOOST_CHECK(cmplx.dimension() == 2);

  std::vector<double> coboundary0;
  coboundary0.push_back(4);
  coboundary0.push_back(12);
  coboundary0.push_back(1);
  coboundary0.push_back(3);
  std::vector<double> coboundary1;
  coboundary1.push_back(5);
  coboundary1.push_back(13);
  std::vector<double> coboundary2;
  coboundary2.push_back(6);
  coboundary2.push_back(14);
  coboundary2.push_back(1);
  coboundary2.push_back(3);
  std::vector<double> coboundary3;
  coboundary3.push_back(7);
  coboundary3.push_back(15);
  std::vector<double> coboundary4;
  coboundary4.push_back(5);
  coboundary4.push_back(7);
  std::vector<double> coboundary5;
  std::vector<double> coboundary6;
  coboundary6.push_back(5);
  coboundary6.push_back(7);
  std::vector<double> coboundary7;
  std::vector<double> coboundary8;
  coboundary8.push_back(4);
  coboundary8.push_back(12);
  coboundary8.push_back(9);
  coboundary8.push_back(11);
  std::vector<double> coboundary9;
  coboundary9.push_back(5);
  coboundary9.push_back(13);
  std::vector<double> coboundary10;
  coboundary10.push_back(6);
  coboundary10.push_back(14);
  coboundary10.push_back(9);
  coboundary10.push_back(11);
  std::vector<double> coboundary11;
  coboundary11.push_back(7);
  coboundary11.push_back(15);
  std::vector<double> coboundary12;
  coboundary12.push_back(13);
  coboundary12.push_back(15);
  std::vector<double> coboundary13;
  std::vector<double> coboundary14;
  coboundary14.push_back(13);
  coboundary14.push_back(15);
  std::vector<double> coboundary15;

  std::vector<std::vector<double> > coboundaries;
  coboundaries.push_back(coboundary0);
  coboundaries.push_back(coboundary1);
  coboundaries.push_back(coboundary2);
  coboundaries.push_back(coboundary3);
  coboundaries.push_back(coboundary4);
  coboundaries.push_back(coboundary5);
  coboundaries.push_back(coboundary6);
  coboundaries.push_back(coboundary7);
  coboundaries.push_back(coboundary8);
  coboundaries.push_back(coboundary9);
  coboundaries.push_back(coboundary10);
  coboundaries.push_back(coboundary11);
  coboundaries.push_back(coboundary12);
  coboundaries.push_back(coboundary13);
  coboundaries.push_back(coboundary14);
  coboundaries.push_back(coboundary15);

  for (size_t i = 0; i != cmplx.size(); ++i) {
    std::vector<size_t> cbd = cmplx.get_coboundary_of_a_cell(i);
    for (size_t j = 0; j != cbd.size(); ++j) {
      BOOST_CHECK(coboundaries[i][j] == cbd[j]);
    }
  }
}

BOOST_AUTO_TEST_CASE(bitmap_2d_with_periodic_bcond_filtration) {
  std::vector<double> filtrationOrg({0, 1, 2, 3});

  std::vector<unsigned> dimensions({2, 2});

  std::vector<bool> periodic_directions({true, true});

  Bitmap_cubical_complex_periodic_boundary_conditions cmplx(dimensions, filtrationOrg, periodic_directions);
  BOOST_CHECK(cmplx.dimension() == 2);

  std::vector<double> filtration;
  filtration.push_back(0);  // 0
  filtration.push_back(0);  // 1
  filtration.push_back(0);  // 2
  filtration.push_back(1);  // 3
  filtration.push_back(0);  // 4
  filtration.push_back(0);  // 5
  filtration.push_back(0);  // 6
  filtration.push_back(1);  // 7
  filtration.push_back(0);  // 8
  filtration.push_back(0);  // 9
  filtration.push_back(0);  // 10
  filtration.push_back(1);  // 11
  filtration.push_back(2);  // 12
  filtration.push_back(2);  // 13
  filtration.push_back(2);  // 14
  filtration.push_back(3);  // 15

  for (size_t i = 0; i != cmplx.size(); ++i) {
    BOOST_CHECK(filtration[i] == cmplx.get_cell_data(i));
  }
}

BOOST_AUTO_TEST_CASE(all_cells_iterator_and_boundary_iterators_in_Bitmap_cubical_complex_base_check) {
  std::vector<double> expected_filtration;
  expected_filtration.push_back(0);
  expected_filtration.push_back(0);
  expected_filtration.push_back(0);
  expected_filtration.push_back(1);
  expected_filtration.push_back(1);
  expected_filtration.push_back(0);
  expected_filtration.push_back(0);
  expected_filtration.push_back(0);
  expected_filtration.push_back(1);
  expected_filtration.push_back(1);
  expected_filtration.push_back(0);
  expected_filtration.push_back(0);
  expected_filtration.push_back(0);
  expected_filtration.push_back(1);
  expected_filtration.push_back(1);
  expected_filtration.push_back(2);
  expected_filtration.push_back(2);
  expected_filtration.push_back(2);
  expected_filtration.push_back(3);
  expected_filtration.push_back(3);
  expected_filtration.push_back(2);
  expected_filtration.push_back(2);
  expected_filtration.push_back(2);
  expected_filtration.push_back(3);
  expected_filtration.push_back(3);

  std::vector<unsigned> expected_dimension;
  expected_dimension.push_back(0);
  expected_dimension.push_back(1);
  expected_dimension.push_back(0);
  expected_dimension.push_back(1);
  expected_dimension.push_back(0);
  expected_dimension.push_back(1);
  expected_dimension.push_back(2);
  expected_dimension.push_back(1);
  expected_dimension.push_back(2);
  expected_dimension.push_back(1);
  expected_dimension.push_back(0);
  expected_dimension.push_back(1);
  expected_dimension.push_back(0);
  expected_dimension.push_back(1);
  expected_dimension.push_back(0);
  expected_dimension.push_back(1);
  expected_dimension.push_back(2);
  expected_dimension.push_back(1);
  expected_dimension.push_back(2);
  expected_dimension.push_back(1);
  expected_dimension.push_back(0);
  expected_dimension.push_back(1);
  expected_dimension.push_back(0);
  expected_dimension.push_back(1);
  expected_dimension.push_back(0);

  std::vector<size_t> expected_boundary;
  expected_boundary.push_back(0);
  expected_boundary.push_back(2);
  expected_boundary.push_back(2);
  expected_boundary.push_back(4);
  expected_boundary.push_back(0);
  expected_boundary.push_back(10);
  expected_boundary.push_back(1);
  expected_boundary.push_back(11);
  expected_boundary.push_back(7);
  expected_boundary.push_back(5);
  expected_boundary.push_back(2);
  expected_boundary.push_back(12);
  expected_boundary.push_back(3);
  expected_boundary.push_back(13);
  expected_boundary.push_back(9);
  expected_boundary.push_back(7);
  expected_boundary.push_back(4);
  expected_boundary.push_back(14);
  expected_boundary.push_back(10);
  expected_boundary.push_back(12);
  expected_boundary.push_back(12);
  expected_boundary.push_back(14);
  expected_boundary.push_back(10);
  expected_boundary.push_back(20);
  expected_boundary.push_back(11);
  expected_boundary.push_back(21);
  expected_boundary.push_back(17);
  expected_boundary.push_back(15);
  expected_boundary.push_back(12);
  expected_boundary.push_back(22);
  expected_boundary.push_back(13);
  expected_boundary.push_back(23);
  expected_boundary.push_back(19);
  expected_boundary.push_back(17);
  expected_boundary.push_back(14);
  expected_boundary.push_back(24);
  expected_boundary.push_back(20);
  expected_boundary.push_back(22);
  expected_boundary.push_back(22);
  expected_boundary.push_back(24);

  std::vector<size_t> expected_coboundary;
  expected_coboundary.push_back(5);
  expected_coboundary.push_back(1);
  expected_coboundary.push_back(6);
  expected_coboundary.push_back(7);
  expected_coboundary.push_back(1);
  expected_coboundary.push_back(3);
  expected_coboundary.push_back(8);
  expected_coboundary.push_back(9);
  expected_coboundary.push_back(3);
  expected_coboundary.push_back(6);
  expected_coboundary.push_back(6);
  expected_coboundary.push_back(8);
  expected_coboundary.push_back(8);
  expected_coboundary.push_back(5);
  expected_coboundary.push_back(15);
  expected_coboundary.push_back(11);
  expected_coboundary.push_back(6);
  expected_coboundary.push_back(16);
  expected_coboundary.push_back(7);
  expected_coboundary.push_back(17);
  expected_coboundary.push_back(11);
  expected_coboundary.push_back(13);
  expected_coboundary.push_back(8);
  expected_coboundary.push_back(18);
  expected_coboundary.push_back(9);
  expected_coboundary.push_back(19);
  expected_coboundary.push_back(13);
  expected_coboundary.push_back(16);
  expected_coboundary.push_back(16);
  expected_coboundary.push_back(18);
  expected_coboundary.push_back(18);
  expected_coboundary.push_back(15);
  expected_coboundary.push_back(21);
  expected_coboundary.push_back(16);
  expected_coboundary.push_back(17);
  expected_coboundary.push_back(21);
  expected_coboundary.push_back(23);
  expected_coboundary.push_back(18);
  expected_coboundary.push_back(19);
  expected_coboundary.push_back(23);

  std::vector<unsigned> sizes(2);
  sizes[0] = 2;
  sizes[1] = 2;

  std::vector<double> data(4);
  data[0] = 0;
  data[1] = 1;
  data[2] = 2;
  data[3] = 3;

  Bitmap_cubical_complex_base ba(sizes, data);
  int i = 0;
  int bd_it = 0;
  int cbd_it = 0;
  for (Bitmap_cubical_complex_base::All_cells_iterator it = ba.all_cells_iterator_begin();
       it != ba.all_cells_iterator_end(); ++it) {
    BOOST_CHECK(expected_filtration[i] == ba.get_cell_data(*it));
    BOOST_CHECK(expected_dimension[i] == ba.get_dimension_of_a_cell(*it));

    Bitmap_cubical_complex_base::Boundary_range bdrange = ba.boundary_range(*it);
    for (Bitmap_cubical_complex_base::Boundary_iterator bd = bdrange.begin(); bd != bdrange.end(); ++bd) {
      BOOST_CHECK(expected_boundary[bd_it] == *bd);
      ++bd_it;
    }

    Bitmap_cubical_complex_base::Coboundary_range cbdrange = ba.coboundary_range(*it);
    for (Bitmap_cubical_complex_base::Coboundary_iterator cbd = cbdrange.begin(); cbd != cbdrange.end(); ++cbd) {
      BOOST_CHECK(expected_coboundary[cbd_it] == *cbd);
      ++cbd_it;
    }
    ++i;
  }
}

BOOST_AUTO_TEST_CASE(all_cells_iterator_and_boundary_iterators_in_Bitmap_cubical_complex_base_check_range_check_2) {
  std::vector<double> expected_filtration;
  expected_filtration.push_back(0);
  expected_filtration.push_back(0);
  expected_filtration.push_back(0);
  expected_filtration.push_back(1);
  expected_filtration.push_back(1);
  expected_filtration.push_back(0);
  expected_filtration.push_back(0);
  expected_filtration.push_back(0);
  expected_filtration.push_back(1);
  expected_filtration.push_back(1);
  expected_filtration.push_back(0);
  expected_filtration.push_back(0);
  expected_filtration.push_back(0);
  expected_filtration.push_back(1);
  expected_filtration.push_back(1);
  expected_filtration.push_back(2);
  expected_filtration.push_back(2);
  expected_filtration.push_back(2);
  expected_filtration.push_back(3);
  expected_filtration.push_back(3);
  expected_filtration.push_back(2);
  expected_filtration.push_back(2);
  expected_filtration.push_back(2);
  expected_filtration.push_back(3);
  expected_filtration.push_back(3);

  std::vector<unsigned> expected_dimension;
  expected_dimension.push_back(0);
  expected_dimension.push_back(1);
  expected_dimension.push_back(0);
  expected_dimension.push_back(1);
  expected_dimension.push_back(0);
  expected_dimension.push_back(1);
  expected_dimension.push_back(2);
  expected_dimension.push_back(1);
  expected_dimension.push_back(2);
  expected_dimension.push_back(1);
  expected_dimension.push_back(0);
  expected_dimension.push_back(1);
  expected_dimension.push_back(0);
  expected_dimension.push_back(1);
  expected_dimension.push_back(0);
  expected_dimension.push_back(1);
  expected_dimension.push_back(2);
  expected_dimension.push_back(1);
  expected_dimension.push_back(2);
  expected_dimension.push_back(1);
  expected_dimension.push_back(0);
  expected_dimension.push_back(1);
  expected_dimension.push_back(0);
  expected_dimension.push_back(1);
  expected_dimension.push_back(0);

  std::vector<size_t> expected_boundary;
  expected_boundary.push_back(0);
  expected_boundary.push_back(2);
  expected_boundary.push_back(2);
  expected_boundary.push_back(4);
  expected_boundary.push_back(0);
  expected_boundary.push_back(10);
  expected_boundary.push_back(1);
  expected_boundary.push_back(11);
  expected_boundary.push_back(7);
  expected_boundary.push_back(5);
  expected_boundary.push_back(2);
  expected_boundary.push_back(12);
  expected_boundary.push_back(3);
  expected_boundary.push_back(13);
  expected_boundary.push_back(9);
  expected_boundary.push_back(7);
  expected_boundary.push_back(4);
  expected_boundary.push_back(14);
  expected_boundary.push_back(10);
  expected_boundary.push_back(12);
  expected_boundary.push_back(12);
  expected_boundary.push_back(14);
  expected_boundary.push_back(10);
  expected_boundary.push_back(20);
  expected_boundary.push_back(11);
  expected_boundary.push_back(21);
  expected_boundary.push_back(17);
  expected_boundary.push_back(15);
  expected_boundary.push_back(12);
  expected_boundary.push_back(22);
  expected_boundary.push_back(13);
  expected_boundary.push_back(23);
  expected_boundary.push_back(19);
  expected_boundary.push_back(17);
  expected_boundary.push_back(14);
  expected_boundary.push_back(24);
  expected_boundary.push_back(20);
  expected_boundary.push_back(22);
  expected_boundary.push_back(22);
  expected_boundary.push_back(24);

  std::vector<size_t> expected_coboundary;
  expected_coboundary.push_back(5);
  expected_coboundary.push_back(1);
  expected_coboundary.push_back(6);
  expected_coboundary.push_back(7);
  expected_coboundary.push_back(1);
  expected_coboundary.push_back(3);
  expected_coboundary.push_back(8);
  expected_coboundary.push_back(9);
  expected_coboundary.push_back(3);
  expected_coboundary.push_back(6);
  expected_coboundary.push_back(6);
  expected_coboundary.push_back(8);
  expected_coboundary.push_back(8);
  expected_coboundary.push_back(5);
  expected_coboundary.push_back(15);
  expected_coboundary.push_back(11);
  expected_coboundary.push_back(6);
  expected_coboundary.push_back(16);
  expected_coboundary.push_back(7);
  expected_coboundary.push_back(17);
  expected_coboundary.push_back(11);
  expected_coboundary.push_back(13);
  expected_coboundary.push_back(8);
  expected_coboundary.push_back(18);
  expected_coboundary.push_back(9);
  expected_coboundary.push_back(19);
  expected_coboundary.push_back(13);
  expected_coboundary.push_back(16);
  expected_coboundary.push_back(16);
  expected_coboundary.push_back(18);
  expected_coboundary.push_back(18);
  expected_coboundary.push_back(15);
  expected_coboundary.push_back(21);
  expected_coboundary.push_back(16);
  expected_coboundary.push_back(17);
  expected_coboundary.push_back(21);
  expected_coboundary.push_back(23);
  expected_coboundary.push_back(18);
  expected_coboundary.push_back(19);
  expected_coboundary.push_back(23);

  std::vector<unsigned> sizes(2);
  sizes[0] = 2;
  sizes[1] = 2;

  std::vector<double> data(4);
  data[0] = 0;
  data[1] = 1;
  data[2] = 2;
  data[3] = 3;

  Bitmap_cubical_complex_base ba(sizes, data);
  int i = 0;
  int bd_it = 0;
  int cbd_it = 0;

  Bitmap_cubical_complex_base::All_cells_range range(&ba);
  for (Bitmap_cubical_complex_base::All_cells_iterator it = range.begin(); it != range.end(); ++it) {
    BOOST_CHECK(expected_filtration[i] == ba.get_cell_data(*it));
    BOOST_CHECK(expected_dimension[i] == ba.get_dimension_of_a_cell(*it));

    Bitmap_cubical_complex_base::Boundary_range bdrange = ba.boundary_range(*it);
    for (Bitmap_cubical_complex_base::Boundary_iterator bd = bdrange.begin(); bd != bdrange.end(); ++bd) {
      BOOST_CHECK(expected_boundary[bd_it] == *bd);
      ++bd_it;
    }

    Bitmap_cubical_complex_base::Coboundary_range cbdrange = ba.coboundary_range(*it);
    for (Bitmap_cubical_complex_base::Coboundary_iterator cbd = cbdrange.begin(); cbd != cbdrange.end(); ++cbd) {
      BOOST_CHECK(expected_coboundary[cbd_it] == *cbd);
      ++cbd_it;
    }
    ++i;
  }
}

BOOST_AUTO_TEST_CASE(all_cells_iterator_and_boundary_iterators_in_Bitmap_cubical_complex_base_check_range_check) {
  std::vector<double> expected_filtration;
  expected_filtration.push_back(0);
  expected_filtration.push_back(0);
  expected_filtration.push_back(0);
  expected_filtration.push_back(1);
  expected_filtration.push_back(1);
  expected_filtration.push_back(0);
  expected_filtration.push_back(0);
  expected_filtration.push_back(0);
  expected_filtration.push_back(1);
  expected_filtration.push_back(1);
  expected_filtration.push_back(0);
  expected_filtration.push_back(0);
  expected_filtration.push_back(0);
  expected_filtration.push_back(1);
  expected_filtration.push_back(1);
  expected_filtration.push_back(2);
  expected_filtration.push_back(2);
  expected_filtration.push_back(2);
  expected_filtration.push_back(3);
  expected_filtration.push_back(3);
  expected_filtration.push_back(2);
  expected_filtration.push_back(2);
  expected_filtration.push_back(2);
  expected_filtration.push_back(3);
  expected_filtration.push_back(3);

  std::vector<unsigned> expected_dimension;
  expected_dimension.push_back(0);
  expected_dimension.push_back(1);
  expected_dimension.push_back(0);
  expected_dimension.push_back(1);
  expected_dimension.push_back(0);
  expected_dimension.push_back(1);
  expected_dimension.push_back(2);
  expected_dimension.push_back(1);
  expected_dimension.push_back(2);
  expected_dimension.push_back(1);
  expected_dimension.push_back(0);
  expected_dimension.push_back(1);
  expected_dimension.push_back(0);
  expected_dimension.push_back(1);
  expected_dimension.push_back(0);
  expected_dimension.push_back(1);
  expected_dimension.push_back(2);
  expected_dimension.push_back(1);
  expected_dimension.push_back(2);
  expected_dimension.push_back(1);
  expected_dimension.push_back(0);
  expected_dimension.push_back(1);
  expected_dimension.push_back(0);
  expected_dimension.push_back(1);
  expected_dimension.push_back(0);

  std::vector<size_t> expected_boundary;
  expected_boundary.push_back(0);
  expected_boundary.push_back(2);
  expected_boundary.push_back(2);
  expected_boundary.push_back(4);
  expected_boundary.push_back(0);
  expected_boundary.push_back(10);
  expected_boundary.push_back(1);
  expected_boundary.push_back(11);
  expected_boundary.push_back(7);
  expected_boundary.push_back(5);
  expected_boundary.push_back(2);
  expected_boundary.push_back(12);
  expected_boundary.push_back(3);
  expected_boundary.push_back(13);
  expected_boundary.push_back(9);
  expected_boundary.push_back(7);
  expected_boundary.push_back(4);
  expected_boundary.push_back(14);
  expected_boundary.push_back(10);
  expected_boundary.push_back(12);
  expected_boundary.push_back(12);
  expected_boundary.push_back(14);
  expected_boundary.push_back(10);
  expected_boundary.push_back(20);
  expected_boundary.push_back(11);
  expected_boundary.push_back(21);
  expected_boundary.push_back(17);
  expected_boundary.push_back(15);
  expected_boundary.push_back(12);
  expected_boundary.push_back(22);
  expected_boundary.push_back(13);
  expected_boundary.push_back(23);
  expected_boundary.push_back(19);
  expected_boundary.push_back(17);
  expected_boundary.push_back(14);
  expected_boundary.push_back(24);
  expected_boundary.push_back(20);
  expected_boundary.push_back(22);
  expected_boundary.push_back(22);
  expected_boundary.push_back(24);

  std::vector<size_t> expected_coboundary;
  expected_coboundary.push_back(5);
  expected_coboundary.push_back(1);
  expected_coboundary.push_back(6);
  expected_coboundary.push_back(7);
  expected_coboundary.push_back(1);
  expected_coboundary.push_back(3);
  expected_coboundary.push_back(8);
  expected_coboundary.push_back(9);
  expected_coboundary.push_back(3);
  expected_coboundary.push_back(6);
  expected_coboundary.push_back(6);
  expected_coboundary.push_back(8);
  expected_coboundary.push_back(8);
  expected_coboundary.push_back(5);
  expected_coboundary.push_back(15);
  expected_coboundary.push_back(11);
  expected_coboundary.push_back(6);
  expected_coboundary.push_back(16);
  expected_coboundary.push_back(7);
  expected_coboundary.push_back(17);
  expected_coboundary.push_back(11);
  expected_coboundary.push_back(13);
  expected_coboundary.push_back(8);
  expected_coboundary.push_back(18);
  expected_coboundary.push_back(9);
  expected_coboundary.push_back(19);
  expected_coboundary.push_back(13);
  expected_coboundary.push_back(16);
  expected_coboundary.push_back(16);
  expected_coboundary.push_back(18);
  expected_coboundary.push_back(18);
  expected_coboundary.push_back(15);
  expected_coboundary.push_back(21);
  expected_coboundary.push_back(16);
  expected_coboundary.push_back(17);
  expected_coboundary.push_back(21);
  expected_coboundary.push_back(23);
  expected_coboundary.push_back(18);
  expected_coboundary.push_back(19);
  expected_coboundary.push_back(23);

  std::vector<unsigned> sizes(2);
  sizes[0] = 2;
  sizes[1] = 2;

  std::vector<double> data(4);
  data[0] = 0;
  data[1] = 1;
  data[2] = 2;
  data[3] = 3;

  Bitmap_cubical_complex_base ba(sizes, data);
  int i = 0;
  int bd_it = 0;
  int cbd_it = 0;

  Bitmap_cubical_complex_base::All_cells_range range = ba.all_cells_range();
  for (Bitmap_cubical_complex_base::All_cells_iterator it = range.begin(); it != range.end(); ++it) {
    BOOST_CHECK(expected_filtration[i] == ba.get_cell_data(*it));
    BOOST_CHECK(expected_dimension[i] == ba.get_dimension_of_a_cell(*it));

    Bitmap_cubical_complex_base::Boundary_range bdrange = ba.boundary_range(*it);
    for (Bitmap_cubical_complex_base::Boundary_iterator bd = bdrange.begin(); bd != bdrange.end(); ++bd) {
      BOOST_CHECK(expected_boundary[bd_it] == *bd);
      ++bd_it;
    }

    Bitmap_cubical_complex_base::Coboundary_range cbdrange = ba.coboundary_range(*it);
    for (Bitmap_cubical_complex_base::Coboundary_iterator cbd = cbdrange.begin(); cbd != cbdrange.end(); ++cbd) {
      BOOST_CHECK(expected_coboundary[cbd_it] == *cbd);
      ++cbd_it;
    }
    ++i;
  }
}

BOOST_AUTO_TEST_CASE(Top_dimensional_cells_iterator_range_check) {
  std::vector<double> expected_filtration;
  expected_filtration.push_back(0);
  expected_filtration.push_back(0);
  expected_filtration.push_back(0);
  expected_filtration.push_back(1);
  expected_filtration.push_back(1);
  expected_filtration.push_back(0);
  expected_filtration.push_back(0);
  expected_filtration.push_back(0);
  expected_filtration.push_back(1);
  expected_filtration.push_back(1);
  expected_filtration.push_back(0);
  expected_filtration.push_back(0);
  expected_filtration.push_back(0);
  expected_filtration.push_back(1);
  expected_filtration.push_back(1);
  expected_filtration.push_back(2);
  expected_filtration.push_back(2);
  expected_filtration.push_back(2);
  expected_filtration.push_back(3);
  expected_filtration.push_back(3);
  expected_filtration.push_back(2);
  expected_filtration.push_back(2);
  expected_filtration.push_back(2);
  expected_filtration.push_back(3);
  expected_filtration.push_back(3);

  std::vector<unsigned> sizes(2);
  sizes[0] = 2;
  sizes[1] = 2;

  std::vector<double> data(4);
  data[0] = 0;
  data[1] = 1;
  data[2] = 2;
  data[3] = 3;

  Bitmap_cubical_complex_base ba(sizes, data);
  int i = 0;

  Bitmap_cubical_complex_base::Top_dimensional_cells_range range = ba.top_dimensional_cells_range();
  for (Bitmap_cubical_complex_base::Top_dimensional_cells_iterator it = range.begin(); it != range.end(); ++it) {
    BOOST_CHECK(data[i] == ba.get_cell_data(*it));
    BOOST_CHECK(ba.get_dimension_of_a_cell(*it) == 2);
    ++i;
  }
}

BOOST_AUTO_TEST_CASE(check_if_boundary_of_boundary_is_zero_non_periodic_case_3_d) {
  std::vector<unsigned> sizes(3);
  sizes[0] = 3;
  sizes[1] = 3;
  sizes[2] = 3;

  std::vector<double> data(27, 0);

  int number_of_all_elements = (2 * sizes[0] + 1) * (2 * sizes[1] + 1) * (2 * sizes[2] + 1);
  std::vector<int> elems_in_boundary(number_of_all_elements, 0);
  Bitmap_cubical_complex_base ba(sizes, data);
  for (Bitmap_cubical_complex_base::All_cells_iterator it = ba.all_cells_iterator_begin();
       it != ba.all_cells_iterator_end(); ++it) {    
    int i = 1;
    Bitmap_cubical_complex_base::Boundary_range bdrange = ba.boundary_range(*it);
    for (Bitmap_cubical_complex_base::Boundary_iterator bd = bdrange.begin(); bd != bdrange.end(); ++bd) {
      Bitmap_cubical_complex_base::Boundary_range second_bdrange = ba.boundary_range(*bd);
      int j = 1;
      for (Bitmap_cubical_complex_base::Boundary_iterator bd2 = second_bdrange.begin(); bd2 != second_bdrange.end();
           ++bd2) {
        elems_in_boundary[*bd2] += i * j;
        j *= -1;
      }
      i *= -1;
    }
    // check if there is anything nonzero in elems_in_boundary
    for (size_t i = 0; i != elems_in_boundary.size(); ++i) {
      BOOST_CHECK(elems_in_boundary[i] == 0);
    }
  }
}

BOOST_AUTO_TEST_CASE(check_if_boundary_of_boundary_is_zero_non_periodic_case_4_d) {
  std::vector<unsigned> sizes(4);
  sizes[0] = 3;
  sizes[1] = 3;
  sizes[2] = 3;
  sizes[3] = 3;

  std::vector<double> data(81, 0);

  int number_of_all_elements = (2 * sizes[0] + 1) * (2 * sizes[1] + 1) * (2 * sizes[2] + 1) * (2 * sizes[3] + 1);
  std::vector<int> elems_in_boundary(number_of_all_elements, 0);
  Bitmap_cubical_complex_base ba(sizes, data);
  for (Bitmap_cubical_complex_base::All_cells_iterator it = ba.all_cells_iterator_begin();
       it != ba.all_cells_iterator_end(); ++it) {    
    int i = 1;
    Bitmap_cubical_complex_base::Boundary_range bdrange = ba.boundary_range(*it);
    for (Bitmap_cubical_complex_base::Boundary_iterator bd = bdrange.begin(); bd != bdrange.end(); ++bd) {
      Bitmap_cubical_complex_base::Boundary_range second_bdrange = ba.boundary_range(*bd);
      int j = 1;
      for (Bitmap_cubical_complex_base::Boundary_iterator bd2 = second_bdrange.begin(); bd2 != second_bdrange.end();
           ++bd2) {
        elems_in_boundary[*bd2] += i * j;
        j *= -1;
      }
      i *= -1;
    }
    // check if there is anything nonzero in elems_in_boundary
    for (size_t i = 0; i != elems_in_boundary.size(); ++i) {
      BOOST_CHECK(elems_in_boundary[i] == 0);
    }
  }
}

BOOST_AUTO_TEST_CASE(check_if_boundary_of_boundary_is_zero_periodic_case_2d) {
  std::vector<unsigned> sizes(2);
  sizes[0] = 3;
  sizes[1] = 3;

  std::vector<bool> directions_of_periodicity(2, true);
  std::vector<double> data(9, 0);

  int number_of_all_elements = (2 * sizes[0]) * (2 * sizes[1]);  // *(2*sizes[2]);
  std::vector<int> elems_in_boundary(number_of_all_elements, 0);
  Bitmap_cubical_complex_periodic_boundary_conditions ba(sizes, data, directions_of_periodicity);
  for (Bitmap_cubical_complex_periodic_boundary_conditions::All_cells_iterator it = ba.all_cells_iterator_begin();
       it != ba.all_cells_iterator_end(); ++it) {    
    int i = 1;

    // std::cout << "Element : " << *it << std::endl;

    Bitmap_cubical_complex_periodic_boundary_conditions_base::Boundary_range bdrange = ba.boundary_range(*it);
    for (Bitmap_cubical_complex_periodic_boundary_conditions::Boundary_iterator bd = bdrange.begin();
         bd != bdrange.end(); ++bd) {
      // std::cout <<  *bd << " ";
      Bitmap_cubical_complex_periodic_boundary_conditions::Boundary_range second_bdrange = ba.boundary_range(*bd);
      int j = 1;
      for (Bitmap_cubical_complex_periodic_boundary_conditions::Boundary_iterator bd2 = second_bdrange.begin();
           bd2 != second_bdrange.end(); ++bd2) {
        elems_in_boundary[*bd2] += i * j;
        j *= -1;
      }
      i *= -1;
    }
    // getchar();

    // check if there is anything nonzero in elems_in_boundary
    for (size_t i = 0; i != elems_in_boundary.size(); ++i) {
      BOOST_CHECK(elems_in_boundary[i] == 0);
    }
  }
}

BOOST_AUTO_TEST_CASE(check_if_boundary_of_boundary_is_zero_periodic_case_3d) {
  std::vector<unsigned> sizes(3);
  sizes[0] = 3;
  sizes[1] = 3;
  sizes[2] = 3;

  std::vector<bool> directions_of_periodicity(3, true);

  std::vector<double> data(27, 0);

  int number_of_all_elements = (2 * sizes[0]) * (2 * sizes[1]) * (2 * sizes[2]);
  Bitmap_cubical_complex_periodic_boundary_conditions ba(sizes, data, directions_of_periodicity);
  std::vector<int> elems_in_boundary(number_of_all_elements, 0);
  for (Bitmap_cubical_complex_periodic_boundary_conditions::All_cells_iterator it = ba.all_cells_iterator_begin();
       it != ba.all_cells_iterator_end(); ++it) {
    // std::cout << "Element : " << *it << std::endl;
    
    int i = 1;

    Bitmap_cubical_complex_periodic_boundary_conditions_base::Boundary_range bdrange = ba.boundary_range(*it);
    for (Bitmap_cubical_complex_periodic_boundary_conditions::Boundary_iterator bd = bdrange.begin();
         bd != bdrange.end(); ++bd) {
      Bitmap_cubical_complex_periodic_boundary_conditions::Boundary_range second_bdrange = ba.boundary_range(*bd);
      // std::cout << *bd << " ";
      int j = 1;
      for (Bitmap_cubical_complex_periodic_boundary_conditions::Boundary_iterator bd2 = second_bdrange.begin();
           bd2 != second_bdrange.end(); ++bd2) {
        elems_in_boundary[*bd2] += i * j;
        j *= -1;
      }
      i *= -1;
    }

    // check if there is anything nonzero in elems_in_boundary
    for (size_t i = 0; i != elems_in_boundary.size(); ++i) {
      BOOST_CHECK(elems_in_boundary[i] == 0);
    }
  }
}

BOOST_AUTO_TEST_CASE(check_if_boundary_of_boundary_is_zero_periodic_case_4d) {
  std::vector<unsigned> sizes(4);
  sizes[0] = 3;
  sizes[1] = 3;
  sizes[2] = 3;
  sizes[3] = 3;

  std::vector<bool> directions_of_periodicity(4, true);

  std::vector<double> data(81, 0);

  int number_of_all_elements = (2 * sizes[0]) * (2 * sizes[1]) * (2 * sizes[2]) * (2 * sizes[3]);
  std::vector<int> elems_in_boundary(number_of_all_elements, 0);
  Bitmap_cubical_complex_periodic_boundary_conditions ba(sizes, data, directions_of_periodicity);
  for (Bitmap_cubical_complex_periodic_boundary_conditions::All_cells_iterator it = ba.all_cells_iterator_begin();
       it != ba.all_cells_iterator_end(); ++it) {    
    int i = 1;

    Bitmap_cubical_complex_periodic_boundary_conditions_base::Boundary_range bdrange = ba.boundary_range(*it);
    for (Bitmap_cubical_complex_periodic_boundary_conditions::Boundary_iterator bd = bdrange.begin();
         bd != bdrange.end(); ++bd) {
      Bitmap_cubical_complex_periodic_boundary_conditions::Boundary_range second_bdrange = ba.boundary_range(*bd);
      int j = 1;
      for (Bitmap_cubical_complex_periodic_boundary_conditions::Boundary_iterator bd2 = second_bdrange.begin();
           bd2 != second_bdrange.end(); ++bd2) {
        elems_in_boundary[*bd2] += i * j;
        j *= -1;
      }
      i *= -1;
    }

    // check if there is anything nonzero in elems_in_boundary
    for (size_t i = 0; i != elems_in_boundary.size(); ++i) {
      BOOST_CHECK(elems_in_boundary[i] == 0);
    }
  }
}

BOOST_AUTO_TEST_CASE(compute_incidence_between_cells_test) {
  std::vector<unsigned> sizes(3);
  sizes[0] = 3;
  sizes[1] = 3;
  sizes[2] = 3;

  std::vector<double> data(27, 0);

  int number_of_all_elements = (2 * sizes[0] + 1) * (2 * sizes[1] + 1) * (2 * sizes[1] + 1);
  Bitmap_cubical_complex_base ba(sizes, data);
  std::vector<int> elems_in_boundary(number_of_all_elements, 0);
  for (Bitmap_cubical_complex_base::All_cells_iterator it = ba.all_cells_iterator_begin();
       it != ba.all_cells_iterator_end(); ++it) {    
    Bitmap_cubical_complex_base::Boundary_range bdrange = ba.boundary_range(*it);
    for (Bitmap_cubical_complex_base::Boundary_iterator bd = bdrange.begin(); bd != bdrange.end(); ++bd) {
      Bitmap_cubical_complex_base::Boundary_range second_bdrange = ba.boundary_range(*bd);
      for (Bitmap_cubical_complex_base::Boundary_iterator bd2 = second_bdrange.begin(); bd2 != second_bdrange.end();
           ++bd2) {
        elems_in_boundary[*bd2] +=
            ba.compute_incidence_between_cells(*it, *bd) * ba.compute_incidence_between_cells(*bd, *bd2);
      }
    }
    // check if there is anything nonzero in elems_in_boundary
    for (size_t i = 0; i != elems_in_boundary.size(); ++i) {
      BOOST_CHECK(elems_in_boundary[i] == 0);
    }
  }
}

BOOST_AUTO_TEST_CASE(compute_incidence_between_cells_test_periodic_boundary_conditions) {
  std::vector<unsigned> sizes(3);
  sizes[0] = 3;
  sizes[1] = 3;
  sizes[2] = 3;

  std::vector<bool> directions_of_periodicity(3, true);
  std::vector<double> data(27, 0);

  int number_of_all_elements = (2 * sizes[0]) * (2 * sizes[1]) * (2 * sizes[2]);
  Bitmap_cubical_complex_periodic_boundary_conditions ba(sizes, data, directions_of_periodicity);
  
  std::vector<int> elems_in_boundary(number_of_all_elements, 0);  
  for (Bitmap_cubical_complex_periodic_boundary_conditions::All_cells_iterator it = ba.all_cells_iterator_begin();
       it != ba.all_cells_iterator_end(); ++it) {    
    Bitmap_cubical_complex_periodic_boundary_conditions_base::Boundary_range bdrange = ba.boundary_range(*it);
    for (Bitmap_cubical_complex_periodic_boundary_conditions::Boundary_iterator bd = bdrange.begin();
         bd != bdrange.end(); ++bd) {
      // std::cout <<  *bd << " ";
      Bitmap_cubical_complex_periodic_boundary_conditions::Boundary_range second_bdrange = ba.boundary_range(*bd);
      for (Bitmap_cubical_complex_periodic_boundary_conditions::Boundary_iterator bd2 = second_bdrange.begin();
           bd2 != second_bdrange.end(); ++bd2) {
        elems_in_boundary[*bd2] +=
            ba.compute_incidence_between_cells(*it, *bd) * ba.compute_incidence_between_cells(*bd, *bd2);
      }
    }
    // check if there is anything nonzero in elems_in_boundary
    for (size_t i = 0; i != elems_in_boundary.size(); ++i) {
      BOOST_CHECK(elems_in_boundary[i] == 0);
    }
  }
}

BOOST_AUTO_TEST_CASE(perseus_file_read) {
  Bitmap_cubical_complex increasing("sinusoid.txt");

  auto it = increasing.top_dimensional_cells_iterator_begin();
  double value = increasing.get_cell_data(*it);
  std::cout << "First value of sinusoid.txt is " << value << std::endl;
  BOOST_CHECK(value == 10.);
  // Next value
  ++it;
  value = increasing.get_cell_data(*it);
  std::cout << "Second value of sinusoid.txt is " << value << std::endl;
  BOOST_CHECK(value == std::numeric_limits<double>::infinity());
}
