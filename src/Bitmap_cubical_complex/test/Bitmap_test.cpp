#include "gudhi/Bitmap_cubical_complex.h"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "cubical_complex"
#include <boost/test/unit_test.hpp>

using namespace std;

BOOST_AUTO_TEST_CASE(check_dimension) {
  std::vector< double > increasingFiltrationOfTopDimensionalCells;
  increasingFiltrationOfTopDimensionalCells.push_back(1);
  increasingFiltrationOfTopDimensionalCells.push_back(2);
  increasingFiltrationOfTopDimensionalCells.push_back(3);
  increasingFiltrationOfTopDimensionalCells.push_back(4);
  increasingFiltrationOfTopDimensionalCells.push_back(5);
  increasingFiltrationOfTopDimensionalCells.push_back(6);
  increasingFiltrationOfTopDimensionalCells.push_back(7);
  increasingFiltrationOfTopDimensionalCells.push_back(8);
  increasingFiltrationOfTopDimensionalCells.push_back(9);

  std::vector<unsigned> dimensions;
  dimensions.push_back(3);
  dimensions.push_back(3);

  Bitmap_cubical_complex<double> increasing(dimensions, increasingFiltrationOfTopDimensionalCells);
  BOOST_CHECK(increasing.dimension() == 2);
}

BOOST_AUTO_TEST_CASE(topDimensionalCellsIterator_test) {
  std::vector< double > expectedFiltrationValues1;
  expectedFiltrationValues1.push_back(0);
  expectedFiltrationValues1.push_back(0);
  expectedFiltrationValues1.push_back(0);
  expectedFiltrationValues1.push_back(0);
  expectedFiltrationValues1.push_back(100);
  expectedFiltrationValues1.push_back(0);
  expectedFiltrationValues1.push_back(0);
  expectedFiltrationValues1.push_back(0);
  expectedFiltrationValues1.push_back(0);

  std::vector< double > expectedFiltrationValues2;
  expectedFiltrationValues2.push_back(1);
  expectedFiltrationValues2.push_back(2);
  expectedFiltrationValues2.push_back(3);
  expectedFiltrationValues2.push_back(4);
  expectedFiltrationValues2.push_back(5);
  expectedFiltrationValues2.push_back(6);
  expectedFiltrationValues2.push_back(7);
  expectedFiltrationValues2.push_back(8);
  expectedFiltrationValues2.push_back(9);

  std::vector< double > increasingFiltrationOfTopDimensionalCells;
  increasingFiltrationOfTopDimensionalCells.push_back(1);
  increasingFiltrationOfTopDimensionalCells.push_back(2);
  increasingFiltrationOfTopDimensionalCells.push_back(3);
  increasingFiltrationOfTopDimensionalCells.push_back(4);
  increasingFiltrationOfTopDimensionalCells.push_back(5);
  increasingFiltrationOfTopDimensionalCells.push_back(6);
  increasingFiltrationOfTopDimensionalCells.push_back(7);
  increasingFiltrationOfTopDimensionalCells.push_back(8);
  increasingFiltrationOfTopDimensionalCells.push_back(9);

  std::vector< double > oneDimensionalCycle;
  oneDimensionalCycle.push_back(0);
  oneDimensionalCycle.push_back(0);
  oneDimensionalCycle.push_back(0);
  oneDimensionalCycle.push_back(0);
  oneDimensionalCycle.push_back(100);
  oneDimensionalCycle.push_back(0);
  oneDimensionalCycle.push_back(0);
  oneDimensionalCycle.push_back(0);
  oneDimensionalCycle.push_back(0);

  std::vector<unsigned> dimensions;
  dimensions.push_back(3);
  dimensions.push_back(3);

  Bitmap_cubical_complex<double> increasing(dimensions, increasingFiltrationOfTopDimensionalCells);
  Bitmap_cubical_complex<double> hole(dimensions, oneDimensionalCycle);


  int i = 0;
  for (Bitmap_cubical_complex<double>::Top_dimensional_cells_iterator it = increasing.top_dimensional_cells_begin(); it != increasing.top_dimensional_cells_end(); ++it) {
    BOOST_CHECK(*it == expectedFiltrationValues2[i]);
    ++i;
  }
  i = 0;
  for (Bitmap_cubical_complex<double>::Top_dimensional_cells_iterator it = hole.top_dimensional_cells_begin(); it != hole.top_dimensional_cells_end(); ++it) {
    BOOST_CHECK(*it == expectedFiltrationValues1[i]);
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
  boundary8.push_back(7);
  boundary8.push_back(9);
  std::vector<double> boundary9;
  boundary9.push_back(2);
  boundary9.push_back(16);
  std::vector<double> boundary10;
  boundary10.push_back(3);
  boundary10.push_back(17);
  boundary10.push_back(9);
  boundary10.push_back(11);
  std::vector<double> boundary11;
  boundary11.push_back(4);
  boundary11.push_back(18);
  std::vector<double> boundary12;
  boundary12.push_back(5);
  boundary12.push_back(19);
  boundary12.push_back(11);
  boundary12.push_back(13);
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
  boundary22.push_back(21);
  boundary22.push_back(23);
  std::vector<double> boundary23;
  boundary23.push_back(16);
  boundary23.push_back(30);
  std::vector<double> boundary24;
  boundary24.push_back(17);
  boundary24.push_back(31);
  boundary24.push_back(23);
  boundary24.push_back(25);
  std::vector<double> boundary25;
  boundary25.push_back(18);
  boundary25.push_back(32);
  std::vector<double> boundary26;
  boundary26.push_back(19);
  boundary26.push_back(33);
  boundary26.push_back(25);
  boundary26.push_back(27);
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
  boundary36.push_back(35);
  boundary36.push_back(37);
  std::vector<double> boundary37;
  boundary37.push_back(30);
  boundary37.push_back(44);
  std::vector<double> boundary38;
  boundary38.push_back(31);
  boundary38.push_back(45);
  boundary38.push_back(37);
  boundary38.push_back(39);
  std::vector<double> boundary39;
  boundary39.push_back(32);
  boundary39.push_back(46);
  std::vector<double> boundary40;
  boundary40.push_back(33);
  boundary40.push_back(47);
  boundary40.push_back(39);
  boundary40.push_back(41);
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
  std::vector< std::vector<double> > boundaries;
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



  std::vector< double > increasingFiltrationOfTopDimensionalCells;
  increasingFiltrationOfTopDimensionalCells.push_back(1);
  increasingFiltrationOfTopDimensionalCells.push_back(2);
  increasingFiltrationOfTopDimensionalCells.push_back(3);
  increasingFiltrationOfTopDimensionalCells.push_back(4);
  increasingFiltrationOfTopDimensionalCells.push_back(5);
  increasingFiltrationOfTopDimensionalCells.push_back(6);
  increasingFiltrationOfTopDimensionalCells.push_back(7);
  increasingFiltrationOfTopDimensionalCells.push_back(8);
  increasingFiltrationOfTopDimensionalCells.push_back(9);

  std::vector<unsigned> dimensions;
  dimensions.push_back(3);
  dimensions.push_back(3);

  Bitmap_cubical_complex<double> increasing(dimensions, increasingFiltrationOfTopDimensionalCells);
  for (size_t i = 0; i != increasing.size_of_bitmap(); ++i) {
    std::vector< size_t > bd = increasing.get_boundary_of_a_cell(i);
    for (size_t j = 0; j != bd.size(); ++j) {
      BOOST_CHECK(boundaries[i][j] == bd[j]);
    }
  }
}

BOOST_AUTO_TEST_CASE(compute_boundary_test_2) {
  std::vector< double > increasingFiltrationOfTopDimensionalCells;
  increasingFiltrationOfTopDimensionalCells.push_back(1);
  increasingFiltrationOfTopDimensionalCells.push_back(2);
  increasingFiltrationOfTopDimensionalCells.push_back(3);
  increasingFiltrationOfTopDimensionalCells.push_back(4);
  increasingFiltrationOfTopDimensionalCells.push_back(5);
  increasingFiltrationOfTopDimensionalCells.push_back(6);
  increasingFiltrationOfTopDimensionalCells.push_back(7);
  increasingFiltrationOfTopDimensionalCells.push_back(8);
  increasingFiltrationOfTopDimensionalCells.push_back(9);

  std::vector<unsigned> dimensions;
  dimensions.push_back(3);
  dimensions.push_back(3);

  Bitmap_cubical_complex<double> increasing(dimensions, increasingFiltrationOfTopDimensionalCells);


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
  for (size_t i = 0; i != increasing.size_of_bitmap(); ++i) {
    std::vector< size_t > bd = increasing.get_coboundary_of_a_cell(i);
    for (size_t j = 0; j != bd.size(); ++j) {
      BOOST_CHECK(coboundaryElements[number] == bd[j]);
      ++number;
    }

  }
}

BOOST_AUTO_TEST_CASE(compute_boundary_test_3) {
  std::vector< double > increasingFiltrationOfTopDimensionalCells;
  increasingFiltrationOfTopDimensionalCells.push_back(1);
  increasingFiltrationOfTopDimensionalCells.push_back(2);
  increasingFiltrationOfTopDimensionalCells.push_back(3);
  increasingFiltrationOfTopDimensionalCells.push_back(4);
  increasingFiltrationOfTopDimensionalCells.push_back(5);
  increasingFiltrationOfTopDimensionalCells.push_back(6);
  increasingFiltrationOfTopDimensionalCells.push_back(7);
  increasingFiltrationOfTopDimensionalCells.push_back(8);
  increasingFiltrationOfTopDimensionalCells.push_back(9);

  std::vector<unsigned> dimensions;
  dimensions.push_back(3);
  dimensions.push_back(3);

  Bitmap_cubical_complex<double> increasing(dimensions, increasingFiltrationOfTopDimensionalCells);

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

  for (size_t i = 0; i != increasing.size_of_bitmap(); ++i) {
    BOOST_CHECK(increasing.get_dimension_of_a_cell(i) == dim[i]);
  }
}

BOOST_AUTO_TEST_CASE(Filtration_simplex_iterator_test) {
  std::vector< double > increasingFiltrationOfTopDimensionalCells;
  increasingFiltrationOfTopDimensionalCells.push_back(1);
  increasingFiltrationOfTopDimensionalCells.push_back(2);
  increasingFiltrationOfTopDimensionalCells.push_back(3);
  increasingFiltrationOfTopDimensionalCells.push_back(4);
  increasingFiltrationOfTopDimensionalCells.push_back(5);
  increasingFiltrationOfTopDimensionalCells.push_back(6);
  increasingFiltrationOfTopDimensionalCells.push_back(7);
  increasingFiltrationOfTopDimensionalCells.push_back(8);
  increasingFiltrationOfTopDimensionalCells.push_back(9);

  std::vector<unsigned> dimensions;
  dimensions.push_back(3);
  dimensions.push_back(3);

  Bitmap_cubical_complex<double> increasing(dimensions, increasingFiltrationOfTopDimensionalCells);

  std::vector< unsigned > dim;
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


  Bitmap_cubical_complex<double>::Filtration_simplex_range range = increasing.filtration_simplex_range();
  size_t position = 0;
  for (Bitmap_cubical_complex<double>::Filtration_simplex_iterator it = range.begin(); it != range.end(); ++it) {
    BOOST_CHECK(increasing.dimension(*it) == dim[position]);
    BOOST_CHECK(increasing.filtration(*it) == fil[position]);
    ++position;
  }
}
