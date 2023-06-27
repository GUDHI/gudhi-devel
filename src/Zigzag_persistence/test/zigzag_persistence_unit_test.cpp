#include <boost/test/tools/old/interface.hpp>
#include <iostream>
#include <list>
#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "zigzag_persistence"
#include <boost/test/unit_test.hpp>

#include <gudhi/Zigzag_persistence.h>
#include <gudhi/Simplex_tree.h>

using namespace Gudhi;
using namespace boost::unit_test;
using ST = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_wide_indexation>;
using ZP = Gudhi::zigzag_persistence::Zigzag_persistence<ST>;
using Vertex_handle = ST::Vertex_handle;
using Filtration_value = ST::Filtration_value;
using interval_index = ZP::interval_index;
using interval_filtration = ZP::interval_filtration;

// TODO:
// void persistence_diagram(std::ostream& os, Filtration_value shortest_interval = 0.);

struct cmp_intervals_by_length {
	cmp_intervals_by_length() {}
	bool operator()(interval_filtration p, interval_filtration q) {
		if (p.length() != q.length()) {
			return p.length() > q.length();
		}
		if (p.dim() != q.dim()) {
			return p.dim() < q.dim();
		}
		if (p.birth() != q.birth()) {
			return p.birth() < q.birth();
		}
		return p.death() < q.death();
	}
};

BOOST_AUTO_TEST_CASE(constructor) {
	BOOST_CHECK_NO_THROW(ZP zp);
	BOOST_CHECK_NO_THROW(ZP zp(28));
	BOOST_CHECK_NO_THROW(ZP zp(28,2));
	ZP zp;
	BOOST_CHECK(zp.persistence_diagram(0).empty());
}

void test_barcode(ZP& zp, std::vector<interval_filtration>& barcode)
{
	std::stable_sort(barcode.begin(), barcode.end(), cmp_intervals_by_length());
	auto it = barcode.begin();
	for (const auto& interval : zp.persistence_diagram()){
		BOOST_CHECK_EQUAL(interval.dim(), it->dim());
		BOOST_CHECK_EQUAL(interval.birth(), it->birth());
		BOOST_CHECK_EQUAL(interval.death(), it->death());
		++it;
	}
	BOOST_CHECK(it == barcode.end());
}

void test_indices(ZP& zp, std::vector<interval_index>& indices, std::vector<Filtration_value>& indexToFil)
{
	auto it = indices.begin();
	for (const auto& interval : zp.index_persistence_diagram()){
		BOOST_CHECK_EQUAL(interval.dim(), it->dim());
		BOOST_CHECK_EQUAL(interval.birth(), it->birth());
		BOOST_CHECK_EQUAL(interval.death(), it->death());
		auto p = zp.index_to_filtration(interval.birth(), interval.death());
		BOOST_CHECK_EQUAL(p.first, indexToFil[interval.birth()]);
		BOOST_CHECK_EQUAL(p.second, indexToFil[interval.death()]);
		++it;
	}
	BOOST_CHECK(it == indices.end());
}

std::vector<std::vector<Vertex_handle> > get_simplices()
{
	return {
		{0}, 
		{1}, 
		{2}, 
		{0,1},
		{0,2},
		{3},
		{1,2},
		{4},
		{3,4},
		{5},
		{0,1,2},
		{4,5},
		{3,5},
		{3,4,5},
		{0,1,2},			//r
		{3,4,5},			//r
		{1,4},
		{0,1,2},
		{2,4},
		{3,4,5},
		{0,4},
		{0,2,4},
		{1,2,4},
		{0,1,4},
		{3,4,5},			//r
		{3,4},					//r
		{3,5},					//r
		{0,1,2,4}};
}

std::vector<Filtration_value> get_filtration_values()
{
	return {
		0, 0, 0,
		1, 1, 1,
		2, 2, 2,
		3, 3, 3, 3,
		4,
		5,
		6, 6, 6,
		7, 7, 7, 7, 7, 7,
		8,
		9, 9, 9
	};
}

BOOST_AUTO_TEST_CASE(zigzag_persistence_single) {
	ZP zp(28);
	std::vector<interval_index> realIndices;
	std::vector<interval_filtration> realBarcode;
	realIndices.reserve(13);
	realBarcode.reserve(9);

	std::vector<std::vector<Vertex_handle> > simplices = get_simplices();
	std::vector<Filtration_value> filValues = get_filtration_values();

	for (unsigned int i = 0; i < 14; ++i){
		zp.insert_simplex(simplices[i], filValues[i]);
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

	for (unsigned int i = 14; i < 16; ++i){
		zp.remove_simplex(simplices[i], filValues[i]);
	}

	for (unsigned int i = 16; i < 24; ++i){
		zp.insert_simplex(simplices[i], filValues[i]);
	}

	realIndices.emplace_back(0, 5, 16);
	realIndices.emplace_back(1, 14, 17);
	realIndices.emplace_back(1, 15, 19);
	realIndices.emplace_back(1, 20, 21);
	realIndices.emplace_back(1, 18, 22);
	
	realBarcode.emplace_back(0, 1, 6);
	realBarcode.emplace_back(1, 5, 6);
	realBarcode.emplace_back(1, 6, 7);

	for (unsigned int i = 24; i < 27; ++i){
		zp.remove_simplex(simplices[i], filValues[i]);
	}

	realIndices.emplace_back(1, 24, 25);
	realBarcode.emplace_back(1, 8, 9);

	zp.insert_simplex(simplices[27], filValues[27]);

	realIndices.emplace_back(2, 23, 27);
	realBarcode.emplace_back(2, 7, 9);

	test_indices(zp, realIndices, filValues);
	test_barcode(zp, realBarcode);
}

BOOST_AUTO_TEST_CASE(zigzag_persistence_single_max1) {
	ZP zp(28, 1);
	std::vector<interval_index> realIndices;
	std::vector<Filtration_value> indexToFil(28);
	std::vector<interval_filtration> realBarcode;
	realIndices.reserve(5);
	realBarcode.reserve(3);

	std::vector<std::vector<Vertex_handle> > simplices = get_simplices();
	std::vector<Filtration_value> filValues = get_filtration_values();
	unsigned int currIndex = 0;

	for (unsigned int i = 0; i < 14; ++i){
		zp.insert_simplex(simplices[i], filValues[i]);
		if (simplices[i].size() < 3){
			indexToFil[currIndex++] = filValues[i];
		}
	}

	realIndices.emplace_back(0, 1, 3);
	realIndices.emplace_back(0, 2, 4);
	realIndices.emplace_back(0, 7, 8);
	realIndices.emplace_back(0, 9, 10);

	realBarcode.emplace_back(0, 0, 1);
	realBarcode.emplace_back(0, 0, 1);

	for (unsigned int i = 14; i < 16; ++i){
		zp.remove_simplex(simplices[i], filValues[i]);
		if (simplices[i].size() < 3){
			indexToFil[currIndex++] = filValues[i];
		}
	}

	for (unsigned int i = 16; i < 24; ++i){
		zp.insert_simplex(simplices[i], filValues[i]);
		if (simplices[i].size() < 3){
			indexToFil[currIndex++] = filValues[i];
		}
	}

	realIndices.emplace_back(0, 5, 12);
	realBarcode.emplace_back(0, 1, 6);

	for (unsigned int i = 24; i < 27; ++i){
		zp.remove_simplex(simplices[i], filValues[i]);
		if (simplices[i].size() < 3){
			indexToFil[currIndex++] = filValues[i];
		}
	}

	zp.insert_simplex(simplices[27], filValues[27]);

	test_indices(zp, realIndices, indexToFil);
	test_barcode(zp, realBarcode);
}

BOOST_AUTO_TEST_CASE(zigzag_persistence_batch_with_iterators) {
	ZP zp;
	std::vector<interval_index> realIndices;
	std::vector<interval_filtration> realBarcode;
	realIndices.reserve(13);
	realBarcode.reserve(9);

	std::vector<std::vector<Vertex_handle> > simplices = get_simplices();
	std::vector<Filtration_value> filValues = get_filtration_values();

	zp.insert_simplices_contiguously(simplices.begin(), 
									 simplices.begin() + 14, 
									 filValues.begin());

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

	zp.remove_simplices_contiguously(simplices.begin() + 14,
									 simplices.begin() + 16, 
									 filValues.begin() + 14);
	zp.insert_simplices_contiguously(simplices.begin() + 16, 
									 simplices.begin() + 24, 
									 filValues.begin() + 16);

	realIndices.emplace_back(0, 5, 16);
	realIndices.emplace_back(1, 14, 17);
	realIndices.emplace_back(1, 15, 19);
	realIndices.emplace_back(1, 20, 21);
	realIndices.emplace_back(1, 18, 22);
	
	realBarcode.emplace_back(0, 1, 6);
	realBarcode.emplace_back(1, 5, 6);
	realBarcode.emplace_back(1, 6, 7);

	zp.remove_simplices_contiguously(simplices.begin() + 24, 
									 simplices.begin() + 27, 
									 filValues.begin() + 24);

	realIndices.emplace_back(1, 24, 25);
	realBarcode.emplace_back(1, 8, 9);

	zp.insert_simplices_contiguously(simplices.begin() + 27, 
									 simplices.begin() + 28, 
									 filValues.begin() + 27);

	realIndices.emplace_back(2, 23, 27);
	realBarcode.emplace_back(2, 7, 9);

	test_indices(zp, realIndices, filValues);
	test_barcode(zp, realBarcode);
}

BOOST_AUTO_TEST_CASE(zigzag_persistence_batch_with_iterators_max1) {
	ZP zp(28, 1);
	std::vector<interval_index> realIndices;
	std::vector<Filtration_value> indexToFil(28);
	std::vector<interval_filtration> realBarcode;
	realIndices.reserve(5);
	realBarcode.reserve(3);

	std::vector<std::vector<Vertex_handle> > simplices = get_simplices();
	std::vector<Filtration_value> filValues = get_filtration_values();
	unsigned int currIndex = 0;

	for (unsigned int i = 0; i < 28; ++i){
		if (simplices[i].size() < 3){
			indexToFil[currIndex++] = filValues[i];
		}
	}

	zp.insert_simplices_contiguously(simplices.begin(), 
									 simplices.begin() + 14, 
									 filValues.begin());

	realIndices.emplace_back(0, 1, 3);
	realIndices.emplace_back(0, 2, 4);
	realIndices.emplace_back(0, 7, 8);
	realIndices.emplace_back(0, 9, 10);

	realBarcode.emplace_back(0, 0, 1);
	realBarcode.emplace_back(0, 0, 1);

	zp.remove_simplices_contiguously(simplices.begin() + 14,
									 simplices.begin() + 16, 
									 filValues.begin() + 14);
	zp.insert_simplices_contiguously(simplices.begin() + 16, 
									 simplices.begin() + 24, 
									 filValues.begin() + 16);

	realIndices.emplace_back(0, 5, 12);
	realBarcode.emplace_back(0, 1, 6);

	zp.remove_simplices_contiguously(simplices.begin() + 24, 
									 simplices.begin() + 27, 
									 filValues.begin() + 24);
	zp.insert_simplices_contiguously(simplices.begin() + 27, 
									 simplices.begin() + 28, 
									 filValues.begin() + 27);

	test_indices(zp, realIndices, indexToFil);
	test_barcode(zp, realBarcode);
}

BOOST_AUTO_TEST_CASE(zigzag_persistence_batch) {
	ZP zp;
	std::vector<interval_index> realIndices;
	std::vector<interval_filtration> realBarcode;
	realIndices.reserve(13);
	realBarcode.reserve(9);

	std::vector<std::vector<Vertex_handle> > simplices = get_simplices();
	std::vector<Filtration_value> filValues = get_filtration_values();

	std::vector<std::vector<Vertex_handle> > subSimplices(simplices.begin(), simplices.begin() + 14);
	std::vector<Filtration_value> subFilValues(filValues.begin(), filValues.begin() + 14);
	zp.insert_simplices_contiguously(subSimplices, subFilValues);

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

	subSimplices = std::vector<std::vector<Vertex_handle> >(simplices.begin() + 14, simplices.begin() + 16);
	subFilValues = std::vector<Filtration_value>(filValues.begin() + 14, filValues.begin() + 16);
	zp.remove_simplices_contiguously(subSimplices, subFilValues);

	subSimplices = std::vector<std::vector<Vertex_handle> >(simplices.begin() + 16, simplices.begin() + 24);
	subFilValues = std::vector<Filtration_value>(filValues.begin() + 16, filValues.begin() + 24);
	zp.insert_simplices_contiguously(subSimplices, subFilValues);

	realIndices.emplace_back(0, 5, 16);
	realIndices.emplace_back(1, 14, 17);
	realIndices.emplace_back(1, 15, 19);
	realIndices.emplace_back(1, 20, 21);
	realIndices.emplace_back(1, 18, 22);
	
	realBarcode.emplace_back(0, 1, 6);
	realBarcode.emplace_back(1, 5, 6);
	realBarcode.emplace_back(1, 6, 7);

	subSimplices = std::vector<std::vector<Vertex_handle> >(simplices.begin() + 24, simplices.begin() + 27);
	subFilValues = std::vector<Filtration_value>(filValues.begin() + 24, filValues.begin() + 27);
	zp.remove_simplices_contiguously(subSimplices, subFilValues);

	realIndices.emplace_back(1, 24, 25);
	realBarcode.emplace_back(1, 8, 9);

	subSimplices = std::vector<std::vector<Vertex_handle> >(simplices.begin() + 27, simplices.begin() + 28);
	subFilValues = std::vector<Filtration_value>(filValues.begin() + 27, filValues.begin() + 28);
	zp.insert_simplices_contiguously(subSimplices, subFilValues);

	realIndices.emplace_back(2, 23, 27);
	realBarcode.emplace_back(2, 7, 9);

	test_indices(zp, realIndices, filValues);
	test_barcode(zp, realBarcode);
}

BOOST_AUTO_TEST_CASE(zigzag_persistence_batch_max1) {
	ZP zp(28, 1);
	std::vector<interval_index> realIndices;
	std::vector<Filtration_value> indexToFil(28);
	std::vector<interval_filtration> realBarcode;
	realIndices.reserve(5);
	realBarcode.reserve(3);

	std::vector<std::vector<Vertex_handle> > simplices = get_simplices();
	std::vector<Filtration_value> filValues = get_filtration_values();
	unsigned int currIndex = 0;

	for (unsigned int i = 0; i < 28; ++i){
		if (simplices[i].size() < 3){
			indexToFil[currIndex++] = filValues[i];
		}
	}

	std::vector<std::vector<Vertex_handle> > subSimplices(simplices.begin(), simplices.begin() + 14);
	std::vector<Filtration_value> subFilValues(filValues.begin(), filValues.begin() + 14);
	zp.insert_simplices_contiguously(subSimplices, subFilValues);

	realIndices.emplace_back(0, 1, 3);
	realIndices.emplace_back(0, 2, 4);
	realIndices.emplace_back(0, 7, 8);
	realIndices.emplace_back(0, 9, 10);

	realBarcode.emplace_back(0, 0, 1);
	realBarcode.emplace_back(0, 0, 1);

	subSimplices = std::vector<std::vector<Vertex_handle> >(simplices.begin() + 14, simplices.begin() + 16);
	subFilValues = std::vector<Filtration_value>(filValues.begin() + 14, filValues.begin() + 16);
	zp.remove_simplices_contiguously(subSimplices, subFilValues);

	subSimplices = std::vector<std::vector<Vertex_handle> >(simplices.begin() + 16, simplices.begin() + 24);
	subFilValues = std::vector<Filtration_value>(filValues.begin() + 16, filValues.begin() + 24);
	zp.insert_simplices_contiguously(subSimplices, subFilValues);

	realIndices.emplace_back(0, 5, 12);
	realBarcode.emplace_back(0, 1, 6);

	subSimplices = std::vector<std::vector<Vertex_handle> >(simplices.begin() + 24, simplices.begin() + 27);
	subFilValues = std::vector<Filtration_value>(filValues.begin() + 24, filValues.begin() + 27);
	zp.remove_simplices_contiguously(subSimplices, subFilValues);

	subSimplices = std::vector<std::vector<Vertex_handle> >(simplices.begin() + 27, simplices.begin() + 28);
	subFilValues = std::vector<Filtration_value>(filValues.begin() + 27, filValues.begin() + 28);
	zp.insert_simplices_contiguously(subSimplices, subFilValues);

	test_indices(zp, realIndices, indexToFil);
	test_barcode(zp, realBarcode);
}
