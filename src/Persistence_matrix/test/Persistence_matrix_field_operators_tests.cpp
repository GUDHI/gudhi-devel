/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "persistence_matrix"
#include <boost/test/unit_test.hpp>

#include "gudhi/Fields/Z2_field_operators.h"
#include "gudhi/Fields/Zp_field_operators.h"
#include "gudhi/Fields/Multi_field_operators.h"
#include "gudhi/Fields/Multi_field_small_operators.h"

using Gudhi::persistence_matrix::Z2_field_operators;
using Gudhi::persistence_matrix::Zp_field_operators;
using Gudhi::persistence_matrix::Multi_field_operators;
using Gudhi::persistence_matrix::Multi_field_operators_with_small_characteristics;

template<class Z2>
void test_z2_standart_field_operators(Z2& op)
{
	using T = typename Z2::element_type;

	unsigned int z21 = 7u;
	unsigned int z22 = 2u;
	T z23 = op.get_value(7u);
	T z24 = op.get_value(2u);

	//e * m + a
	BOOST_CHECK_EQUAL(op.multiply_and_add(z21, z22, 3u), 1);
	BOOST_CHECK_EQUAL(op.multiply_and_add(z23, z24, T(3)), 1);

	//(e + a) * m
	BOOST_CHECK_EQUAL(op.add_and_multiply(z21, z22, 3u), 1);
	BOOST_CHECK_EQUAL(op.add_and_multiply(z23, z24, T(3)), 1);

	//+
	BOOST_CHECK_EQUAL(op.add(z21, z22), 1);
	BOOST_CHECK_EQUAL(op.add(z23, z24), 1);
	BOOST_CHECK_EQUAL(op.add(z23, T(3)), 0);
	BOOST_CHECK_EQUAL(op.add(z24, T(3)), 1);
	BOOST_CHECK_EQUAL(op.add(op.get_value(6u), z23), 1);
	BOOST_CHECK_EQUAL(op.add(op.get_value(6u), z24), 0);
	z21 = op.add(z21, 3u);
	BOOST_CHECK_EQUAL(z21, 0);
	z21 = op.add(z21, z22);
	BOOST_CHECK_EQUAL(z21, 0);
	z23 = op.get_value(z21);

	//-
	BOOST_CHECK_EQUAL(op.substract(z21, z22), 0);
	BOOST_CHECK_EQUAL(op.substract(z23, z24), 0);
	BOOST_CHECK_EQUAL(op.substract(z23, T(3)), 1);
	BOOST_CHECK_EQUAL(op.substract(z24, T(3)), 1);
	BOOST_CHECK_EQUAL(op.substract(op.get_value(6u), z23), 0);
	BOOST_CHECK_EQUAL(op.substract(op.get_value(6u), z24), 0);
	z21 = op.substract(z21, 3u);
	BOOST_CHECK_EQUAL(z21, 1);
	z21 = op.substract(z21, z22);
	BOOST_CHECK_EQUAL(z21, 1);
	z23 = op.get_value(z21);

	//*
	BOOST_CHECK_EQUAL(op.multiply(z21, z22), 0);
	BOOST_CHECK_EQUAL(op.multiply(z23, z24), 0);
	BOOST_CHECK_EQUAL(op.multiply(z23, T(3)), 1);
	BOOST_CHECK_EQUAL(op.multiply(z24, T(3)), 0);
	BOOST_CHECK_EQUAL(op.multiply(op.get_value(6u), z23), 0);
	BOOST_CHECK_EQUAL(op.multiply(op.get_value(6u), z24), 0);
	z21 = op.multiply(z21, 3u);
	BOOST_CHECK_EQUAL(z21, 1);
	z21 = op.multiply(z21, z22);
	BOOST_CHECK_EQUAL(z21, 0);
	z23 = op.get_value(z21);

	//==
	BOOST_CHECK(op.are_equal(z21, 0u));
	BOOST_CHECK(op.are_equal(z21, z22));
	BOOST_CHECK(op.are_equal(op.get_value(z21), z23));
	BOOST_CHECK(op.are_equal(op.get_value(z21), z24));
}

template<class Z5>
void test_z5_standart_field_operators(Z5& op){
	unsigned int z51 = 7;
	unsigned int z52 = 3;
	auto z53 = op.get_value(7);
	auto z54 = op.get_value(3);

	//e * m + a
	BOOST_CHECK_EQUAL(op.multiply_and_add(z51, z52, 3), 4);
	BOOST_CHECK_EQUAL(op.multiply_and_add(z53, z54, 3), 4);

	//(e + a) * m
	BOOST_CHECK_EQUAL(op.add_and_multiply(z51, z52, 3), 0);
	BOOST_CHECK_EQUAL(op.add_and_multiply(z53, z54, 3), 0);

	//+
	BOOST_CHECK_EQUAL(op.add(z51, z52), 0);
	BOOST_CHECK_EQUAL(op.add(z53, z54), 0);
	BOOST_CHECK_EQUAL(op.add(z53, 3), 0);
	BOOST_CHECK_EQUAL(op.add(z54, 3), 1);
	BOOST_CHECK_EQUAL(op.add(7, z53), 4);
	BOOST_CHECK_EQUAL(op.add(7, z54), 0);
	z51 = op.add(z51, 3);
	BOOST_CHECK_EQUAL(z51, 0);
	z51 = op.add(z51, z52);
	BOOST_CHECK_EQUAL(z51, 3);
	z53 = op.get_value(z51);

	//-
	BOOST_CHECK_EQUAL(op.substract(z51, z52), 0);
	BOOST_CHECK_EQUAL(op.substract(z53, z54), 0);
	BOOST_CHECK_EQUAL(op.substract(z53, 3), 0);
	BOOST_CHECK_EQUAL(op.substract(z54, 3), 0);
	BOOST_CHECK_EQUAL(op.substract(7, z53), 4);
	BOOST_CHECK_EQUAL(op.substract(7, z54), 4);
	z51 = op.substract(z51, 3);
	BOOST_CHECK_EQUAL(z51, 0);
	z51 = op.substract(z51, z52);
	BOOST_CHECK_EQUAL(z51, 2);
	z53 = op.get_value(z51);

	//*
	BOOST_CHECK_EQUAL(op.multiply(z51, z52), 1);
	BOOST_CHECK_EQUAL(op.multiply(z53, z54), 1);
	BOOST_CHECK_EQUAL(op.multiply(z53, 3), 1);
	BOOST_CHECK_EQUAL(op.multiply(z54, 3), 4);
	BOOST_CHECK_EQUAL(op.multiply(7, z53), 4);
	BOOST_CHECK_EQUAL(op.multiply(7, z54), 1);
	z51 = op.multiply(z51, 3);
	BOOST_CHECK_EQUAL(z51, 1);
	z51 = op.multiply(z51, z52);
	BOOST_CHECK_EQUAL(z51, 3);
	z53 = op.get_value(z51);

	//==
	BOOST_CHECK(op.are_equal(z51, 3));
	BOOST_CHECK(op.are_equal(z51, z52));
	BOOST_CHECK(op.are_equal(z51, z53));
	BOOST_CHECK(op.are_equal(z51, z54));
}

template<class Z2>
void test_z2_standart_field_properties(Z2& op){
	unsigned int z21 = 7;
	unsigned int z22 = 2;

	BOOST_CHECK_EQUAL(op.get_inverse(z21), 1);
	BOOST_CHECK_EQUAL(op.get_inverse(z22), 0);
	BOOST_CHECK(op.get_partial_inverse(z21, 35) == std::make_pair(typename Z2::element_type(1), 35u));
	BOOST_CHECK(op.get_partial_inverse(z22, 35) == std::make_pair(typename Z2::element_type(0), 35u));

	BOOST_CHECK_EQUAL(op.get_additive_identity(), 0);
	BOOST_CHECK_EQUAL(op.get_multiplicative_identity(), 1);
	BOOST_CHECK_EQUAL(op.get_partial_multiplicative_identity(35), 1);

	BOOST_CHECK_EQUAL(op.get_characteristic(), 2);

	BOOST_CHECK_EQUAL(op.get_value(z21), 1);
	BOOST_CHECK_EQUAL(op.get_value(z22), 0);
}

template<class Z5>
void test_z5_standart_field_properties(Z5& op){
	unsigned int z51 = 7;
	unsigned int z52 = 3;

	BOOST_CHECK_EQUAL(op.get_inverse(z51), 3);
	BOOST_CHECK_EQUAL(op.get_inverse(z52), 2);
	BOOST_CHECK(op.get_partial_inverse(z51, 35) == std::make_pair(3u, 35u));
	BOOST_CHECK(op.get_partial_inverse(z52, 35) == std::make_pair(2u, 35u));

	BOOST_CHECK_EQUAL(op.get_additive_identity(), 0);

	BOOST_CHECK_EQUAL(op.get_multiplicative_identity(), 1);
	BOOST_CHECK_EQUAL(op.get_partial_multiplicative_identity(35), 1);

	BOOST_CHECK_EQUAL(op.get_characteristic(), 5);

	BOOST_CHECK_EQUAL(op.get_value(z51), 2);
	BOOST_CHECK_EQUAL(op.get_value(z52), 3);
}

template<class Z7>
void test_z7_standart_field_properties(Z7& op){
	unsigned int z71 = 8;
	unsigned int z72 = 3;

	BOOST_CHECK_EQUAL(op.get_inverse(z71), 1);
	BOOST_CHECK_EQUAL(op.get_inverse(z72), 5);
	BOOST_CHECK(op.get_partial_inverse(z71, 35) == std::make_pair(1u, 35u));
	BOOST_CHECK(op.get_partial_inverse(z72, 35) == std::make_pair(5u, 35u));

	BOOST_CHECK_EQUAL(op.get_additive_identity(), 0);

	BOOST_CHECK_EQUAL(op.get_multiplicative_identity(), 1);
	BOOST_CHECK_EQUAL(op.get_partial_multiplicative_identity(35), 1);

	BOOST_CHECK_EQUAL(op.get_characteristic(), 7);

	BOOST_CHECK_EQUAL(op.get_value(z71), 1);
	BOOST_CHECK_EQUAL(op.get_value(z72), 3);
}

BOOST_AUTO_TEST_CASE(Field_operators_operation)
{
	Z2_field_operators z2op;
	test_z2_standart_field_operators(z2op);

	Zp_field_operators zpop;
	zpop.set_characteristic(2);
	test_z2_standart_field_operators(zpop);
	zpop.set_characteristic(5);
	test_z5_standart_field_operators(zpop);
}

BOOST_AUTO_TEST_CASE(Field_operators_properties)
{
	Z2_field_operators z2op;
	test_z2_standart_field_properties(z2op);

	Zp_field_operators zpop;
	zpop.set_characteristic(2);
	test_z2_standart_field_properties(zpop);
	zpop.set_characteristic(5);
	test_z5_standart_field_properties(zpop);
	zpop.set_characteristic(7);
	test_z7_standart_field_properties(zpop);
}

template<class MF>
void test_multi_field_operators(MF& op){
	using T = typename MF::element_type;

	T m1(5005);
	T m2(5007);

	//e * m + a
	BOOST_CHECK_EQUAL(op.multiply_and_add(m1, m2, 3), 3);

	//(e + a) * m
	BOOST_CHECK_EQUAL(op.add_and_multiply(m1, m2, 3), 6);

	//+
	BOOST_CHECK_EQUAL(op.add(m1, m2), T(2));
	BOOST_CHECK_EQUAL(op.add(m1, 3), T(3));
	BOOST_CHECK_EQUAL(op.add(m2, 3), T(5));
	BOOST_CHECK_EQUAL(op.add(6, m1), T(6));
	BOOST_CHECK_EQUAL(op.add(6, m2), T(8));
	m1 = op.add(m1, 3);
	BOOST_CHECK_EQUAL(m1, T(3));
	m1 = op.add(m1, m2);
	BOOST_CHECK_EQUAL(m1, T(5));

	//-
	BOOST_CHECK_EQUAL(op.substract(m1, m2), T(3));
	BOOST_CHECK_EQUAL(op.substract(m1, 3), T(2));
	BOOST_CHECK_EQUAL(op.substract(m2, 3), T(5004));
	BOOST_CHECK_EQUAL(op.substract(6, m1), T(1));
	BOOST_CHECK_EQUAL(op.substract(6, m2), T(4));
	m2 = op.substract(m2, 3);
	BOOST_CHECK_EQUAL(m2, T(5004));
	m2 = op.substract(m2, m1);
	BOOST_CHECK_EQUAL(m2, T(4999));

	//*
	BOOST_CHECK_EQUAL(op.multiply(m1, m2), T(4975));
	BOOST_CHECK_EQUAL(op.multiply(m1, 3), T(15));
	BOOST_CHECK_EQUAL(op.multiply(m2, 3), T(4987));
	BOOST_CHECK_EQUAL(op.multiply(6, m1), T(30));
	BOOST_CHECK_EQUAL(op.multiply(6, m2), T(4969));
	m1 = op.multiply(m1, 3);
	BOOST_CHECK_EQUAL(m1, T(15));
	m1 = op.multiply(m1, m2);
	BOOST_CHECK_EQUAL(m1, T(4915));

	//==
	BOOST_CHECK(!op.are_equal(m1, m2));
	m2 = op.substract(m2, 84);
	BOOST_CHECK(op.are_equal(m1, T(4915)));
	BOOST_CHECK(op.are_equal(m1, m2));
}

template<class MF>
void test_multi_field_properties(MF& op){
	using T = typename MF::element_type;

	T m1(1);
	T m2(7);

	BOOST_CHECK_EQUAL(op.get_inverse(m1), T(1));
	BOOST_CHECK_EQUAL(op.get_inverse(m2), T(2758));
	BOOST_CHECK(op.get_partial_inverse(m1, 35) == std::make_pair(T(1716), T(35)));
	BOOST_CHECK(op.get_partial_inverse(m2, 35) == std::make_pair(T(3003), T(5)));

	BOOST_CHECK_EQUAL(op.get_additive_identity(), T(0));
	BOOST_CHECK_EQUAL(op.get_multiplicative_identity(), T(1));
	BOOST_CHECK_EQUAL(op.get_partial_multiplicative_identity(7), T(715));

	BOOST_CHECK_EQUAL(op.get_characteristic(), 5005);

	BOOST_CHECK_EQUAL(op.get_value(m1), 1);
	BOOST_CHECK_EQUAL(op.get_value(m2), 7);
}

BOOST_AUTO_TEST_CASE(Multi_Field_operators_operation)
{
	Multi_field_operators mfop;
	mfop.set_characteristic(5, 13);
	test_multi_field_operators(mfop);

	Multi_field_operators_with_small_characteristics smfop;
	smfop.set_characteristic(5, 13);
	test_multi_field_operators(smfop);
}

BOOST_AUTO_TEST_CASE(Multi_Field_operators_properties)
{
	Multi_field_operators mfop;
	mfop.set_characteristic(5, 13);
	test_multi_field_properties(mfop);

	Multi_field_operators_with_small_characteristics smfop;
	smfop.set_characteristic(5, 13);
	test_multi_field_properties(smfop);

	mfop.set_characteristic(3, 30);
	smfop.set_characteristic(3, 30);

	BOOST_CHECK_EQUAL(mfop.get_characteristic(), smfop.get_characteristic());	// == 3234846615
	BOOST_CHECK_EQUAL(mfop.get_partial_inverse(2, 35).first, smfop.get_partial_inverse(2, 35).first);	// == 2033332158
	BOOST_CHECK_EQUAL(mfop.get_partial_inverse(2, 35).second, smfop.get_partial_inverse(2, 35).second);	// == 2033332158
	BOOST_CHECK_EQUAL(mfop.get_partial_multiplicative_identity(35), smfop.get_partial_multiplicative_identity(35));
}


