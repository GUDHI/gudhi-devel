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

#include "gudhi/utilities/Z2_field.h"
#include "gudhi/utilities/Zp_field.h"

using Gudhi::persistence_matrix::Z2_field_element;
using Gudhi::persistence_matrix::Zp_field_element;

BOOST_AUTO_TEST_CASE(Field_constructors)
{
	//default constructor
	Z2_field_element z2_d;
	BOOST_CHECK_EQUAL(z2_d, 0u);
	Zp_field_element<5> z5_d;
	BOOST_CHECK_EQUAL(z5_d, 0u);
	Zp_field_element<13> z13_d;
	BOOST_CHECK_EQUAL(z13_d, 0u);

	//value constructor
	Z2_field_element z2_v(5);
	BOOST_CHECK_EQUAL(z2_v, 1u);
	Zp_field_element<5> z5_v(7);
	BOOST_CHECK_EQUAL(z5_v, 2u);
	Zp_field_element<13> z13_v(5);
	BOOST_CHECK_EQUAL(z13_v, 5u);

	//copy constructor
	Z2_field_element z2_c1(5);
	Z2_field_element z2_c2 = z2_c1;
	BOOST_CHECK_EQUAL(z2_c2, 1u);
	Z2_field_element z2_c3(z2_c2);
	BOOST_CHECK_EQUAL(z2_c3, 1u);
	Zp_field_element<5> z5_c1(7);
	Zp_field_element<5> z5_c2 = z5_c1;
	BOOST_CHECK_EQUAL(z5_c2, 2u);
	Zp_field_element<5> z5_c3(z5_c2);
	BOOST_CHECK_EQUAL(z5_c3, 2u);
	Zp_field_element<13> z13_c1(5);
	Zp_field_element<13> z13_c2 = z13_c1;
	BOOST_CHECK_EQUAL(z13_c2, 5u);
	Zp_field_element<13> z13_c3(z13_c2);
	BOOST_CHECK_EQUAL(z13_c3, 5u);

	//move constructor
	Z2_field_element z2_m1(5);
	Z2_field_element z2_m2(std::move(z2_m1));
	BOOST_CHECK_EQUAL(z2_m2, 1u);
	BOOST_CHECK_EQUAL(z2_m1, 0u);
	Zp_field_element<5> z5_m1(7);
	Zp_field_element<5> z5_m2(std::move(z5_m1));
	BOOST_CHECK_EQUAL(z5_m2, 2u);
	BOOST_CHECK_EQUAL(z5_m1, 0u);
	Zp_field_element<5> z13_m1(5);
	Zp_field_element<5> z13_m2(std::move(z13_m1));
	BOOST_CHECK_EQUAL(z13_m2, 5u);
	BOOST_CHECK_EQUAL(z13_m1, 0u);

	//swap
	Z2_field_element z2_s1(5);
	Z2_field_element z2_s2(8);
	swap(z2_s1, z2_s2);
	BOOST_CHECK_EQUAL(z2_s2, 1u);
	BOOST_CHECK_EQUAL(z2_s1, 0u);
	Zp_field_element<5> z5_s1(4);
	Zp_field_element<5> z5_s2(8);
	swap(z5_s1, z5_s2);
	BOOST_CHECK_EQUAL(z5_s2, 4u);
	BOOST_CHECK_EQUAL(z5_s1, 3u);
	Zp_field_element<13> z13_s1(4);
	Zp_field_element<13> z13_s2(22);
	swap(z13_s1, z13_s2);
	BOOST_CHECK_EQUAL(z13_s2, 4u);
	BOOST_CHECK_EQUAL(z13_s1, 9u);
}

BOOST_AUTO_TEST_CASE(Field_operators)
{
	Z2_field_element z21(7);
	Z2_field_element z22(2);
	Zp_field_element<5> z51(7);
	Zp_field_element<5> z52(3);

	//+
	BOOST_CHECK_EQUAL(z21 + z22, 1u);
	BOOST_CHECK_EQUAL(z21 + 3u, 0u);
	BOOST_CHECK_EQUAL(z22 + 3u, 1u);
	BOOST_CHECK_EQUAL(6u + z21, 1u);
	BOOST_CHECK_EQUAL(6u + z22, 0u);
	BOOST_CHECK_EQUAL(z51 + z52, 0u);
	BOOST_CHECK_EQUAL(z51 + 3u, 0u);
	BOOST_CHECK_EQUAL(z52 + 3u, 1u);
	BOOST_CHECK_EQUAL(7u + z51, 4u);
	BOOST_CHECK_EQUAL(7u + z52, 0u);
	z21 += 3u;
	BOOST_CHECK_EQUAL(z21, 0u);
	z21 += z22;
	BOOST_CHECK_EQUAL(z21, 0u);
	z51 += 3u;
	BOOST_CHECK_EQUAL(z51, 0u);
	z51 += z52;
	BOOST_CHECK_EQUAL(z51, 3u);

	//-
	BOOST_CHECK_EQUAL(z21 - z22, 0u);
	BOOST_CHECK_EQUAL(z21 - 3u, 1u);
	BOOST_CHECK_EQUAL(z22 - 3u, 1u);
	BOOST_CHECK_EQUAL(6u - z21, 0u);
	BOOST_CHECK_EQUAL(6u - z22, 0u);
	BOOST_CHECK_EQUAL(z51 - z52, 0u);
	BOOST_CHECK_EQUAL(z51 - 3u, 0u);
	BOOST_CHECK_EQUAL(z52 - 3u, 0u);
	BOOST_CHECK_EQUAL(7u - z51, 4u);
	BOOST_CHECK_EQUAL(7u - z52, 4u);
	z21 -= 3u;
	BOOST_CHECK_EQUAL(z21, 1u);
	z21 -= z22;
	BOOST_CHECK_EQUAL(z21, 1u);
	z51 -= 3u;
	BOOST_CHECK_EQUAL(z51, 0u);
	z51 -= z52;
	BOOST_CHECK_EQUAL(z51, 2u);

	//*
	BOOST_CHECK_EQUAL(z21 * z22, 0u);
	BOOST_CHECK_EQUAL(z21 * 3u, 1u);
	BOOST_CHECK_EQUAL(z22 * 3u, 0u);
	BOOST_CHECK_EQUAL(6u * z21, 0u);
	BOOST_CHECK_EQUAL(6u * z22, 0u);
	BOOST_CHECK_EQUAL(z51 * z52, 1u);
	BOOST_CHECK_EQUAL(z51 * 3u, 1u);
	BOOST_CHECK_EQUAL(z52 * 3u, 4u);
	BOOST_CHECK_EQUAL(7u * z51, 4u);
	BOOST_CHECK_EQUAL(7u * z52, 1u);
	z21 *= 3u;
	BOOST_CHECK_EQUAL(z21, 1u);
	z21 *= z22;
	BOOST_CHECK_EQUAL(z21, 0u);
	z51 *= 3u;
	BOOST_CHECK_EQUAL(z51, 1u);
	z51 *= z52;
	BOOST_CHECK_EQUAL(z51, 3u);

	//==
	BOOST_CHECK(z21 == z22);
	BOOST_CHECK(z21 == 0u);
	BOOST_CHECK(0u == z21);
	BOOST_CHECK(z22 == 0u);
	BOOST_CHECK(0u == z22);
	BOOST_CHECK(!(z21 == 1u));
	BOOST_CHECK(!(3u == z21));
	BOOST_CHECK(!(z22 == 1u));
	BOOST_CHECK(!(3u == z22));
	BOOST_CHECK(z51 == z52);
	BOOST_CHECK(z51 == 3u);
	BOOST_CHECK(3u == z51);
	BOOST_CHECK(z52 == 3u);
	BOOST_CHECK(3u == z52);
	BOOST_CHECK(!(z51 == 1u));
	BOOST_CHECK(!(7u == z51));
	BOOST_CHECK(!(z52 == 7u));
	BOOST_CHECK(!(1u == z52));
}

BOOST_AUTO_TEST_CASE(Field_other)
{
	Z2_field_element z21(7);
	Z2_field_element z22(2);
	Zp_field_element<5> z51(7);
	Zp_field_element<5> z52(3);
	Zp_field_element<7> z71(8);
	Zp_field_element<7> z72(3);

	BOOST_CHECK_EQUAL(z21.get_inverse(), 1u);
	BOOST_CHECK_EQUAL(z22.get_inverse(), 0u);
	BOOST_CHECK_EQUAL(z51.get_inverse(), 3u);
	BOOST_CHECK_EQUAL(z52.get_inverse(), 2u);
	BOOST_CHECK_EQUAL(z71.get_inverse(), 1u);
	BOOST_CHECK_EQUAL(z72.get_inverse(), 5u);

	BOOST_CHECK_EQUAL(z21.get_additive_identity(), 0u);
	BOOST_CHECK_EQUAL(z22.get_additive_identity(), 0u);
	BOOST_CHECK_EQUAL(z51.get_additive_identity(), 0u);
	BOOST_CHECK_EQUAL(z52.get_additive_identity(), 0u);
	BOOST_CHECK_EQUAL(z71.get_additive_identity(), 0u);
	BOOST_CHECK_EQUAL(z72.get_additive_identity(), 0u);
	BOOST_CHECK_EQUAL(z21.get_multiplicative_identity(), 1u);
	BOOST_CHECK_EQUAL(z22.get_multiplicative_identity(), 1u);
	BOOST_CHECK_EQUAL(z51.get_multiplicative_identity(), 1u);
	BOOST_CHECK_EQUAL(z52.get_multiplicative_identity(), 1u);
	BOOST_CHECK_EQUAL(z71.get_multiplicative_identity(), 1u);
	BOOST_CHECK_EQUAL(z72.get_multiplicative_identity(), 1u);

	BOOST_CHECK_EQUAL(z21.get_characteristic(), 2);
	BOOST_CHECK_EQUAL(z22.get_characteristic(), 2);
	BOOST_CHECK_EQUAL(z51.get_characteristic(), 5);
	BOOST_CHECK_EQUAL(z52.get_characteristic(), 5);
	BOOST_CHECK_EQUAL(z71.get_characteristic(), 7);
	BOOST_CHECK_EQUAL(z72.get_characteristic(), 7);

	BOOST_CHECK_EQUAL(z21.get_value(), 1);
	BOOST_CHECK_EQUAL(z22.get_value(), 0);
	BOOST_CHECK_EQUAL(z51.get_value(), 2);
	BOOST_CHECK_EQUAL(z52.get_value(), 3);
	BOOST_CHECK_EQUAL(z71.get_value(), 1);
	BOOST_CHECK_EQUAL(z72.get_value(), 3);
}
