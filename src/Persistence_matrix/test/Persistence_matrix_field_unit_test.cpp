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

#include "gudhi/Fields/Z2_field.h"
#include "gudhi/Fields/Zp_field.h"
#include "gudhi/Fields/Zp_field_shared.h"
#include "gudhi/Fields/Multi_field.h"
#include "gudhi/Fields/Multi_field_small.h"
#include "gudhi/Fields/Multi_field_shared.h"
#include "gudhi/Fields/Multi_field_small_shared.h"

using Gudhi::persistence_matrix::Z2_field_element;
using Gudhi::persistence_matrix::Zp_field_element;
using Gudhi::persistence_matrix::Shared_Zp_field_element;
using Gudhi::persistence_matrix::Multi_field_element;
using Gudhi::persistence_matrix::Multi_field_element_with_small_characteristics;
using Gudhi::persistence_matrix::Shared_multi_field_element;
using Gudhi::persistence_matrix::Shared_multi_field_element_with_small_characteristics;

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
	Zp_field_element<13> z13_m1(5);
	Zp_field_element<13> z13_m2(std::move(z13_m1));
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

BOOST_AUTO_TEST_CASE(Shared_Field_constructors)
{
	Shared_Zp_field_element::initialize(2);

	//default constructor
	Shared_Zp_field_element e;
	BOOST_CHECK_EQUAL(e, 0u);

	//value constructor
	Shared_Zp_field_element e2(5);
	BOOST_CHECK_EQUAL(e2, 1u);
	e = 5;
	BOOST_CHECK_EQUAL(e, 1u);

	//copy constructor
	Shared_Zp_field_element e3(e);
	Shared_Zp_field_element e4 = e;
	BOOST_CHECK_EQUAL(e3, 1u);
	BOOST_CHECK_EQUAL(e4, 1u);

	//move constructor
	Shared_Zp_field_element e5(std::move(e2));
	BOOST_CHECK_EQUAL(e5, 1u);
	BOOST_CHECK_EQUAL(e2, 0u);

	//swap
	swap(e2, e3);
	BOOST_CHECK_EQUAL(e2, 1u);
	BOOST_CHECK_EQUAL(e3, 0u);

	Shared_Zp_field_element::initialize(5);

	//default constructor
	Shared_Zp_field_element f;
	BOOST_CHECK_EQUAL(f, 0u);

	//value constructor
	Shared_Zp_field_element f2(7);
	BOOST_CHECK_EQUAL(f2, 2u);
	f = 7;
	BOOST_CHECK_EQUAL(f, 2u);

	//copy constructor
	Shared_Zp_field_element f3(f);
	Shared_Zp_field_element f4 = f;
	BOOST_CHECK_EQUAL(f3, 2u);
	BOOST_CHECK_EQUAL(f4, 2u);

	//move constructor
	Shared_Zp_field_element f5(std::move(f2));
	BOOST_CHECK_EQUAL(f5, 2u);
	BOOST_CHECK_EQUAL(f2, 0u);

	//swap
	swap(f2, f3);
	BOOST_CHECK_EQUAL(f2, 2u);
	BOOST_CHECK_EQUAL(f3, 0u);

	Shared_Zp_field_element::initialize(13);

	//default constructor
	Shared_Zp_field_element g;
	BOOST_CHECK_EQUAL(g, 0u);

	//value constructor
	Shared_Zp_field_element g2(5);
	BOOST_CHECK_EQUAL(g2, 5u);
	g = 5;
	BOOST_CHECK_EQUAL(g, 5u);

	//copy constructor
	Shared_Zp_field_element g3(g);
	Shared_Zp_field_element g4 = g;
	BOOST_CHECK_EQUAL(g3, 5u);
	BOOST_CHECK_EQUAL(g4, 5u);

	//move constructor
	Shared_Zp_field_element g5(std::move(g2));
	BOOST_CHECK_EQUAL(g5, 5u);
	BOOST_CHECK_EQUAL(g2, 0u);

	//swap
	swap(g2, g3);
	BOOST_CHECK_EQUAL(g2, 5u);
	BOOST_CHECK_EQUAL(g3, 0u);
}

BOOST_AUTO_TEST_CASE(Shared_Field_operators)
{
	Shared_Zp_field_element::initialize(2);

	Shared_Zp_field_element z21(7);
	Shared_Zp_field_element z22(2);

	//+
	BOOST_CHECK_EQUAL(z21 + z22, 1u);
	BOOST_CHECK_EQUAL(z21 + 3u, 0u);
	BOOST_CHECK_EQUAL(z22 + 3u, 1u);
	BOOST_CHECK_EQUAL(6u + z21, 1u);
	BOOST_CHECK_EQUAL(6u + z22, 0u);
	z21 += 3u;
	BOOST_CHECK_EQUAL(z21, 0u);
	z21 += z22;
	BOOST_CHECK_EQUAL(z21, 0u);

	//-
	BOOST_CHECK_EQUAL(z21 - z22, 0u);
	BOOST_CHECK_EQUAL(z21 - 3u, 1u);
	BOOST_CHECK_EQUAL(z22 - 3u, 1u);
	BOOST_CHECK_EQUAL(6u - z21, 0u);
	BOOST_CHECK_EQUAL(6u - z22, 0u);
	z21 -= 3u;
	BOOST_CHECK_EQUAL(z21, 1u);
	z21 -= z22;
	BOOST_CHECK_EQUAL(z21, 1u);

	//*
	BOOST_CHECK_EQUAL(z21 * z22, 0u);
	BOOST_CHECK_EQUAL(z21 * 3u, 1u);
	BOOST_CHECK_EQUAL(z22 * 3u, 0u);
	BOOST_CHECK_EQUAL(6u * z21, 0u);
	BOOST_CHECK_EQUAL(6u * z22, 0u);
	z21 *= 3u;
	BOOST_CHECK_EQUAL(z21, 1u);
	z21 *= z22;
	BOOST_CHECK_EQUAL(z21, 0u);

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

	Shared_Zp_field_element::initialize(5);

	Shared_Zp_field_element z51(7);
	Shared_Zp_field_element z52(3);

	//+
	BOOST_CHECK_EQUAL(z51 + z52, 0u);
	BOOST_CHECK_EQUAL(z51 + 3u, 0u);
	BOOST_CHECK_EQUAL(z52 + 3u, 1u);
	BOOST_CHECK_EQUAL(7u + z51, 4u);
	BOOST_CHECK_EQUAL(7u + z52, 0u);
	z51 += 3u;
	BOOST_CHECK_EQUAL(z51, 0u);
	z51 += z52;
	BOOST_CHECK_EQUAL(z51, 3u);

	//-
	BOOST_CHECK_EQUAL(z51 - z52, 0u);
	BOOST_CHECK_EQUAL(z51 - 3u, 0u);
	BOOST_CHECK_EQUAL(z52 - 3u, 0u);
	BOOST_CHECK_EQUAL(7u - z51, 4u);
	BOOST_CHECK_EQUAL(7u - z52, 4u);
	z51 -= 3u;
	BOOST_CHECK_EQUAL(z51, 0u);
	z51 -= z52;
	BOOST_CHECK_EQUAL(z51, 2u);

	//*
	BOOST_CHECK_EQUAL(z51 * z52, 1u);
	BOOST_CHECK_EQUAL(z51 * 3u, 1u);
	BOOST_CHECK_EQUAL(z52 * 3u, 4u);
	BOOST_CHECK_EQUAL(7u * z51, 4u);
	BOOST_CHECK_EQUAL(7u * z52, 1u);
	z51 *= 3u;
	BOOST_CHECK_EQUAL(z51, 1u);
	z51 *= z52;
	BOOST_CHECK_EQUAL(z51, 3u);

	//==
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

BOOST_AUTO_TEST_CASE(Shared_Field_other)
{
	Shared_Zp_field_element::initialize(2);

	Shared_Zp_field_element z21(7);
	Shared_Zp_field_element z22(2);

	BOOST_CHECK_EQUAL(z21.get_inverse(), 1u);
	BOOST_CHECK_EQUAL(z22.get_inverse(), 0u);

	BOOST_CHECK_EQUAL(z21.get_additive_identity(), 0u);
	BOOST_CHECK_EQUAL(z22.get_additive_identity(), 0u);
	BOOST_CHECK_EQUAL(z21.get_multiplicative_identity(), 1u);
	BOOST_CHECK_EQUAL(z22.get_multiplicative_identity(), 1u);

	BOOST_CHECK_EQUAL(z21.get_characteristic(), 2);
	BOOST_CHECK_EQUAL(z22.get_characteristic(), 2);

	BOOST_CHECK_EQUAL(z21.get_value(), 1);
	BOOST_CHECK_EQUAL(z22.get_value(), 0);

	Shared_Zp_field_element::initialize(5);

	Shared_Zp_field_element z51(7);
	Shared_Zp_field_element z52(3);

	BOOST_CHECK_EQUAL(z51.get_inverse(), 3u);
	BOOST_CHECK_EQUAL(z52.get_inverse(), 2u);

	BOOST_CHECK_EQUAL(z51.get_additive_identity(), 0u);
	BOOST_CHECK_EQUAL(z52.get_additive_identity(), 0u);
	BOOST_CHECK_EQUAL(z51.get_multiplicative_identity(), 1u);
	BOOST_CHECK_EQUAL(z52.get_multiplicative_identity(), 1u);

	BOOST_CHECK_EQUAL(z51.get_characteristic(), 5);
	BOOST_CHECK_EQUAL(z52.get_characteristic(), 5);

	BOOST_CHECK_EQUAL(z51.get_value(), 2);
	BOOST_CHECK_EQUAL(z52.get_value(), 3);

	Shared_Zp_field_element::initialize(7);

	Shared_Zp_field_element z71(8);
	Shared_Zp_field_element z72(3);

	BOOST_CHECK_EQUAL(z71.get_inverse(), 1u);
	BOOST_CHECK_EQUAL(z72.get_inverse(), 5u);

	BOOST_CHECK_EQUAL(z71.get_additive_identity(), 0u);
	BOOST_CHECK_EQUAL(z72.get_additive_identity(), 0u);
	BOOST_CHECK_EQUAL(z71.get_multiplicative_identity(), 1u);
	BOOST_CHECK_EQUAL(z72.get_multiplicative_identity(), 1u);

	BOOST_CHECK_EQUAL(z71.get_characteristic(), 7);
	BOOST_CHECK_EQUAL(z72.get_characteristic(), 7);

	BOOST_CHECK_EQUAL(z71.get_value(), 1);
	BOOST_CHECK_EQUAL(z72.get_value(), 3);
}

BOOST_AUTO_TEST_CASE(Multi_Field_constructors)
{
	//default constructor
	Multi_field_element<5,13> m_d;
	BOOST_CHECK_EQUAL(m_d, mpz_class(0));

	//value constructor
	Multi_field_element<5,13> m_v(5006);
	BOOST_CHECK_EQUAL(m_v, mpz_class(1));

	//copy constructor
	Multi_field_element<5,13> m_c1(5006);
	Multi_field_element<5,13> m_c2 = m_c1;
	BOOST_CHECK_EQUAL(m_c2, mpz_class(1));
	Multi_field_element<5,13> m_c3(m_c2);
	BOOST_CHECK_EQUAL(m_c3, mpz_class(1));

	//move constructor
	Multi_field_element<5,13> m_m1(5006);
	Multi_field_element<5,13> m_m2(std::move(m_m1));
	BOOST_CHECK_EQUAL(m_m2, mpz_class(1));
	BOOST_CHECK_EQUAL(m_m1, mpz_class(0));

	//swap
	Multi_field_element<5,13> m_s1(5006);
	Multi_field_element<5,13> m_s2(5005);
	swap(m_s1, m_s2);
	BOOST_CHECK_EQUAL(m_s2, mpz_class(1));
	BOOST_CHECK_EQUAL(m_s1, mpz_class(0));

	//default constructor
	Multi_field_element_with_small_characteristics<5,13> ms_d;
	BOOST_CHECK_EQUAL(ms_d, 0u);

	//value constructor
	Multi_field_element_with_small_characteristics<5,13> ms_v(5006);
	BOOST_CHECK_EQUAL(ms_v, 1u);

	//copy constructor
	Multi_field_element_with_small_characteristics<5,13> ms_c1(5006);
	Multi_field_element_with_small_characteristics<5,13> ms_c2 = ms_c1;
	BOOST_CHECK_EQUAL(ms_c2, 1u);
	Multi_field_element_with_small_characteristics<5,13> ms_c3(m_c2);
	BOOST_CHECK_EQUAL(ms_c3, 1u);

	//move constructor
	Multi_field_element_with_small_characteristics<5,13> ms_m1(5006);
	Multi_field_element_with_small_characteristics<5,13> ms_m2(std::move(ms_m1));
	BOOST_CHECK_EQUAL(ms_m2, 1u);
	BOOST_CHECK_EQUAL(ms_m1, 0u);

	//swap
	Multi_field_element_with_small_characteristics<5,13> ms_s1(5006);
	Multi_field_element_with_small_characteristics<5,13> ms_s2(5005);
	swap(ms_s1, ms_s2);
	BOOST_CHECK_EQUAL(ms_s2, 1u);
	BOOST_CHECK_EQUAL(ms_s1, 0u);
}

BOOST_AUTO_TEST_CASE(Multi_Field_operators)
{
	Multi_field_element<5,13> m1(5005);
	Multi_field_element<5,13> m2(5007);

	//+
	BOOST_CHECK_EQUAL(m1 + m2, mpz_class(2));
	BOOST_CHECK_EQUAL(m1 + mpz_class(3), mpz_class(3));
	BOOST_CHECK_EQUAL(m2 + mpz_class(3), mpz_class(5));
	BOOST_CHECK_EQUAL(mpz_class(6) + m1, mpz_class(6));
	BOOST_CHECK_EQUAL(mpz_class(6) + m2, mpz_class(8));
	m1 += mpz_class(3);
	BOOST_CHECK_EQUAL(m1, mpz_class(3));
	m1 += m2;
	BOOST_CHECK_EQUAL(m1, mpz_class(5));

	//-
	BOOST_CHECK_EQUAL(m1 - m2, mpz_class(3));
	BOOST_CHECK_EQUAL(m1 - mpz_class(3), mpz_class(2));
	BOOST_CHECK_EQUAL(m2 - mpz_class(3), mpz_class(5004));
	BOOST_CHECK_EQUAL(mpz_class(6) - m1, mpz_class(1));
	BOOST_CHECK_EQUAL(mpz_class(6) - m2, mpz_class(4));
	m2 -= mpz_class(3);
	BOOST_CHECK_EQUAL(m2, mpz_class(5004));
	m2 -= m1;
	BOOST_CHECK_EQUAL(m2, mpz_class(4999));

	//*
	BOOST_CHECK_EQUAL(m1 * m2, mpz_class(4975));
	BOOST_CHECK_EQUAL(m1 * mpz_class(3), mpz_class(15));
	BOOST_CHECK_EQUAL(m2 * mpz_class(3), mpz_class(4987));
	BOOST_CHECK_EQUAL(mpz_class(6) * m1, mpz_class(30));
	BOOST_CHECK_EQUAL(mpz_class(6) * m2, mpz_class(4969));
	m1 *= mpz_class(3);
	BOOST_CHECK_EQUAL(m1, mpz_class(15));
	m1 *= m2;
	BOOST_CHECK_EQUAL(m1, mpz_class(4915));

	//==
	BOOST_CHECK(!(m1 == m2));
	m2 -= mpz_class(84);
	BOOST_CHECK(m1 == m2);
	BOOST_CHECK(m1 == mpz_class(4915));
	BOOST_CHECK(mpz_class(4915) == m1);
	BOOST_CHECK(m2 == mpz_class(4915));
	BOOST_CHECK(mpz_class(4915) == m2);
	BOOST_CHECK(!(m1 == mpz_class(1)));
	BOOST_CHECK(!(mpz_class(3) == m1));
	BOOST_CHECK(!(m2 == mpz_class(1)));
	BOOST_CHECK(!(mpz_class(3) == m2));

	Multi_field_element_with_small_characteristics<5,13> ms1(5005);
	Multi_field_element_with_small_characteristics<5,13> ms2(5007);

	//+
	BOOST_CHECK_EQUAL(ms1 + ms2, 2u);
	BOOST_CHECK_EQUAL(ms1 + 3u, 3u);
	BOOST_CHECK_EQUAL(ms2 + 3u, 5u);
	BOOST_CHECK_EQUAL(6u + ms1, 6u);
	BOOST_CHECK_EQUAL(6u + ms2, 8u);
	ms1 += 3u;
	BOOST_CHECK_EQUAL(ms1, 3u);
	ms1 += ms2;
	BOOST_CHECK_EQUAL(ms1, 5u);

	//-
	BOOST_CHECK_EQUAL(ms1 - ms2, 3u);
	BOOST_CHECK_EQUAL(ms1 - 3u, 2u);
	BOOST_CHECK_EQUAL(ms2 - 3u, 5004u);
	BOOST_CHECK_EQUAL(6u - ms1, 1u);
	BOOST_CHECK_EQUAL(6u - ms2, 4u);
	ms2 -= 3u;
	BOOST_CHECK_EQUAL(ms2, 5004u);
	ms2 -= ms1;
	BOOST_CHECK_EQUAL(ms2, 4999u);

	//*
	BOOST_CHECK_EQUAL(ms1 * ms2, 4975u);
	BOOST_CHECK_EQUAL(ms1 * 3u, 15u);
	BOOST_CHECK_EQUAL(ms2 * 3u, 4987u);
	BOOST_CHECK_EQUAL(6u * ms1, 30u);
	BOOST_CHECK_EQUAL(6u * ms2, 4969u);
	ms1 *= 3u;
	BOOST_CHECK_EQUAL(ms1, 15u);
	ms1 *= ms2;
	BOOST_CHECK_EQUAL(ms1, 4915u);

	//==
	BOOST_CHECK(!(ms1 == ms2));
	ms2 -= 84u;
	BOOST_CHECK(ms1 == ms2);
	BOOST_CHECK(ms1 == 4915u);
	BOOST_CHECK(4915u == ms1);
	BOOST_CHECK(ms2 == 4915u);
	BOOST_CHECK(4915u == ms2);
	BOOST_CHECK(!(ms1 == 1u));
	BOOST_CHECK(!(3u == ms1));
	BOOST_CHECK(!(ms2 == 1u));
	BOOST_CHECK(!(3u == ms2));
}

BOOST_AUTO_TEST_CASE(Multi_Field_other)
{
	Multi_field_element<5,13> m1(1);
	Multi_field_element<5,13> m2(7);

	BOOST_CHECK_EQUAL(m1.get_inverse(), mpz_class(1));
	BOOST_CHECK_EQUAL(m2.get_inverse(), mpz_class(2758));
//	auto inv = m1.get_partial_inverse(35);
//	std::cout << inv.first << ", " << inv.second << "\n";
//	inv = m2.get_partial_inverse(35);
//	std::cout << inv.first << ", " << inv.second << "\n";
	BOOST_CHECK(m1.get_partial_inverse(35) == std::make_pair(Multi_field_element<5,13>(1716), mpz_class(35)));
	BOOST_CHECK(m2.get_partial_inverse(35) == std::make_pair(Multi_field_element<5,13>(3003), mpz_class(5)));

	BOOST_CHECK_EQUAL(m1.get_additive_identity(), mpz_class(0));
	BOOST_CHECK_EQUAL(m2.get_additive_identity(), mpz_class(0));
	BOOST_CHECK_EQUAL(m1.get_multiplicative_identity(), mpz_class(1));
	BOOST_CHECK_EQUAL(m2.get_multiplicative_identity(), mpz_class(1));
	BOOST_CHECK_EQUAL(m1.get_partial_multiplicative_identity(), mpz_class(0));
	BOOST_CHECK_EQUAL(m2.get_partial_multiplicative_identity(), mpz_class(715));

	BOOST_CHECK_EQUAL(m1.get_characteristic(), 5005);
	BOOST_CHECK_EQUAL(m2.get_characteristic(), 5005);

	BOOST_CHECK_EQUAL(m1.get_value(), 1);
	BOOST_CHECK_EQUAL(m2.get_value(), 7);

	Multi_field_element_with_small_characteristics<5,13> ms1(1);
	Multi_field_element_with_small_characteristics<5,13> ms2(7);

	BOOST_CHECK_EQUAL(ms1.get_inverse(), 1u);
	BOOST_CHECK_EQUAL(ms2.get_inverse(), 2758u);
	BOOST_CHECK(ms1.get_partial_inverse(35) == std::make_pair(Multi_field_element_with_small_characteristics<5,13>(1716), 35u));
	BOOST_CHECK(ms2.get_partial_inverse(35) == std::make_pair(Multi_field_element_with_small_characteristics<5,13>(3003), 5u));

	BOOST_CHECK_EQUAL(ms1.get_additive_identity(), 0u);
	BOOST_CHECK_EQUAL(ms2.get_additive_identity(), 0u);
	BOOST_CHECK_EQUAL(ms1.get_multiplicative_identity(), 1u);
	BOOST_CHECK_EQUAL(ms2.get_multiplicative_identity(), 1u);
	BOOST_CHECK_EQUAL(ms1.get_partial_multiplicative_identity(), 0u);
	BOOST_CHECK_EQUAL(ms2.get_partial_multiplicative_identity(), 715u);

	BOOST_CHECK_EQUAL(ms1.get_characteristic(), 5005u);
	BOOST_CHECK_EQUAL(ms2.get_characteristic(), 5005u);

	BOOST_CHECK_EQUAL(ms1.get_value(), 1u);
	BOOST_CHECK_EQUAL(ms2.get_value(), 7u);

	Multi_field_element<3,30> mb1(2);
	Multi_field_element_with_small_characteristics<3,30> mb2(2);

	BOOST_CHECK_EQUAL(mb1.get_characteristic(), mb2.get_characteristic());	// == 3234846615
	BOOST_CHECK_EQUAL(mb1.get_partial_inverse(35).first.get_value(), mb2.get_partial_inverse(35).first.get_value());	// == 2033332158
}

BOOST_AUTO_TEST_CASE(Shared_Multi_Field_constructors)
{
	Shared_multi_field_element::initialize(5, 13);

	//default constructor
	Shared_multi_field_element m_d;
	BOOST_CHECK_EQUAL(m_d, mpz_class(0));

	//value constructor
	Shared_multi_field_element m_v(5006);
	BOOST_CHECK_EQUAL(m_v, mpz_class(1));

	//copy constructor
	Shared_multi_field_element m_c1(5006);
	Shared_multi_field_element m_c2 = m_c1;
	BOOST_CHECK_EQUAL(m_c2, mpz_class(1));
	Shared_multi_field_element m_c3(m_c2);
	BOOST_CHECK_EQUAL(m_c3, mpz_class(1));

	//move constructor
	Shared_multi_field_element m_m1(5006);
	Shared_multi_field_element m_m2(std::move(m_m1));
	BOOST_CHECK_EQUAL(m_m2, mpz_class(1));
	BOOST_CHECK_EQUAL(m_m1, mpz_class(0));

	//swap
	Shared_multi_field_element m_s1(5006);
	Shared_multi_field_element m_s2(5005);
	swap(m_s1, m_s2);
	BOOST_CHECK_EQUAL(m_s2, mpz_class(1));
	BOOST_CHECK_EQUAL(m_s1, mpz_class(0));

	Shared_multi_field_element_with_small_characteristics::initialize(5, 13);

	//default constructor
	Shared_multi_field_element_with_small_characteristics ms_d;
	BOOST_CHECK_EQUAL(ms_d, 0u);

	//value constructor
	Shared_multi_field_element_with_small_characteristics ms_v(5006);
	BOOST_CHECK_EQUAL(ms_v, 1u);

	//copy constructor
	Shared_multi_field_element_with_small_characteristics ms_c1(5006);
	Shared_multi_field_element_with_small_characteristics ms_c2 = ms_c1;
	BOOST_CHECK_EQUAL(ms_c2, 1u);
	Shared_multi_field_element_with_small_characteristics ms_c3(m_c2);
	BOOST_CHECK_EQUAL(ms_c3, 1u);

	//move constructor
	Shared_multi_field_element_with_small_characteristics ms_m1(5006);
	Shared_multi_field_element_with_small_characteristics ms_m2(std::move(ms_m1));
	BOOST_CHECK_EQUAL(ms_m2, 1u);
	BOOST_CHECK_EQUAL(ms_m1, 0u);

	//swap
	Shared_multi_field_element_with_small_characteristics ms_s1(5006);
	Shared_multi_field_element_with_small_characteristics ms_s2(5005);
	swap(ms_s1, ms_s2);
	BOOST_CHECK_EQUAL(ms_s2, 1u);
	BOOST_CHECK_EQUAL(ms_s1, 0u);
}

BOOST_AUTO_TEST_CASE(Shared_Multi_Field_operators)
{
	Shared_multi_field_element::initialize(5, 13);

	Shared_multi_field_element m1(5005);
	Shared_multi_field_element m2(5007);

	//+
	BOOST_CHECK_EQUAL(m1 + m2, mpz_class(2));
	BOOST_CHECK_EQUAL(m1 + mpz_class(3), mpz_class(3));
	BOOST_CHECK_EQUAL(m2 + mpz_class(3), mpz_class(5));
	BOOST_CHECK_EQUAL(mpz_class(6) + m1, mpz_class(6));
	BOOST_CHECK_EQUAL(mpz_class(6) + m2, mpz_class(8));
	m1 += mpz_class(3);
	BOOST_CHECK_EQUAL(m1, mpz_class(3));
	m1 += m2;
	BOOST_CHECK_EQUAL(m1, mpz_class(5));

	//-
	BOOST_CHECK_EQUAL(m1 - m2, mpz_class(3));
	BOOST_CHECK_EQUAL(m1 - mpz_class(3), mpz_class(2));
	BOOST_CHECK_EQUAL(m2 - mpz_class(3), mpz_class(5004));
	BOOST_CHECK_EQUAL(mpz_class(6) - m1, mpz_class(1));
	BOOST_CHECK_EQUAL(mpz_class(6) - m2, mpz_class(4));
	m2 -= mpz_class(3);
	BOOST_CHECK_EQUAL(m2, mpz_class(5004));
	m2 -= m1;
	BOOST_CHECK_EQUAL(m2, mpz_class(4999));

	//*
	BOOST_CHECK_EQUAL(m1 * m2, mpz_class(4975));
	BOOST_CHECK_EQUAL(m1 * mpz_class(3), mpz_class(15));
	BOOST_CHECK_EQUAL(m2 * mpz_class(3), mpz_class(4987));
	BOOST_CHECK_EQUAL(mpz_class(6) * m1, mpz_class(30));
	BOOST_CHECK_EQUAL(mpz_class(6) * m2, mpz_class(4969));
	m1 *= mpz_class(3);
	BOOST_CHECK_EQUAL(m1, mpz_class(15));
	m1 *= m2;
	BOOST_CHECK_EQUAL(m1, mpz_class(4915));

	//==
	BOOST_CHECK(!(m1 == m2));
	m2 -= mpz_class(84);
	BOOST_CHECK(m1 == m2);
	BOOST_CHECK(m1 == mpz_class(4915));
	BOOST_CHECK(mpz_class(4915) == m1);
	BOOST_CHECK(m2 == mpz_class(4915));
	BOOST_CHECK(mpz_class(4915) == m2);
	BOOST_CHECK(!(m1 == mpz_class(1)));
	BOOST_CHECK(!(mpz_class(3) == m1));
	BOOST_CHECK(!(m2 == mpz_class(1)));
	BOOST_CHECK(!(mpz_class(3) == m2));

	Shared_multi_field_element_with_small_characteristics::initialize(5, 13);

	Shared_multi_field_element_with_small_characteristics ms1(5005);
	Shared_multi_field_element_with_small_characteristics ms2(5007);

	//+
	BOOST_CHECK_EQUAL(ms1 + ms2, 2u);
	BOOST_CHECK_EQUAL(ms1 + 3u, 3u);
	BOOST_CHECK_EQUAL(ms2 + 3u, 5u);
	BOOST_CHECK_EQUAL(6u + ms1, 6u);
	BOOST_CHECK_EQUAL(6u + ms2, 8u);
	ms1 += 3u;
	BOOST_CHECK_EQUAL(ms1, 3u);
	ms1 += ms2;
	BOOST_CHECK_EQUAL(ms1, 5u);

	//-
	BOOST_CHECK_EQUAL(ms1 - ms2, 3u);
	BOOST_CHECK_EQUAL(ms1 - 3u, 2u);
	BOOST_CHECK_EQUAL(ms2 - 3u, 5004u);
	BOOST_CHECK_EQUAL(6u - ms1, 1u);
	BOOST_CHECK_EQUAL(6u - ms2, 4u);
	ms2 -= 3u;
	BOOST_CHECK_EQUAL(ms2, 5004u);
	ms2 -= ms1;
	BOOST_CHECK_EQUAL(ms2, 4999u);

	//*
	BOOST_CHECK_EQUAL(ms1 * ms2, 4975u);
	BOOST_CHECK_EQUAL(ms1 * 3u, 15u);
	BOOST_CHECK_EQUAL(ms2 * 3u, 4987u);
	BOOST_CHECK_EQUAL(6u * ms1, 30u);
	BOOST_CHECK_EQUAL(6u * ms2, 4969u);
	ms1 *= 3u;
	BOOST_CHECK_EQUAL(ms1, 15u);
	ms1 *= ms2;
	BOOST_CHECK_EQUAL(ms1, 4915u);

	//==
	BOOST_CHECK(!(ms1 == ms2));
	ms2 -= 84u;
	BOOST_CHECK(ms1 == ms2);
	BOOST_CHECK(ms1 == 4915u);
	BOOST_CHECK(4915u == ms1);
	BOOST_CHECK(ms2 == 4915u);
	BOOST_CHECK(4915u == ms2);
	BOOST_CHECK(!(ms1 == 1u));
	BOOST_CHECK(!(3u == ms1));
	BOOST_CHECK(!(ms2 == 1u));
	BOOST_CHECK(!(3u == ms2));
}

BOOST_AUTO_TEST_CASE(Shared_Multi_Field_other)
{
	Shared_multi_field_element::initialize(5, 13);

	Shared_multi_field_element m1(1);
	Shared_multi_field_element m2(7);

	BOOST_CHECK_EQUAL(m1.get_inverse(), mpz_class(1));
	BOOST_CHECK_EQUAL(m2.get_inverse(), mpz_class(2758));
//	auto inv = m1.get_partial_inverse(35);
//	std::cout << inv.first << ", " << inv.second << "\n";
//	inv = m2.get_partial_inverse(35);
//	std::cout << inv.first << ", " << inv.second << "\n";
	BOOST_CHECK(m1.get_partial_inverse(35) == std::make_pair(Shared_multi_field_element(1716), mpz_class(35)));
	BOOST_CHECK(m2.get_partial_inverse(35) == std::make_pair(Shared_multi_field_element(3003), mpz_class(5)));

	BOOST_CHECK_EQUAL(m1.get_additive_identity(), mpz_class(0));
	BOOST_CHECK_EQUAL(m2.get_additive_identity(), mpz_class(0));
	BOOST_CHECK_EQUAL(m1.get_multiplicative_identity(), mpz_class(1));
	BOOST_CHECK_EQUAL(m2.get_multiplicative_identity(), mpz_class(1));
	BOOST_CHECK_EQUAL(m1.get_partial_multiplicative_identity(), mpz_class(0));
	BOOST_CHECK_EQUAL(m2.get_partial_multiplicative_identity(), mpz_class(715));

	BOOST_CHECK_EQUAL(m1.get_characteristic(), 5005);
	BOOST_CHECK_EQUAL(m2.get_characteristic(), 5005);

	BOOST_CHECK_EQUAL(m1.get_value(), 1);
	BOOST_CHECK_EQUAL(m2.get_value(), 7);

	Shared_multi_field_element_with_small_characteristics::initialize(5, 13);

	Shared_multi_field_element_with_small_characteristics ms1(1);
	Shared_multi_field_element_with_small_characteristics ms2(7);

	BOOST_CHECK_EQUAL(ms1.get_inverse(), 1u);
	BOOST_CHECK_EQUAL(ms2.get_inverse(), 2758u);
	BOOST_CHECK(ms1.get_partial_inverse(35) == std::make_pair(Shared_multi_field_element_with_small_characteristics(1716), 35u));
	BOOST_CHECK(ms2.get_partial_inverse(35) == std::make_pair(Shared_multi_field_element_with_small_characteristics(3003), 5u));

	BOOST_CHECK_EQUAL(ms1.get_additive_identity(), 0u);
	BOOST_CHECK_EQUAL(ms2.get_additive_identity(), 0u);
	BOOST_CHECK_EQUAL(ms1.get_multiplicative_identity(), 1u);
	BOOST_CHECK_EQUAL(ms2.get_multiplicative_identity(), 1u);
	BOOST_CHECK_EQUAL(ms1.get_partial_multiplicative_identity(), 0u);
	BOOST_CHECK_EQUAL(ms2.get_partial_multiplicative_identity(), 715u);

	BOOST_CHECK_EQUAL(ms1.get_characteristic(), 5005u);
	BOOST_CHECK_EQUAL(ms2.get_characteristic(), 5005u);

	BOOST_CHECK_EQUAL(ms1.get_value(), 1u);
	BOOST_CHECK_EQUAL(ms2.get_value(), 7u);

	Shared_multi_field_element mb1(2);
	Shared_multi_field_element_with_small_characteristics mb2(2);

	BOOST_CHECK_EQUAL(mb1.get_characteristic(), mb2.get_characteristic());	// == 3234846615
	BOOST_CHECK_EQUAL(mb1.get_partial_inverse(35).first.get_value(), mb2.get_partial_inverse(35).first.get_value());	// == 2033332158
}
