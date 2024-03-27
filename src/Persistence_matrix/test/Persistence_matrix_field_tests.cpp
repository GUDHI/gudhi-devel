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

template<class Z2>
void test_z2_standart_field_constructors(){
	//default constructor
	Z2 z2_d;
	BOOST_CHECK_EQUAL(z2_d, 0);

	//value constructor
	Z2 z2_v(5);
	BOOST_CHECK_EQUAL(z2_v, 1);
	z2_d = 5;
	BOOST_CHECK_EQUAL(z2_d, 1);

	//copy constructor
	Z2 z2_c1(z2_d);
	BOOST_CHECK_EQUAL(z2_c1, 1);
	Z2 z2_c2 = z2_c1;
	BOOST_CHECK_EQUAL(z2_c2, 1);

	//move constructor
	Z2 z2_m1(5);
	Z2 z2_m2(std::move(z2_m1));
	BOOST_CHECK_EQUAL(z2_m2, 1);
	BOOST_CHECK_EQUAL(z2_m1, 0);

	//swap
	Z2 z2_s1(5);
	Z2 z2_s2(8);
	swap(z2_s1, z2_s2);
	BOOST_CHECK_EQUAL(z2_s2, 1);
	BOOST_CHECK_EQUAL(z2_s1, 0);
}

template<class Z5>
void test_z5_standart_field_constructors(){
	//default constructor
	Z5 z5_d;
	BOOST_CHECK_EQUAL(z5_d, 0);

	//value constructor
	Z5 z5_v(7);
	BOOST_CHECK_EQUAL(z5_v, 2);
	z5_d = 7;
	BOOST_CHECK_EQUAL(z5_d, 2);

	//copy constructor
	Z5 z5_c1(7);
	Z5 z5_c2 = z5_c1;
	BOOST_CHECK_EQUAL(z5_c2, 2);
	Z5 z5_c3(z5_c2);
	BOOST_CHECK_EQUAL(z5_c3, 2);

	//move constructor
	Z5 z5_m1(7);
	Z5 z5_m2(std::move(z5_m1));
	BOOST_CHECK_EQUAL(z5_m2, 2);
	BOOST_CHECK_EQUAL(z5_m1, 0);

	//swap
	Z5 z5_s1(4);
	Z5 z5_s2(8);
	swap(z5_s1, z5_s2);
	BOOST_CHECK_EQUAL(z5_s2, 4);
	BOOST_CHECK_EQUAL(z5_s1, 3);
}

template<class Z13>
void test_z13_standart_field_constructors(){
	//default constructor
	Z13 z13_d;
	BOOST_CHECK_EQUAL(z13_d, 0);

	//value constructor
	Z13 z13_v(5);
	BOOST_CHECK_EQUAL(z13_v, 5);
	z13_d = 5;
	BOOST_CHECK_EQUAL(z13_d, 5);

	//copy constructor
	Z13 z13_c1(5);
	Z13 z13_c2 = z13_c1;
	BOOST_CHECK_EQUAL(z13_c2, 5);
	Z13 z13_c3(z13_c2);
	BOOST_CHECK_EQUAL(z13_c3, 5);

	//move constructor
	Z13 z13_m1(5);
	Z13 z13_m2(std::move(z13_m1));
	BOOST_CHECK_EQUAL(z13_m2, 5);
	BOOST_CHECK_EQUAL(z13_m1, 0);

	//swap
	Z13 z13_s1(4);
	Z13 z13_s2(22);
	swap(z13_s1, z13_s2);
	BOOST_CHECK_EQUAL(z13_s2, 4);
	BOOST_CHECK_EQUAL(z13_s1, 9);
}

template<class Z2>
void test_z2_standart_field_operators(){
	Z2 z21(7);
	Z2 z22(2);

	//+
	BOOST_CHECK_EQUAL(z21 + z22, 1);
	BOOST_CHECK_EQUAL(z21 + 3, 0);
	BOOST_CHECK_EQUAL(z22 + 3, 1);
	BOOST_CHECK_EQUAL(6 + z21, 1);
	BOOST_CHECK_EQUAL(6 + z22, 0);
	z21 += 3;
	BOOST_CHECK_EQUAL(z21, 0);
	z21 += z22;
	BOOST_CHECK_EQUAL(z21, 0);

	//-
	BOOST_CHECK_EQUAL(z21 - z22, 0);
	BOOST_CHECK_EQUAL(z21 - 3, 1);
	BOOST_CHECK_EQUAL(z22 - 3, 1);
	BOOST_CHECK_EQUAL(6 - z21, 0);
	BOOST_CHECK_EQUAL(6 - z22, 0);
	z21 -= 3;
	BOOST_CHECK_EQUAL(z21, 1);
	z21 -= z22;
	BOOST_CHECK_EQUAL(z21, 1);

	//*
	BOOST_CHECK_EQUAL(z21 * z22, 0);
	BOOST_CHECK_EQUAL(z21 * 3, 1);
	BOOST_CHECK_EQUAL(z22 * 3, 0);
	BOOST_CHECK_EQUAL(6 * z21, 0);
	BOOST_CHECK_EQUAL(6 * z22, 0);
	z21 *= 3;
	BOOST_CHECK_EQUAL(z21, 1);
	z21 *= z22;
	BOOST_CHECK_EQUAL(z21, 0);

	//==
	BOOST_CHECK(z21 == z22);
	BOOST_CHECK(z21 == 0);
	BOOST_CHECK(0 == z21);
	BOOST_CHECK(z22 == 0);
	BOOST_CHECK(0 == z22);
	BOOST_CHECK(z21 != 1);
	BOOST_CHECK(3 != z21);
	BOOST_CHECK(z22 != 1);
	BOOST_CHECK(3 != z22);
}

template<class Z5>
void test_z5_standart_field_operators(){
	Z5 z51(7);
	Z5 z52(3);

	//+
	BOOST_CHECK_EQUAL(z51 + z52, 0);
	BOOST_CHECK_EQUAL(z51 + 3, 0);
	BOOST_CHECK_EQUAL(z52 + 3, 1);
	BOOST_CHECK_EQUAL(7 + z51, 4);
	BOOST_CHECK_EQUAL(7 + z52, 0);
	z51 += 3;
	BOOST_CHECK_EQUAL(z51, 0);
	z51 += z52;
	BOOST_CHECK_EQUAL(z51, 3);

	//-
	BOOST_CHECK_EQUAL(z51 - z52, 0);
	BOOST_CHECK_EQUAL(z51 - 3, 0);
	BOOST_CHECK_EQUAL(z52 - 3, 0);
	BOOST_CHECK_EQUAL(7 - z51, 4);
	BOOST_CHECK_EQUAL(7 - z52, 4);
	z51 -= 3;
	BOOST_CHECK_EQUAL(z51, 0);
	z51 -= z52;
	BOOST_CHECK_EQUAL(z51, 2);

	//*
	BOOST_CHECK_EQUAL(z51 * z52, 1);
	BOOST_CHECK_EQUAL(z51 * 3, 1);
	BOOST_CHECK_EQUAL(z52 * 3, 4);
	BOOST_CHECK_EQUAL(7 * z51, 4);
	BOOST_CHECK_EQUAL(7 * z52, 1);
	z51 *= 3;
	BOOST_CHECK_EQUAL(z51, 1);
	z51 *= z52;
	BOOST_CHECK_EQUAL(z51, 3);

	//==
	BOOST_CHECK(z51 == z52);
	BOOST_CHECK(z51 == 3);
	BOOST_CHECK(3 == z51);
	BOOST_CHECK(z52 == 3);
	BOOST_CHECK(3 == z52);
	BOOST_CHECK(z51 != 1);
	BOOST_CHECK(7 != z51);
	BOOST_CHECK(z52 != 7);
	BOOST_CHECK(1 != z52);
}

template<class Z2>
void test_z2_standart_field_properties(){
	Z2 z21(7);
	Z2 z22(2);

	BOOST_CHECK_EQUAL(z21.get_inverse(), 1);
	BOOST_CHECK_EQUAL(z22.get_inverse(), 0);
	BOOST_CHECK(z21.get_partial_inverse(35) == std::make_pair(Z2(1), 35u));
	BOOST_CHECK(z22.get_partial_inverse(35) == std::make_pair(Z2(0), 35u));

	BOOST_CHECK_EQUAL(z21.get_additive_identity(), 0);
	BOOST_CHECK_EQUAL(z22.get_additive_identity(), 0);

	BOOST_CHECK_EQUAL(z21.get_multiplicative_identity(), 1);
	BOOST_CHECK_EQUAL(z22.get_multiplicative_identity(), 1);
	BOOST_CHECK_EQUAL(z21.get_partial_multiplicative_identity(7), 1);
	BOOST_CHECK_EQUAL(z22.get_partial_multiplicative_identity(7), 1);

	BOOST_CHECK_EQUAL(z21.get_characteristic(), 2);
	BOOST_CHECK_EQUAL(z22.get_characteristic(), 2);

	BOOST_CHECK_EQUAL(z21.get_value(), 1);
	BOOST_CHECK_EQUAL(z22.get_value(), 0);
}

template<class Z5>
void test_z5_standart_field_properties(){
	Z5 z51(7);
	Z5 z52(3);

	BOOST_CHECK_EQUAL(z51.get_inverse(), 3);
	BOOST_CHECK_EQUAL(z52.get_inverse(), 2);
	BOOST_CHECK(z51.get_partial_inverse(35) == std::make_pair(Z5(3), 35u));
	BOOST_CHECK(z52.get_partial_inverse(35) == std::make_pair(Z5(2), 35u));

	BOOST_CHECK_EQUAL(z51.get_additive_identity(), 0);
	BOOST_CHECK_EQUAL(z52.get_additive_identity(), 0);

	BOOST_CHECK_EQUAL(z51.get_multiplicative_identity(), 1);
	BOOST_CHECK_EQUAL(z52.get_multiplicative_identity(), 1);
	BOOST_CHECK_EQUAL(z51.get_partial_multiplicative_identity(7), 1);
	BOOST_CHECK_EQUAL(z52.get_partial_multiplicative_identity(7), 1);

	BOOST_CHECK_EQUAL(z51.get_characteristic(), 5);
	BOOST_CHECK_EQUAL(z52.get_characteristic(), 5);

	BOOST_CHECK_EQUAL(z51.get_value(), 2);
	BOOST_CHECK_EQUAL(z52.get_value(), 3);
}

template<class Z7>
void test_z7_standart_field_properties(){
	Z7 z71(8);
	Z7 z72(3);

	BOOST_CHECK_EQUAL(z71.get_inverse(), 1);
	BOOST_CHECK_EQUAL(z72.get_inverse(), 5);
	BOOST_CHECK(z71.get_partial_inverse(35) == std::make_pair(Z7(1), 35u));
	BOOST_CHECK(z72.get_partial_inverse(35) == std::make_pair(Z7(5), 35u));

	BOOST_CHECK_EQUAL(z71.get_additive_identity(), 0);
	BOOST_CHECK_EQUAL(z72.get_additive_identity(), 0);

	BOOST_CHECK_EQUAL(z71.get_multiplicative_identity(), 1);
	BOOST_CHECK_EQUAL(z72.get_multiplicative_identity(), 1);
	BOOST_CHECK_EQUAL(z71.get_partial_multiplicative_identity(7), 1);
	BOOST_CHECK_EQUAL(z72.get_partial_multiplicative_identity(7), 1);

	BOOST_CHECK_EQUAL(z71.get_characteristic(), 7);
	BOOST_CHECK_EQUAL(z72.get_characteristic(), 7);

	BOOST_CHECK_EQUAL(z71.get_value(), 1);
	BOOST_CHECK_EQUAL(z72.get_value(), 3);
}

BOOST_AUTO_TEST_CASE(Field_constructors)
{
	test_z2_standart_field_constructors<Z2_field_element>();
	test_z2_standart_field_constructors<Zp_field_element<2> >();
	test_z5_standart_field_constructors<Zp_field_element<5> >();
	test_z13_standart_field_constructors<Zp_field_element<13> >();
}

BOOST_AUTO_TEST_CASE(Field_operators)
{
	test_z2_standart_field_operators<Z2_field_element>();
	test_z2_standart_field_operators<Zp_field_element<2> >();
	test_z5_standart_field_operators<Zp_field_element<5> >();
}

BOOST_AUTO_TEST_CASE(Field_properties)
{
	test_z2_standart_field_properties<Z2_field_element>();
	test_z2_standart_field_properties<Zp_field_element<2> >();
	test_z5_standart_field_properties<Zp_field_element<5> >();
	test_z7_standart_field_properties<Zp_field_element<7> >();
}

BOOST_AUTO_TEST_CASE(Shared_Field_constructors)
{
	Shared_Zp_field_element::initialize(2);
	test_z2_standart_field_constructors<Shared_Zp_field_element>();

	Shared_Zp_field_element::initialize(5);
	test_z5_standart_field_constructors<Shared_Zp_field_element>();

	Shared_Zp_field_element::initialize(13);
	test_z13_standart_field_constructors<Shared_Zp_field_element>();
}

BOOST_AUTO_TEST_CASE(Shared_Field_operators)
{
	Shared_Zp_field_element::initialize(2);
	test_z2_standart_field_operators<Shared_Zp_field_element>();

	Shared_Zp_field_element::initialize(5);
	test_z5_standart_field_operators<Shared_Zp_field_element>();
}

BOOST_AUTO_TEST_CASE(Shared_Field_properties)
{
	Shared_Zp_field_element::initialize(2);
	test_z2_standart_field_properties<Shared_Zp_field_element>();

	Shared_Zp_field_element::initialize(5);
	test_z5_standart_field_properties<Shared_Zp_field_element>();

	Shared_Zp_field_element::initialize(7);
	test_z7_standart_field_properties<Shared_Zp_field_element>();
}

template<class MF>
void test_multi_field_constructors(){
	using T = typename MF::element_type;

	//default constructor
	Multi_field_element<5,13> m_d;
	BOOST_CHECK_EQUAL(m_d, T(0));

	//value constructor
	Multi_field_element<5,13> m_v(5006);
	BOOST_CHECK_EQUAL(m_v, T(1));

	//copy constructor
	Multi_field_element<5,13> m_c1(5006);
	Multi_field_element<5,13> m_c2 = m_c1;
	BOOST_CHECK_EQUAL(m_c2, T(1));
	Multi_field_element<5,13> m_c3(m_c2);
	BOOST_CHECK_EQUAL(m_c3, T(1));

	//move constructor
	Multi_field_element<5,13> m_m1(5006);
	Multi_field_element<5,13> m_m2(std::move(m_m1));
	BOOST_CHECK_EQUAL(m_m2, T(1));
	BOOST_CHECK_EQUAL(m_m1, T(0));

	//swap
	Multi_field_element<5,13> m_s1(5006);
	Multi_field_element<5,13> m_s2(5005);
	swap(m_s1, m_s2);
	BOOST_CHECK_EQUAL(m_s2, T(1));
	BOOST_CHECK_EQUAL(m_s1, T(0));
}

template<class MF>
void test_multi_field_operators(){
	using T = typename MF::element_type;

	MF m1(5005);
	MF m2(5007);

	//+
	BOOST_CHECK_EQUAL(m1 + m2, T(2));
	BOOST_CHECK_EQUAL(m1 + T(3), T(3));
	BOOST_CHECK_EQUAL(m2 + T(3), T(5));
	BOOST_CHECK_EQUAL(T(6) + m1, T(6));
	BOOST_CHECK_EQUAL(T(6) + m2, T(8));
	m1 += T(3);
	BOOST_CHECK_EQUAL(m1, T(3));
	m1 += m2;
	BOOST_CHECK_EQUAL(m1, T(5));

	//-
	BOOST_CHECK_EQUAL(m1 - m2, T(3));
	BOOST_CHECK_EQUAL(m1 - T(3), T(2));
	BOOST_CHECK_EQUAL(m2 - T(3), T(5004));
	BOOST_CHECK_EQUAL(T(6) - m1, T(1));
	BOOST_CHECK_EQUAL(T(6) - m2, T(4));
	m2 -= T(3);
	BOOST_CHECK_EQUAL(m2, T(5004));
	m2 -= m1;
	BOOST_CHECK_EQUAL(m2, T(4999));

	//*
	BOOST_CHECK_EQUAL(m1 * m2, T(4975));
	BOOST_CHECK_EQUAL(m1 * T(3), T(15));
	BOOST_CHECK_EQUAL(m2 * T(3), T(4987));
	BOOST_CHECK_EQUAL(T(6) * m1, T(30));
	BOOST_CHECK_EQUAL(T(6) * m2, T(4969));
	m1 *= T(3);
	BOOST_CHECK_EQUAL(m1, T(15));
	m1 *= m2;
	BOOST_CHECK_EQUAL(m1, T(4915));

	//==
	BOOST_CHECK(m1 != m2);
	m2 -= T(84);
	BOOST_CHECK(m1 == m2);
	BOOST_CHECK(m1 == T(4915));
	BOOST_CHECK(T(4915) == m1);
	BOOST_CHECK(m2 == T(4915));
	BOOST_CHECK(T(4915) == m2);
	BOOST_CHECK(m1 != T(1));
	BOOST_CHECK(T(3) != m1);
	BOOST_CHECK(m2 != T(1));
	BOOST_CHECK(T(3) != m2);
}

template<class MF>
void test_multi_field_properties(){
	using T = typename MF::element_type;

	MF m1(1);
	MF m2(7);

	BOOST_CHECK_EQUAL(m1.get_inverse(), T(1));
	BOOST_CHECK_EQUAL(m2.get_inverse(), T(2758));
	BOOST_CHECK(m1.get_partial_inverse(35) == std::make_pair(MF(1716), T(35)));
	BOOST_CHECK(m2.get_partial_inverse(35) == std::make_pair(MF(3003), T(5)));

	BOOST_CHECK_EQUAL(m1.get_additive_identity(), T(0));
	BOOST_CHECK_EQUAL(m2.get_additive_identity(), T(0));
	BOOST_CHECK_EQUAL(m1.get_multiplicative_identity(), T(1));
	BOOST_CHECK_EQUAL(m2.get_multiplicative_identity(), T(1));
	BOOST_CHECK_EQUAL(m1.get_partial_multiplicative_identity(7), T(715));
	BOOST_CHECK_EQUAL(m2.get_partial_multiplicative_identity(7), T(715));

	BOOST_CHECK_EQUAL(m1.get_characteristic(), 5005);
	BOOST_CHECK_EQUAL(m2.get_characteristic(), 5005);

	BOOST_CHECK_EQUAL(m1.get_value(), 1);
	BOOST_CHECK_EQUAL(m2.get_value(), 7);
}

BOOST_AUTO_TEST_CASE(Multi_Field_constructors)
{
	test_multi_field_constructors<Multi_field_element<5,13> >();
	test_multi_field_constructors<Multi_field_element_with_small_characteristics<5,13> >();
}

BOOST_AUTO_TEST_CASE(Multi_Field_operators)
{
	test_multi_field_operators<Multi_field_element<5,13> >();
	test_multi_field_operators<Multi_field_element_with_small_characteristics<5,13> >();
}

BOOST_AUTO_TEST_CASE(Multi_Field_properties)
{
	test_multi_field_properties<Multi_field_element<5,13> >();
	test_multi_field_properties<Multi_field_element_with_small_characteristics<5,13> >();

	Multi_field_element<3,30> mb1(2);
	Multi_field_element_with_small_characteristics<3,30> mb2(2);

	BOOST_CHECK_EQUAL(mb1.get_characteristic(), mb2.get_characteristic());	// == 3234846615
	BOOST_CHECK_EQUAL(mb1.get_partial_inverse(35).first.get_value(), mb2.get_partial_inverse(35).first.get_value());	// == 2033332158
}

BOOST_AUTO_TEST_CASE(Shared_Multi_Field_constructors)
{
	Shared_multi_field_element::initialize(5, 13);
	test_multi_field_constructors<Shared_multi_field_element>();

	Shared_multi_field_element_with_small_characteristics::initialize(5, 13);
	test_multi_field_constructors<Shared_multi_field_element>();
}

BOOST_AUTO_TEST_CASE(Shared_Multi_Field_operators)
{
	Shared_multi_field_element::initialize(5, 13);
	test_multi_field_operators<Shared_multi_field_element>();

	Shared_multi_field_element_with_small_characteristics::initialize(5, 13);
	test_multi_field_operators<Shared_multi_field_element>();
}

BOOST_AUTO_TEST_CASE(Shared_Multi_Field_properties)
{
	Shared_multi_field_element::initialize(5, 13);
	test_multi_field_properties<Shared_multi_field_element>();

	Shared_multi_field_element_with_small_characteristics::initialize(5, 13);
	test_multi_field_properties<Shared_multi_field_element>();

	Shared_multi_field_element::initialize(3, 30);
	Shared_multi_field_element_with_small_characteristics::initialize(3, 30);
	Shared_multi_field_element mb1(2);
	Shared_multi_field_element_with_small_characteristics mb2(2);

	BOOST_CHECK_EQUAL(mb1.get_characteristic(), mb2.get_characteristic());	// == 3234846615
	BOOST_CHECK_EQUAL(mb1.get_partial_inverse(35).first.get_value(), mb2.get_partial_inverse(35).first.get_value());	// == 2033332158
}
