/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <iostream>
#include <random>
#include <vector>
#include <utility>

#include "gudhi/matrix.h"
#include "gudhi/options.h"
#include "gudhi/utilities/Z2_field.h"
#include "gudhi/utilities/Zp_field.h"
#include "gudhi/utilities/utilities.h"

using Gudhi::persistence_matrix::Z2_field_element;
using Gudhi::persistence_matrix::Zp_field_element;
using Gudhi::persistence_matrix::Matrix;
using Gudhi::persistence_matrix::Representative_cycles_options;
using Gudhi::persistence_matrix::Default_options;
using Gudhi::persistence_matrix::Zigzag_options;
using Gudhi::persistence_matrix::Multi_persistence_options;
using Gudhi::persistence_matrix::Cohomology_persistence_options;
using Gudhi::persistence_matrix::Column_types;
using Gudhi::persistence_matrix::Bar;

using boundary_type = std::vector<unsigned int>;
template<class Field_type = Zp_field_element<5> >
using field_boundary_type = std::vector<std::pair<unsigned int,Field_type> >;

template<Column_types column_type = Column_types::SET, bool separated_by_dimension = false, bool parallelizable = false>
struct test_options1 : Default_options<Z2_field_element, column_type, separated_by_dimension, parallelizable>{
	static const bool has_column_pairings = true;
	static const bool has_vine_update = true;
};

template<Column_types column_type = Column_types::SET, bool separated_by_dimension = false, bool parallelizable = false>
struct test_options2 : Default_options<Z2_field_element, column_type, separated_by_dimension, parallelizable>{
	static const bool has_row_access = true;
	static const bool has_column_pairings = true;
	static const bool has_vine_update = true;
	static const bool is_of_boundary_type = false;
	static const bool has_removable_columns = true;
	static const bool is_indexed_by_position = true;
};

template<class Matrix_type>
void test_comp_zp(Matrix_type m)
{
	std::vector<field_boundary_type<> > ob;
	field_boundary_type<> fb;

	m.insert_boundary(fb);
	m.get_column(0);
	m.get_max_dimension();
	m.get_number_of_columns();
	m.get_column_dimension(0);
	m.add_to(0, 1);
	m.is_zero_cell(0, 0);
	m.is_zero_column(0);
	m.get_pivot(0);
	swap(m, m);
	m.print();
	m = Matrix_type(ob);
}

template<class Matrix_type>
void test_comp_z2(Matrix_type m)
{
	std::vector<boundary_type> ob;
	boundary_type fb;

	m.insert_boundary(fb);
	m.get_column(0);
	m.get_max_dimension();
	m.get_number_of_columns();
	m.get_column_dimension(0);
	m.add_to(0, 1);
	m.is_zero_cell(0, 0);
	m.is_zero_column(0);
	m.get_pivot(0);
	swap(m, m);
	m.print();
	m = Matrix_type(ob);
}

int main(int argc, char* const argv[]) {
	Zp_field_element<5> f(3);
	Zp_field_element<5> f2(7);

	std::clog << "== : " << (f == f2) << " " << (f == 3u) << " " << (f2 == 3u) << " " << (f == 7u) << "\n";

	std::clog << "+ : " << (f + f2) << " " << (f + 3u) << " " << (f2 + 3u) << " " << (7u + f) << "\n";
	std::clog << "- : " << (f - f2) << " " << (f - 3u) << " " << (f2 - 3u) << " " << (7u - f) << "\n";
	std::clog << "* : " << (f * f2) << " " << (f * 3u) << " " << (f2 * 3u) << " " << (7u * f) << "\n";

	f += f2;
	f2 += 3u;

	std::clog << "+= : " << f << " " << f2 << "\n";

	unsigned int a = 3;

	a = f;
	std::clog << "= : " << f << " " << a << "\n";

	std::vector<boundary_type> orderedBoundaries1;
	boundary_type b;
	orderedBoundaries1.emplace_back();
	orderedBoundaries1.emplace_back();
	orderedBoundaries1.emplace_back();
	orderedBoundaries1.push_back(boundary_type{0,1});
	orderedBoundaries1.push_back(boundary_type{1,2});

	std::vector<field_boundary_type<> > orderedBoundaries2;
	field_boundary_type<> fb;
	std::pair<unsigned int,Zp_field_element<5>> p;
	orderedBoundaries2.emplace_back();
	orderedBoundaries2.emplace_back();
	orderedBoundaries2.emplace_back();
	orderedBoundaries2.push_back(field_boundary_type<>{{0,3},{1,2}});
	orderedBoundaries2.push_back(field_boundary_type<>{{1,3},{2,2}});

	Matrix<Representative_cycles_options<Zp_field_element<5> > > m1(orderedBoundaries2);
	Matrix<Representative_cycles_options<Zp_field_element<2> > > m2(orderedBoundaries1);
	Matrix<Representative_cycles_options<Zp_field_element<5>,Column_types::LIST> > m3(orderedBoundaries2);
	Matrix<Representative_cycles_options<Zp_field_element<2>,Column_types::LIST> > m4(orderedBoundaries1);
	Matrix<Representative_cycles_options<Zp_field_element<5>,Column_types::UNORDERED_SET> > m5(orderedBoundaries2);
	Matrix<Representative_cycles_options<Zp_field_element<2>,Column_types::UNORDERED_SET> > m6(orderedBoundaries1);
	Matrix<Representative_cycles_options<Zp_field_element<5>,Column_types::VECTOR> > m7(orderedBoundaries2);
	Matrix<Representative_cycles_options<Zp_field_element<2>,Column_types::VECTOR> > m8(orderedBoundaries1);
	Matrix<Representative_cycles_options<Zp_field_element<2>,Column_types::HEAP> > m10(orderedBoundaries1);

	Matrix<Default_options<Zp_field_element<5> > > m11(orderedBoundaries2);
	Matrix<Default_options<Zp_field_element<2> > > m12(orderedBoundaries1);
	Matrix<Default_options<Zp_field_element<5>,Column_types::LIST> > m13(orderedBoundaries2);
	Matrix<Default_options<Zp_field_element<2>,Column_types::LIST> > m14(orderedBoundaries1);
	Matrix<Default_options<Zp_field_element<5>,Column_types::UNORDERED_SET> > m15(orderedBoundaries2);
	Matrix<Default_options<Zp_field_element<2>,Column_types::UNORDERED_SET> > m16(orderedBoundaries1);
	Matrix<Default_options<Zp_field_element<5>,Column_types::VECTOR> > m17(orderedBoundaries2);
	Matrix<Default_options<Zp_field_element<2>,Column_types::VECTOR> > m18(orderedBoundaries1);
	Matrix<Default_options<Zp_field_element<2>,Column_types::HEAP> > m20(orderedBoundaries1);

	Matrix<Multi_persistence_options<> > m21(orderedBoundaries1);
	Matrix<Multi_persistence_options<Column_types::LIST> > m22(orderedBoundaries1);
	Matrix<Multi_persistence_options<Column_types::UNORDERED_SET> > m23(orderedBoundaries1);
	Matrix<Multi_persistence_options<Column_types::VECTOR> > m24(orderedBoundaries1);
	Matrix<Multi_persistence_options<Column_types::HEAP> > m25(orderedBoundaries1);

	Matrix<Zigzag_options<> > m31(orderedBoundaries1);
	Matrix<Zigzag_options<Column_types::LIST> > m32(orderedBoundaries1);

	Matrix<Cohomology_persistence_options<Zp_field_element<5> > > m41(orderedBoundaries2);
	Matrix<Cohomology_persistence_options<Zp_field_element<2> > > m42(orderedBoundaries1);

	Matrix<test_options1<> > m51(orderedBoundaries1);
	Matrix<test_options1<Column_types::LIST> > m52(orderedBoundaries1);
	Matrix<test_options1<Column_types::UNORDERED_SET> > m53(orderedBoundaries1);
	Matrix<test_options1<Column_types::VECTOR> > m54(orderedBoundaries1);
	Matrix<test_options1<Column_types::HEAP> > m55(orderedBoundaries1);

	Matrix<test_options2<> > m61(orderedBoundaries1);
	Matrix<test_options2<Column_types::LIST> > m62(orderedBoundaries1);

	test_comp_zp(m1);
	test_comp_z2(m2);
	test_comp_zp(m3);
	test_comp_z2(m4);
	test_comp_zp(m5);
	test_comp_z2(m6);
	test_comp_zp(m7);
	test_comp_z2(m8);
	test_comp_z2(m10);

	auto& dgm1 = m1.get_current_barcode();
	m1.update_representative_cycles();
	m1.get_representative_cycles();
	m1.get_representative_cycle(dgm1.front());
	auto& dgm2 = m2.get_current_barcode();
	m2.update_representative_cycles();
	m2.get_representative_cycles();
	m2.get_representative_cycle(dgm2.front());
	auto& dgm3 = m3.get_current_barcode();
	m3.update_representative_cycles();
	m3.get_representative_cycles();
	m3.get_representative_cycle(dgm3.front());
	auto& dgm4 = m4.get_current_barcode();
	m4.update_representative_cycles();
	m4.get_representative_cycles();
	m4.get_representative_cycle(dgm4.front());
	auto& dgm5 = m5.get_current_barcode();
	m5.update_representative_cycles();
	m5.get_representative_cycles();
	m5.get_representative_cycle(dgm5.front());
	auto& dgm6 = m6.get_current_barcode();
	m6.update_representative_cycles();
	m6.get_representative_cycles();
	m6.get_representative_cycle(dgm6.front());
	auto& dgm7 = m7.get_current_barcode();
	m7.update_representative_cycles();
	m7.get_representative_cycles();
	m7.get_representative_cycle(dgm7.front());
	auto& dgm8 = m8.get_current_barcode();
	m8.update_representative_cycles();
	m8.get_representative_cycles();
	m8.get_representative_cycle(dgm8.front());
	auto& dgm10 = m10.get_current_barcode();
	m10.update_representative_cycles();
	m10.get_representative_cycles();
	m10.get_representative_cycle(dgm10.front());

	test_comp_zp(m11);
	test_comp_z2(m12);
	test_comp_zp(m13);
	test_comp_z2(m14);
	test_comp_zp(m15);
	test_comp_z2(m16);
	test_comp_zp(m17);
	test_comp_z2(m18);
	test_comp_z2(m20);

	m11.zero_cell(0, 0);
	m11.zero_column(0);
	m12.zero_cell(0, 0);
	m12.zero_column(0);
	m13.zero_cell(0, 0);
	m13.zero_column(0);
	m14.zero_cell(0, 0);
	m14.zero_column(0);
	m15.zero_cell(0, 0);
	m15.zero_column(0);
	m16.zero_cell(0, 0);
	m16.zero_column(0);
	m17.zero_cell(0, 0);
	m17.zero_column(0);
	m18.zero_cell(0, 0);
	m18.zero_column(0);
	m20.zero_cell(0, 0);
	m20.zero_column(0);

//	std::cout << "number of columns2: " << m21.get_number_of_columns() << "\n";
	test_comp_z2(m21);
//	std::cout << "number of columns3: " << m21.get_number_of_columns() << "\n";
	test_comp_z2(m22);
	test_comp_z2(m23);
	test_comp_z2(m24);
	test_comp_z2(m25);

//	std::cout << "number of columns4: " << m21.get_number_of_columns() << "\n";
	m21.get_column_with_pivot(0);
//	std::cout << "number of columns5: " << m21.get_number_of_columns() << "\n";
	m22.get_column_with_pivot(0);
	m23.get_column_with_pivot(0);
	m24.get_column_with_pivot(0);
	m25.get_column_with_pivot(0);
	m21.get_current_barcode();
	m21.vine_swap_with_z_eq_1_case(3);
	m21.vine_swap(3);
	m22.get_current_barcode();
	m22.vine_swap_with_z_eq_1_case(3);
	m22.vine_swap(3);
	m23.get_current_barcode();
	m23.vine_swap_with_z_eq_1_case(3);
	m23.vine_swap(3);
	m24.get_current_barcode();
	m24.vine_swap_with_z_eq_1_case(3);
	m24.vine_swap(3);
	m25.get_current_barcode();
	m25.vine_swap_with_z_eq_1_case(3);
	m25.vine_swap(3);

	test_comp_z2(m31);
	test_comp_z2(m32);

	m31.get_column_with_pivot(0);
	m32.get_column_with_pivot(0);
	m31.get_current_barcode();
	m31.vine_swap_with_z_eq_1_case(3,4);
	m31.vine_swap(3,4);
	m32.get_current_barcode();
	m32.vine_swap_with_z_eq_1_case(3,4);
	m32.vine_swap(3,4);
	m31.get_row(0);
	m31.erase_last();
	m32.get_row(0);
	m32.erase_last();

	test_comp_zp(m41);
	test_comp_z2(m42);

//	std::cout << "number of columns2: " << m51.get_number_of_columns() << "\n";
	test_comp_z2(m51);
//	std::cout << "number of columns3: " << m51.get_number_of_columns() << "\n";
	test_comp_z2(m52);
	test_comp_z2(m53);
	test_comp_z2(m54);
	test_comp_z2(m55);
//std::cout << "number of columns2: " << m51.get_number_of_columns() << "\n";
	m51.get_column_with_pivot(0);
	m52.get_column_with_pivot(0);
	m53.get_column_with_pivot(0);
	m54.get_column_with_pivot(0);
	m55.get_column_with_pivot(0);
//	std::cout << "number of columns3: " << m51.get_number_of_columns() << "\n";
	m51.get_current_barcode();
//	std::cout << "number of columns4: " << m51.get_number_of_columns() << "\n";
	m51.vine_swap_with_z_eq_1_case(3,4);
	m51.vine_swap(3,4);
	m52.get_current_barcode();
	m52.vine_swap_with_z_eq_1_case(3,4);
	m52.vine_swap(3,4);
	m53.get_current_barcode();
	m53.vine_swap_with_z_eq_1_case(3,4);
	m53.vine_swap(3,4);
	m54.get_current_barcode();
	m54.vine_swap_with_z_eq_1_case(3,4);
	m54.vine_swap(3,4);
	m55.get_current_barcode();
	m55.vine_swap_with_z_eq_1_case(3,4);
	m55.vine_swap(3,4);

	test_comp_z2(m61);
	test_comp_z2(m62);

	m61.get_column_with_pivot(0);
	m62.get_column_with_pivot(0);
	m61.get_current_barcode();
	m61.vine_swap_with_z_eq_1_case(3);
	m61.vine_swap(3);
	m62.get_current_barcode();
	m62.vine_swap_with_z_eq_1_case(3);
	m62.vine_swap(3);
	m61.get_row(0);
	m61.erase_last();
	m62.get_row(0);
	m62.erase_last();

	return 0;
}
