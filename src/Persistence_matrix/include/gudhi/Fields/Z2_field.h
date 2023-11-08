/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef MATRIX_FIELD_Z2_H_
#define MATRIX_FIELD_Z2_H_

#include <utility>

namespace Gudhi {
namespace persistence_matrix {

class Z2_field_element {
public:
	using element_type = unsigned int;
	template <class T>
	using isInteger = std::enable_if_t<std::is_integral_v<T> >;

	Z2_field_element();
	Z2_field_element(unsigned int element);
	Z2_field_element(int element);
	Z2_field_element(const Z2_field_element& toCopy);
	Z2_field_element(Z2_field_element&& toMove) noexcept;

	friend void operator+=(Z2_field_element& f1, Z2_field_element const& f2){
		f1.element_ = (f1.element_ != f2.element_);
	}
	friend Z2_field_element operator+(Z2_field_element f1, Z2_field_element const& f2){
		f1 += f2;
		return f1;
	}
	friend void operator+=(Z2_field_element& f, unsigned int const v){
		f.element_ = (f.element_ != (v % 2));
	}
	//v will be casted into an unsigned int
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend Z2_field_element operator+(Z2_field_element f, const Integer_type v){
		f += v;
		return f;
	}
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend Integer_type operator+(const Integer_type v, Z2_field_element const& f){
		return f.element_ != (v % 2);
	}

	friend void operator-=(Z2_field_element& f1, Z2_field_element const& f2){
		f1.element_ = (f1.element_ != f2.element_);
	}
	friend Z2_field_element operator-(Z2_field_element f1, Z2_field_element const& f2){
		f1 -= f2;
		return f1;
	}
	friend void operator-=(Z2_field_element& f, unsigned int const v){
		f.element_ = (f.element_ != (v % 2));
	}
	//v will be casted into an unsigned int
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend Z2_field_element operator-(Z2_field_element f, const Integer_type v){
		f -= v;
		return f;
	}
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend Integer_type operator-(const Integer_type v, Z2_field_element const& f){
		return f.element_ != (v % 2);
	}

	friend void operator*=(Z2_field_element& f1, Z2_field_element const& f2){
		f1.element_ = (f1.element_ && f2.element_);
	}
	friend Z2_field_element operator*(Z2_field_element f1, Z2_field_element const& f2){
		f1 *= f2;
		return f1;
	}
	friend void operator*=(Z2_field_element& f, unsigned int const v){
		f.element_ = (f.element_ && (v % 2));
	}
	//v will be casted into an unsigned int
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend Z2_field_element operator*(Z2_field_element f, const Integer_type v){
		f *= v;
		return f;
	}
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend Integer_type operator*(const Integer_type v, Z2_field_element const& f){
		return f.element_ && (v % 2);
	}

	friend bool operator==(const Z2_field_element& f1, const Z2_field_element& f2){
		return f1.element_ == f2.element_;
	}
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend bool operator==(const Integer_type v, const Z2_field_element& f){
		return (v % 2) == f.element_;
	}
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend bool operator==(const Z2_field_element& f, const Integer_type v){
		return (v % 2) == f.element_;
	}
	friend bool operator!=(const Z2_field_element& f1, const Z2_field_element& f2){
		return !(f1 == f2);
	}
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend bool operator!=(const Integer_type v, const Z2_field_element& f){
		return !(v == f);
	}
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend bool operator!=(const Z2_field_element& f, const Integer_type v){
		return !(v == f);
	}

	Z2_field_element& operator=(Z2_field_element other);
	Z2_field_element& operator=(const unsigned int value);
	operator unsigned int() const;
	friend void swap(Z2_field_element& f1, Z2_field_element& f2){
		std::swap(f1.element_, f2.element_);
	}

	Z2_field_element get_inverse() const;
	std::pair<Z2_field_element, unsigned int> get_partial_inverse(unsigned int product_of_characteristics) const;

	static Z2_field_element get_additive_identity();
	static Z2_field_element get_multiplicative_identity();
	static Z2_field_element get_partial_multiplicative_identity();
	static constexpr unsigned int get_characteristic();

	unsigned int get_value() const;

	static constexpr bool handles_only_z2(){
		return true;
	}

private:
	bool element_;
};

inline Z2_field_element::Z2_field_element()
	: element_(false)
{}

inline Z2_field_element::Z2_field_element(unsigned int element)
	: element_(element & 1U)
{}

inline Z2_field_element::Z2_field_element(int element)
	: element_(element & 1U)
{}

inline Z2_field_element::Z2_field_element(const Z2_field_element &toCopy)
	: element_(toCopy.element_)
{}

inline Z2_field_element::Z2_field_element(Z2_field_element &&toMove) noexcept
	: element_(std::exchange(toMove.element_, 0))
{}

//inline Z2_field_element &Z2_field_element::operator+=(Z2_field_element const &f)
//{
//	element_ = (element_ != f.element_);
//	return *this;
//}

//inline Z2_field_element &Z2_field_element::operator+=(unsigned int const v)
//{
//	element_ = (element_ != (v % 2));
//	return *this;
//}

//inline Z2_field_element &Z2_field_element::operator-=(const Z2_field_element &f)
//{
//	element_ = (element_ != f.element_);
//	return *this;
//}

//inline Z2_field_element &Z2_field_element::operator-=(unsigned int const v)
//{
//	element_ = (element_ != (v % 2));
//	return *this;
//}

//inline Z2_field_element &Z2_field_element::operator*=(const Z2_field_element &f)
//{
//	element_ = (element_ && f.element_);
//	return *this;
//}

//inline Z2_field_element &Z2_field_element::operator*=(unsigned int const v)
//{
//	element_ = (element_ && (v % 2));
//	return *this;
//}

inline Z2_field_element &Z2_field_element::operator=(Z2_field_element other)
{
	std::swap(element_, other.element_);
	return *this;
}

inline Z2_field_element &Z2_field_element::operator=(unsigned int const value)
{
	element_ = value & 1U;
	return *this;
}

inline Z2_field_element::operator unsigned int() const
{
	return element_;
}

inline Z2_field_element Z2_field_element::get_inverse() const
{
	return element_ ? Z2_field_element(1) : Z2_field_element();
}

inline std::pair<Z2_field_element, unsigned int>
Z2_field_element::get_partial_inverse(unsigned int product_of_characteristics) const
{
	return {get_inverse(), product_of_characteristics};
}

inline Z2_field_element Z2_field_element::get_additive_identity()
{
	return Z2_field_element();
}

inline Z2_field_element Z2_field_element::get_multiplicative_identity()
{
	return Z2_field_element(1);
}

inline Z2_field_element Z2_field_element::get_partial_multiplicative_identity()
{
	return Z2_field_element(1);
}

inline constexpr unsigned int Z2_field_element::get_characteristic()
{
	return 2;
}

inline unsigned int Z2_field_element::get_value() const
{
	return element_;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif  // MATRIX_FIELD_Z2_H_
