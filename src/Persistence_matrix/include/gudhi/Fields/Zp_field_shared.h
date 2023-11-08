/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef MATRIX_FIELD_ZP_VAR_H_
#define MATRIX_FIELD_ZP_VAR_H_

#include <utility>
#include <vector>
#include <limits.h>
#include <iostream>

namespace Gudhi {
namespace persistence_matrix {

class Shared_Zp_field_element {
public:
	using element_type = unsigned int;
	template <class T>
	using isInteger = std::enable_if_t<std::is_integral_v<T> >;

	Shared_Zp_field_element();
	Shared_Zp_field_element(unsigned int element);
	Shared_Zp_field_element(int element);	//only works if characteristic can be contained in an int
	Shared_Zp_field_element(const Shared_Zp_field_element& toCopy);
	Shared_Zp_field_element(Shared_Zp_field_element&& toMove) noexcept;

	static void initialize(unsigned int characteristic);

	friend void operator+=(Shared_Zp_field_element& f1, Shared_Zp_field_element const &f2){
		f1.element_ = Shared_Zp_field_element::_add(f1.element_, f2.element_);
	}
	friend Shared_Zp_field_element operator+(Shared_Zp_field_element f1, Shared_Zp_field_element const& f2){
		f1 += f2;
		return f1;
	}
	friend void operator+=(Shared_Zp_field_element& f, unsigned int const v){
		f.element_ = Shared_Zp_field_element::_add(f.element_, v < characteristic_ ? v : v % characteristic_);
	}
	//v is assumed to be positive and will be casted into an unsigned int
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend Shared_Zp_field_element operator+(Shared_Zp_field_element f, const Integer_type v){
		f += v;
		return f;
	}
	//v is assumed to be positive
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend Integer_type operator+(Integer_type v, Shared_Zp_field_element const& f){
		v += f.element_;
		if (v >= characteristic_) v %= characteristic_;
		return v;
	}

	friend void operator-=(Shared_Zp_field_element& f1, Shared_Zp_field_element const &f2){
		f1.element_ = Shared_Zp_field_element::_substract(f1.element_, f2.element_);
	}
	friend Shared_Zp_field_element operator-(Shared_Zp_field_element f1, Shared_Zp_field_element const& f2){
		f1 -= f2;
		return f1;
	}
	friend void operator-=(Shared_Zp_field_element& f, unsigned int const v){
		f.element_ = Shared_Zp_field_element::_substract(f.element_, v < characteristic_ ? v : v % characteristic_);
	}
	//v is assumed to be positive and will be casted into an unsigned int
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend Shared_Zp_field_element operator-(Shared_Zp_field_element f, const Integer_type v){
		f -= v;
		return f;
	}
	//v is assumed to be positive
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend Integer_type operator-(Integer_type v, Shared_Zp_field_element const& f){
		if (v >= characteristic_) v %= characteristic_;
		if (f.element_ > v) v += characteristic_;
		v -= f.element_;
		return v;
	}

	friend void operator*=(Shared_Zp_field_element& f1, Shared_Zp_field_element const &f2){
		f1.element_ = Shared_Zp_field_element::_multiply(f1.element_, f2.element_);
	}
	friend Shared_Zp_field_element operator*(Shared_Zp_field_element f1, Shared_Zp_field_element const& f2){
		f1 *= f2;
		return f1;
	}
	friend void operator*=(Shared_Zp_field_element& f, unsigned int const v){
		f.element_ = Shared_Zp_field_element::_multiply(f.element_, v < characteristic_ ? v : v % characteristic_);
	}
	//v is assumed to be positive and will be casted into an unsigned int
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend Shared_Zp_field_element operator*(Shared_Zp_field_element f, const Integer_type v){
		f *= v;
		return f;
	}
	//uses bitwise operations on v, so be carefull with signed integers
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend Integer_type operator*(Integer_type v, Shared_Zp_field_element const& f){
		unsigned int b = f.element_;
		unsigned int res = 0;
		unsigned int temp_b;

		while (v != 0) {
			if (v & 1) {
				if (b >= characteristic_ - res)
					res -= characteristic_;
				res += b;
			}
			v >>= 1;

			temp_b = b;
			if (b >= characteristic_ - b)
				temp_b -= characteristic_;
			b += temp_b;
		}

		return res;
	}

	friend bool operator==(const Shared_Zp_field_element& f1, const Shared_Zp_field_element& f2){
		return f1.element_ == f2.element_;
	}
	//v is assumed to be positive
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend bool operator==(const Integer_type v, const Shared_Zp_field_element& f){
		if (v < characteristic_) return v == f.element_;
		return (v % characteristic_) == f.element_;
	}
	//v is assumed to be positive
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend bool operator==(const Shared_Zp_field_element& f, const Integer_type v){
		if (v < characteristic_) return v == f.element_;
		return (v % characteristic_) == f.element_;
	}
	friend bool operator!=(const Shared_Zp_field_element& f1, const Shared_Zp_field_element& f2){
		return !(f1 == f2);
	}
	//v is assumed to be positive
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend bool operator!=(const Integer_type v, const Shared_Zp_field_element& f){
		return !(v == f);
	}
	//v is assumed to be positive
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend bool operator!=(const Shared_Zp_field_element& f, const Integer_type v){
		return !(v == f);
	}

	Shared_Zp_field_element& operator=(Shared_Zp_field_element other);
	Shared_Zp_field_element& operator=(const unsigned int value);
	operator unsigned int() const;
	friend void swap(Shared_Zp_field_element& f1, Shared_Zp_field_element& f2){
		std::swap(f1.element_, f2.element_);
	}

	Shared_Zp_field_element get_inverse() const;
	std::pair<Shared_Zp_field_element, unsigned int> get_partial_inverse(unsigned int product_of_characteristics) const;

	static Shared_Zp_field_element get_additive_identity();
	static Shared_Zp_field_element get_multiplicative_identity();
	static Shared_Zp_field_element get_partial_multiplicative_identity();
	static unsigned int get_characteristic();

	unsigned int get_value() const;

	static constexpr bool handles_only_z2(){
		return false;
	}

private:
	unsigned int element_;
	static inline unsigned int characteristic_;
	static inline std::vector<unsigned int> inverse_;

	static unsigned int _add(unsigned int element, unsigned int v);
	static unsigned int _substract(unsigned int element, unsigned int v);
	static unsigned int _multiply(unsigned int element, unsigned int v);
};

//unsigned int Shared_Zp_field_element::characteristic_;
//std::vector<unsigned int> Shared_Zp_field_element::inverse_;

inline Shared_Zp_field_element::Shared_Zp_field_element()
	: element_(0)
{}

inline Shared_Zp_field_element::Shared_Zp_field_element(unsigned int element)
	: element_(element < characteristic_ ? element : element % characteristic_)
{
//	if (characteristic_ != 0)
//		element_ %= characteristic_;
}

inline Shared_Zp_field_element::Shared_Zp_field_element(int element)
	: element_()
{
	int res = element < static_cast<int>(characteristic_) ? element : element % static_cast<int>(characteristic_);
	if (res < 0) res += characteristic_;
	element_ = res;
}

inline Shared_Zp_field_element::Shared_Zp_field_element(const Shared_Zp_field_element &toCopy)
	: element_(toCopy.element_)
{}

inline Shared_Zp_field_element::Shared_Zp_field_element(Shared_Zp_field_element &&toMove) noexcept
	: element_(std::exchange(toMove.element_, 0))
{}

inline void Shared_Zp_field_element::initialize(unsigned int characteristic)
{
//	std::cout << "charac: " << characteristic << " ";
	if (characteristic <= 1)
		throw std::invalid_argument("Characteristic must be strictly positive and a prime number.");

	inverse_.resize(characteristic);
	inverse_[0] = 0;
	for (unsigned int i = 1; i < characteristic; ++i) {
		unsigned int inv = 1;
		unsigned int mult = inv * i;
		while ((mult % characteristic) != 1) {
			++inv;
			if (mult == characteristic)
				throw std::invalid_argument("Characteristic must be a prime number.");
			mult = inv * i;
		}
		inverse_[i] = inv;
	}

	characteristic_ = characteristic;
//	std::cout << characteristic_ << "\n";
}

//inline Shared_Zp_field_element &Shared_Zp_field_element::operator+=(Shared_Zp_field_element const &f)
//{
//	element_ = _add(element_, f.element_);
//	return *this;
//}

//inline Shared_Zp_field_element &Shared_Zp_field_element::operator+=(unsigned int const v)
//{
//	element_ = _add(element_, v % characteristic_);
//	return *this;
//}

//inline Shared_Zp_field_element &Shared_Zp_field_element::operator-=(Shared_Zp_field_element const &f)
//{
//	element_ = _substract(element_, f.element_);
//	return *this;
//}

//inline Shared_Zp_field_element &Shared_Zp_field_element::operator-=(unsigned int const v)
//{
//	element_ = _substract(element_, v % characteristic_);
//	return *this;
//}

//inline Shared_Zp_field_element &Shared_Zp_field_element::operator*=(Shared_Zp_field_element const &f)
//{
//	element_ = _multiply(element_, f.element_);
//	return *this;
//}

//inline Shared_Zp_field_element &Shared_Zp_field_element::operator*=(unsigned int const v)
//{
//	element_ = _multiply(element_, v % characteristic_);
//	return *this;
//}

inline Shared_Zp_field_element &Shared_Zp_field_element::operator=(Shared_Zp_field_element other)
{
	std::swap(element_, other.element_);
	return *this;
}

inline Shared_Zp_field_element &Shared_Zp_field_element::operator=(unsigned int const value)
{
	element_ = value < characteristic_ ? value : value % characteristic_;
	return *this;
}

inline Shared_Zp_field_element::operator unsigned int() const
{
	return element_;
}

inline Shared_Zp_field_element Shared_Zp_field_element::get_inverse() const
{
	return Shared_Zp_field_element(inverse_[element_]);
}

inline std::pair<Shared_Zp_field_element, unsigned int>
Shared_Zp_field_element::get_partial_inverse(unsigned int product_of_characteristics) const
{
	return {get_inverse(), product_of_characteristics};
}

inline Shared_Zp_field_element Shared_Zp_field_element::get_additive_identity()
{
	return Shared_Zp_field_element();
}

inline Shared_Zp_field_element Shared_Zp_field_element::get_multiplicative_identity()
{
	return Shared_Zp_field_element(1);
}

inline Shared_Zp_field_element Shared_Zp_field_element::get_partial_multiplicative_identity()
{
	return Shared_Zp_field_element(1);
}

inline unsigned int Shared_Zp_field_element::get_characteristic()
{
	return characteristic_;
}

inline unsigned int Shared_Zp_field_element::get_value() const
{
	return element_;
}

inline unsigned int Shared_Zp_field_element::_add(unsigned int element, unsigned int v)
{
	if (UINT_MAX - element < v) {
		//automatic unsigned integer overflow behaviour will make it work
		element += v;
		element -= characteristic_;
		return element;
	}

	element += v;
	if (element >= characteristic_) element -= characteristic_;

	return element;
}

inline unsigned int Shared_Zp_field_element::_substract(unsigned int element, unsigned int v)
{
	if (element < v){
		element += characteristic_;
	}
	element -= v;

	return element;
}

inline unsigned int Shared_Zp_field_element::_multiply(unsigned int element, unsigned int v)
{
	unsigned int a = element;
	element = 0;
	unsigned int temp_b;

	while (a != 0) {
		if (a & 1) {
			if (v >= characteristic_ - element)
				element -= characteristic_;
			element += v;
		}
		a >>= 1;

		temp_b = v;
		if (v >= characteristic_ - v)
			temp_b -= characteristic_;
		v += temp_b;
	}

	return element;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif  // MATRIX_FIELD_ZP_VAR_H_
