/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef MATRIX_FIELD_ZP_H_
#define MATRIX_FIELD_ZP_H_

#include <utility>
#include <array>
#include <limits.h>
#include <iostream>

namespace Gudhi {
namespace persistence_matrix {

template<unsigned int characteristic>
class Zp_field_element {
public:
	using element_type = unsigned int;
	template <class T>
	using isInteger = std::enable_if_t<std::is_integral_v<T> >;

	Zp_field_element();
	Zp_field_element(unsigned int element);
	Zp_field_element(int element);
	Zp_field_element(const Zp_field_element& toCopy);
	Zp_field_element(Zp_field_element&& toMove) noexcept;

	friend void operator+=(Zp_field_element& f1, Zp_field_element const& f2){
		f1.element_ = Zp_field_element::_add(f1.element_, f2.element_);
	}
	friend Zp_field_element operator+(Zp_field_element f1, Zp_field_element const& f2){
		f1 += f2;
		return f1;
	}
	friend void operator+=(Zp_field_element& f, const unsigned int v){
		f.element_ = Zp_field_element::_add(f.element_, v < characteristic ? v : v % characteristic);
	}
	//v is assumed to be positive and will be casted into an unsigned int
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend Zp_field_element operator+(Zp_field_element f, const Integer_type v){
		f += v;
		return f;
	}
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend Integer_type operator+(Integer_type v, Zp_field_element const& f){
		v += f.element_;
		if (v >= characteristic) v %= characteristic;
		return v;
	}

	friend void operator-=(Zp_field_element& f1, Zp_field_element const& f2){
		f1.element_ = Zp_field_element::_substract(f1.element_, f2.element_);
	}
	friend Zp_field_element operator-(Zp_field_element f1, Zp_field_element const& f2){
		f1 -= f2;
		return f1;
	}
	friend void operator-=(Zp_field_element& f, unsigned int const v){
		f.element_ = Zp_field_element::_substract(f.element_, v < characteristic ? v : v % characteristic);
	}
	//v is assumed to be positive and will be casted into an unsigned int
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend Zp_field_element operator-(Zp_field_element f, const Integer_type v){
		f -= v;
		return f;
	}
	//v is assumed to be positive
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend unsigned int operator-(Integer_type v, Zp_field_element const& f){
		if (v >= characteristic) v %= characteristic;
		if (f.element_ > v) v += characteristic;
		v -= f.element_;
		return v;
	}

	friend void operator*=(Zp_field_element& f1, Zp_field_element const& f2){
		f1.element_ = Zp_field_element::_multiply(f1.element_, f2.element_);
	}
	friend Zp_field_element operator*(Zp_field_element f1, Zp_field_element const& f2){
		f1 *= f2;
		return f1;
	}
	friend void operator*=(Zp_field_element& f, unsigned int const v){
		f.element_ = Zp_field_element::_multiply(f.element_, v < characteristic ? v : v % characteristic);
	}
	//v is assumed to be positive and will be casted into an unsigned int
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend Zp_field_element operator*(Zp_field_element f, const Integer_type v){
		f *= v;
		return f;
	}
	//uses bitwise operations on v, so be carefull with signed integers
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend unsigned int operator*(Integer_type v, Zp_field_element const& f){
		unsigned int b = f.element_;
		unsigned int res = 0;
		unsigned int temp_b;

		while (v != 0) {
			if (v & 1) {
				if (b >= characteristic - res)
					res -= characteristic;
				res += b;
			}
			v >>= 1;

			temp_b = b;
			if (b >= characteristic - b)
				temp_b -= characteristic;
			b += temp_b;
		}

		return res;
	}

	friend bool operator==(const Zp_field_element& f1, const Zp_field_element& f2){
		return f1.element_ == f2.element_;
	}
	//v is assumed to be positive
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend bool operator==(const Integer_type v, const Zp_field_element& f){
		if (v < characteristic) return v == f.element_;
		return (v % characteristic) == f.element_;
	}
	//v is assumed to be positive
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend bool operator==(const Zp_field_element& f, const Integer_type v){
		if (v < characteristic) return v == f.element_;
		return (v % characteristic) == f.element_;
	}
	friend bool operator!=(const Zp_field_element& f1, const Zp_field_element& f2){
		return !(f1 == f2);
	}
	//v is assumed to be positive
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend bool operator!=(const Integer_type v, const Zp_field_element& f){
		return !(v == f);
	}
	//v is assumed to be positive
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend bool operator!=(const Zp_field_element& f, const Integer_type v){
		return !(v == f);
	}

	Zp_field_element& operator=(Zp_field_element other);
	Zp_field_element& operator=(const unsigned int value);
	operator unsigned int() const;
	friend void swap(Zp_field_element& f1, Zp_field_element& f2){
		std::swap(f1.element_, f2.element_);
	}

	Zp_field_element get_inverse() const;
	std::pair<Zp_field_element, unsigned int> get_partial_inverse(unsigned int product_of_characteristics) const;

	static Zp_field_element get_additive_identity();
	static Zp_field_element get_multiplicative_identity();
	static Zp_field_element get_partial_multiplicative_identity();
	static constexpr unsigned int get_characteristic();

	unsigned int get_value() const;

	static constexpr bool handles_only_z2(){
		return false;
	}

private:
	unsigned int element_;
	static inline std::array<unsigned int,characteristic> inverse_;

	static unsigned int _add(unsigned int element, unsigned int v);
	static unsigned int _substract(unsigned int element, unsigned int v);
	static unsigned int _multiply(unsigned int element, unsigned int v);
	static int _get_inverse(unsigned int element);

	static constexpr bool _is_prime();
};

template<unsigned int characteristic>
inline Zp_field_element<characteristic>::Zp_field_element()
	: element_(0)
{
	static_assert(_is_prime(), "Characteristic has to be a prime number.");
}

template<unsigned int characteristic>
inline Zp_field_element<characteristic>::Zp_field_element(unsigned int element)
	: element_(element < characteristic ? element : element % characteristic)
{
	static_assert(_is_prime(), "Characteristic has to be a prime number.");
}

template<unsigned int characteristic>
inline Zp_field_element<characteristic>::Zp_field_element(int element)
{
	int res = element < static_cast<int>(characteristic) ? element : element % static_cast<int>(characteristic);
	if (res < 0) res += characteristic;
	element_ = res;

	static_assert(_is_prime(), "Characteristic has to be a prime number.");
}

template<unsigned int characteristic>
inline Zp_field_element<characteristic>::Zp_field_element(const Zp_field_element<characteristic> &toCopy)
	: element_(toCopy.element_)
{}

template<unsigned int characteristic>
inline Zp_field_element<characteristic>::Zp_field_element(Zp_field_element<characteristic> &&toMove) noexcept
	: element_(std::exchange(toMove.element_, 0))
{}

template<unsigned int characteristic>
inline Zp_field_element<characteristic> &Zp_field_element<characteristic>::operator=(Zp_field_element other)
{
	std::swap(element_, other.element_);
	return *this;
}

template<unsigned int characteristic>
inline Zp_field_element<characteristic> &Zp_field_element<characteristic>::operator=(unsigned int const value)
{
	element_ = value < characteristic ? value : value % characteristic;
	return *this;
}

template<unsigned int characteristic>
inline Zp_field_element<characteristic>::operator unsigned int() const
{
	return element_;
}

template<unsigned int characteristic>
inline Zp_field_element<characteristic> Zp_field_element<characteristic>::get_inverse() const
{
	if (element_ != 0 && inverse_[element_] == 0) {		//initialize everything at instanciation instead?
		inverse_[element_] = _get_inverse(element_);
	}

	return Zp_field_element<characteristic>(inverse_[element_]);
}

template<unsigned int characteristic>
inline std::pair<Zp_field_element<characteristic>, unsigned int>
Zp_field_element<characteristic>::get_partial_inverse(unsigned int product_of_characteristics) const
{
	return {get_inverse(), product_of_characteristics};
}

template<unsigned int characteristic>
inline Zp_field_element<characteristic> Zp_field_element<characteristic>::get_additive_identity()
{
	return Zp_field_element<characteristic>();
}

template<unsigned int characteristic>
inline Zp_field_element<characteristic> Zp_field_element<characteristic>::get_multiplicative_identity()
{
	return Zp_field_element<characteristic>(1);
}

template<unsigned int characteristic>
inline Zp_field_element<characteristic> Zp_field_element<characteristic>::get_partial_multiplicative_identity()
{
	return Zp_field_element<characteristic>(1);
}

template<unsigned int characteristic>
inline constexpr unsigned int Zp_field_element<characteristic>::get_characteristic()
{
	return characteristic;
}

template<unsigned int characteristic>
inline unsigned int Zp_field_element<characteristic>::get_value() const
{
	return element_;
}

template<unsigned int characteristic>
inline unsigned int Zp_field_element<characteristic>::_add(unsigned int element, unsigned int v)
{
	if (UINT_MAX - element < v) {
		//automatic unsigned integer overflow behaviour will make it work
		element += v;
		element -= characteristic;
		return element;
	}

	element += v;
	if (element >= characteristic) element -= characteristic;

	return element;
}

template<unsigned int characteristic>
inline unsigned int Zp_field_element<characteristic>::_substract(unsigned int element, unsigned int v)
{
	if (element < v){
		element += characteristic;
	}
	element -= v;

	return element;
}

template<unsigned int characteristic>
inline unsigned int Zp_field_element<characteristic>::_multiply(unsigned int element, unsigned int v)
{
	unsigned int a = element;
	element = 0;
	unsigned int temp_b;

	while (a != 0) {
		if (a & 1) {
			if (v >= characteristic - element)
				element -= characteristic;
			element += v;
		}
		a >>= 1;

		temp_b = v;
		if (v >= characteristic - v)
			temp_b -= characteristic;
		v += temp_b;
	}

	return element;
}

template<unsigned int characteristic>
inline int Zp_field_element<characteristic>::_get_inverse(unsigned int element)
{
	//to solve: Ax + My = 1
	int M = characteristic;
	int A = element;
	int y = 0, x = 1;
	//extended euclidien division
	while (A > 1) {
		int quotient = A / M;
		int temp = M;

		M = A % M, A = temp;
		temp = y;

		y = x - quotient * y;
		x = temp;
	}

	if (x < 0)
		x += characteristic;

	return x;
}

template<unsigned int characteristic>
inline constexpr bool Zp_field_element<characteristic>::_is_prime()
{
	if (characteristic <= 1) return false;
	if (characteristic <= 3) return true;
	if (characteristic % 2 == 0 || characteristic % 3 == 0) return false;

	for (long i = 5; i * i <= characteristic; i = i + 6)
		if (characteristic % i == 0 || characteristic % (i + 2) == 0)
			return false;

	return true;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif  // MATRIX_FIELD_ZP_H_
