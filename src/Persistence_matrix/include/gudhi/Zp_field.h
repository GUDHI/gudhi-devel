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
	Zp_field_element();
	Zp_field_element(unsigned int element);
	Zp_field_element(Zp_field_element& toCopy);
	Zp_field_element(Zp_field_element&& toMove) noexcept;

	Zp_field_element& operator+=(Zp_field_element const &f);
	template<unsigned int friendCharacteristic>
	friend Zp_field_element<friendCharacteristic> operator+(Zp_field_element<friendCharacteristic> f1, Zp_field_element<friendCharacteristic> const& f2);
	Zp_field_element& operator+=(unsigned int const &v);
	template<unsigned int friendCharacteristic>
	friend Zp_field_element<friendCharacteristic> operator+(Zp_field_element<friendCharacteristic> f, unsigned int const& v);
	template<unsigned int friendCharacteristic>
	friend unsigned int operator+(unsigned int v, Zp_field_element<friendCharacteristic> const& f);

	Zp_field_element& operator-=(Zp_field_element const &f);
	template<unsigned int friendCharacteristic>
	friend Zp_field_element<friendCharacteristic> operator-(Zp_field_element<friendCharacteristic> f1, Zp_field_element<friendCharacteristic> const& f2);
	Zp_field_element& operator-=(unsigned int const &v);
	template<unsigned int friendCharacteristic>
	friend Zp_field_element<friendCharacteristic> operator-(Zp_field_element<friendCharacteristic> f, unsigned int const& v);
	template<unsigned int friendCharacteristic>
	friend unsigned int operator-(unsigned int v, Zp_field_element<friendCharacteristic> const& f);

	Zp_field_element& operator*=(Zp_field_element const &f);
	template<unsigned int friendCharacteristic>
	friend Zp_field_element<friendCharacteristic> operator*(Zp_field_element<friendCharacteristic> f1, Zp_field_element<friendCharacteristic> const& f2);
	Zp_field_element& operator*=(unsigned int const &v);
	template<unsigned int friendCharacteristic>
	friend Zp_field_element<friendCharacteristic> operator*(Zp_field_element<friendCharacteristic> f, unsigned int const& v);
	template<unsigned int friendCharacteristic>
	friend unsigned int operator*(unsigned int const& v, Zp_field_element<friendCharacteristic> const& f);

	template<unsigned int friendCharacteristic>
	friend bool operator==(const Zp_field_element<friendCharacteristic>& f1, const Zp_field_element<friendCharacteristic>& f2);
	template<unsigned int friendCharacteristic>
	friend bool operator==(const unsigned int& v, const Zp_field_element<friendCharacteristic>& f);
	template<unsigned int friendCharacteristic>
	friend bool operator==(const Zp_field_element<friendCharacteristic>& f, const unsigned int& v);
	Zp_field_element& operator=(Zp_field_element other);
	operator unsigned int() const;
	template<unsigned int friendCharacteristic>
	friend void swap(Zp_field_element<friendCharacteristic>& f1, Zp_field_element<friendCharacteristic>& f2);

	Zp_field_element get_inverse() const;

	static Zp_field_element get_additive_identity();
	static Zp_field_element get_multiplicative_identity();
	static constexpr int get_characteristic();

	unsigned int get_value() const;

private:
	unsigned int element_;
	static std::array<unsigned int,characteristic> inverse_;

	void add_(unsigned int v);
	void substract_(unsigned int v);
	void multiply_(unsigned int v);
};

template<unsigned int characteristic>
inline Zp_field_element<characteristic>::Zp_field_element()
	: element_(0)
{}

template<unsigned int characteristic>
inline Zp_field_element<characteristic>::Zp_field_element(unsigned int element)
	: element_(element % characteristic)
{}

template<unsigned int characteristic>
inline Zp_field_element<characteristic>::Zp_field_element(Zp_field_element<characteristic> &toCopy)
	: element_(toCopy.element_)
{}

template<unsigned int characteristic>
inline Zp_field_element<characteristic>::Zp_field_element(Zp_field_element<characteristic> &&toMove) noexcept
	: element_(std::exchange(toMove.element_, 0))
{}

template<unsigned int characteristic>
inline Zp_field_element<characteristic> &Zp_field_element<characteristic>::operator+=(Zp_field_element<characteristic> const &f)
{
	add_(f.element_);
	return *this;
}

template<unsigned int characteristic>
inline Zp_field_element<characteristic> &Zp_field_element<characteristic>::operator+=(unsigned int const &v)
{
	add_(v % characteristic);
	return *this;
}

template<unsigned int characteristic>
inline Zp_field_element<characteristic> &Zp_field_element<characteristic>::operator-=(Zp_field_element<characteristic> const &f)
{
	substract_(f.element_);
	return *this;
}

template<unsigned int characteristic>
inline Zp_field_element<characteristic> &Zp_field_element<characteristic>::operator-=(unsigned int const &v)
{
	substract_(v % characteristic);
	return *this;
}

template<unsigned int characteristic>
inline Zp_field_element<characteristic> &Zp_field_element<characteristic>::operator*=(Zp_field_element<characteristic> const &f)
{
	multiply_(f.element_);
	return *this;
}

template<unsigned int characteristic>
inline Zp_field_element<characteristic> &Zp_field_element<characteristic>::operator*=(unsigned int const &v)
{
	multiply_(v % characteristic);
	return *this;
}

template<unsigned  characteristic>
inline Zp_field_element<characteristic> &Zp_field_element<characteristic>::operator=(Zp_field_element other)
{
	std::swap(element_, other.element_);
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
	if (element_ != 0 && inverse_[element_] == 0) {
		unsigned int inv = 1;
		unsigned int mult = inv * element_;
		while ((mult % characteristic) != 1) {
			++inv;
			if (mult == characteristic)
				throw std::invalid_argument("Characteristic must be a prime number.");
			mult = inv * element_;
		}
		inverse_[element_] = inv;
	}

	return Zp_field_element<characteristic>(inverse_[element_]);
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
inline constexpr int Zp_field_element<characteristic>::get_characteristic()
{
	return characteristic;
}

template<unsigned int characteristic>
inline unsigned int Zp_field_element<characteristic>::get_value() const
{
	return element_;
}

template<unsigned int characteristic>
inline void Zp_field_element<characteristic>::add_(unsigned int v)
{
	if (UINT_MAX - element_ < v) {
		unsigned long int sum = static_cast<unsigned long int>(element_) + v;
		element_ = sum - characteristic;
		return;
	}

	element_ += v;
	if (element_ >= characteristic) element_ -= characteristic;
}

template<unsigned int characteristic>
inline void Zp_field_element<characteristic>::substract_(unsigned int v)
{
	if (element_ < v){
		if (UINT_MAX - element_ < characteristic) {
			unsigned long int sum = static_cast<unsigned long int>(element_) + characteristic;
			element_ = sum - v;
			return;
		}
		element_ += characteristic;
	}
	element_ -= v;
}

template<unsigned int characteristic>
inline void Zp_field_element<characteristic>::multiply_(unsigned int v)
{
	unsigned int a = element_;
	element_ = 0;
	unsigned int temp_b;

	while (a != 0) {
		if (a & 1) {
			if (v >= characteristic - element_)
				element_ -= characteristic;
			element_ += v;
		}
		a >>= 1;

		temp_b = v;
		if (v >= characteristic - v)
			temp_b -= characteristic;
		v += temp_b;
	}
}

template<unsigned int friendCharacteristic>
Zp_field_element<friendCharacteristic> operator+(Zp_field_element<friendCharacteristic> f1, Zp_field_element<friendCharacteristic> const& f2)
{
	f1 += f2;
	return f1;
}

template<unsigned int friendCharacteristic>
Zp_field_element<friendCharacteristic> operator+(Zp_field_element<friendCharacteristic> f, unsigned int const& v)
{
	f += v;
	return f;
}

template<unsigned int friendCharacteristic>
unsigned int operator+(unsigned int v, Zp_field_element<friendCharacteristic> const& f)
{
	v += f.element_;
	v %= friendCharacteristic;
	return v;
}

template<unsigned int friendCharacteristic>
Zp_field_element<friendCharacteristic> operator-(Zp_field_element<friendCharacteristic> f1, Zp_field_element<friendCharacteristic> const& f2)
{
	f1 -= f2;
	return f1;
}

template<unsigned int friendCharacteristic>
Zp_field_element<friendCharacteristic> operator-(Zp_field_element<friendCharacteristic> f, unsigned int const& v)
{
	f -= v;
	return f;
}

template<unsigned int friendCharacteristic>
unsigned int operator-(unsigned int v, Zp_field_element<friendCharacteristic> const& f)
{
	if (f.element_ > v) v += friendCharacteristic;
	else if (v >= friendCharacteristic) v %= friendCharacteristic;
	v -= f.element_;
	return v;
}

template<unsigned int friendCharacteristic>
Zp_field_element<friendCharacteristic> operator*(Zp_field_element<friendCharacteristic> f1, Zp_field_element<friendCharacteristic> const& f2)
{
	f1 *= f2;
	return f1;
}

template<unsigned int friendCharacteristic>
Zp_field_element<friendCharacteristic> operator*(Zp_field_element<friendCharacteristic> f, unsigned int const& v)
{
	f *= v;
	return f;
}

template<unsigned int friendCharacteristic>
unsigned int operator*(unsigned int const& v, Zp_field_element<friendCharacteristic> const& f)
{
	unsigned int a = v;
	unsigned int b = f.element_;
	unsigned int res = 0;
	unsigned int temp_b;

	while (a != 0) {
		if (a & 1) {
			if (b >= friendCharacteristic - res)
				res -= friendCharacteristic;
			res += b;
		}
		a >>= 1;

		temp_b = b;
		if (b >= friendCharacteristic - b)
			temp_b -= friendCharacteristic;
		b += temp_b;
	}

	return res;
}

template<unsigned int friendCharacteristic>
bool operator==(const Zp_field_element<friendCharacteristic>& f1, const Zp_field_element<friendCharacteristic>& f2)
{
	return f1.element_ == f2.element_;
}

template<unsigned int friendCharacteristic>
bool operator==(const unsigned int& v, const Zp_field_element<friendCharacteristic>& f)
{
	if (v < friendCharacteristic) return v == f.element_;
	return (v % friendCharacteristic) == f.element_;
}

template<unsigned int friendCharacteristic>
bool operator==(const Zp_field_element<friendCharacteristic>& f, const unsigned int& v)
{
	if (v < friendCharacteristic) return v == f.element_;
	return (v % friendCharacteristic) == f.element_;
}

template<unsigned int friendCharacteristic>
void swap(Zp_field_element<friendCharacteristic>& f1, Zp_field_element<friendCharacteristic>& f2)
{
	std::swap(f1.element_, f2.element_);
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif  // MATRIX_FIELD_ZP_H_
