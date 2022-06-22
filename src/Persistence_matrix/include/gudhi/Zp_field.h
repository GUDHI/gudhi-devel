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

template<unsigned int characteristic>
class Zp_field_element {
public:
	Zp_field_element();
	Zp_field_element(unsigned int element);
	Zp_field_element(Zp_field_element& toCopy);
	Zp_field_element(Zp_field_element&& toMove) noexcept;

	Zp_field_element& operator+=(Zp_field_element const &f);
	friend Zp_field_element operator+(Zp_field_element f1, Zp_field_element const& f2);

	Zp_field_element& operator-=(Zp_field_element const &f);
	friend Zp_field_element operator-(Zp_field_element f1, Zp_field_element const& f2);

	Zp_field_element& operator*=(Zp_field_element const &f);
	friend Zp_field_element operator*(Zp_field_element f1, Zp_field_element const& f2);

	Zp_field_element& operator=(Zp_field_element other);
	template<unsigned int friendCharacteristic>
	friend void swap(Zp_field_element& f1, Zp_field_element& f2);

	Zp_field_element get_inverse() const;

	static Zp_field_element get_additive_identity();
	static Zp_field_element get_multiplicative_identity();
	static constexpr int get_characteristic();

	unsigned int get_value() const;

private:
	unsigned int element_;
	static std::array<unsigned int,characteristic> inverse_;
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
inline Zp_field_element<characteristic> &Zp_field_element<characteristic>::operator+=(Zp_field_element const &f)
{
	if (UINT_MAX - element_ < f.element_) {
		unsigned long int sum = static_cast<unsigned long int>(element_) + f.element_;
		element_ = sum - characteristic;
		return *this;
	}

	element_ += f.element_;
	if (element_ >= characteristic) element_ -= characteristic;
	return *this;
}

template<unsigned int characteristic>
inline Zp_field_element<characteristic> &Zp_field_element<characteristic>::operator-=(const Zp_field_element &f)
{
	if (element_ < f.element_){
		if (UINT_MAX - element_ < characteristic) {
			unsigned long int sum = static_cast<unsigned long int>(element_) + characteristic;
			element_ = sum - f.element_;
			return *this;
		}
		element_ += characteristic;
	}
	element_ -= f.element_;
	return *this;
}

template<unsigned int characteristic>
inline Zp_field_element<characteristic> &Zp_field_element<characteristic>::operator*=(const Zp_field_element<characteristic> &f)
{
	unsigned int a = element_;
	unsigned int b = f.element_;
	element_ = 0;
	unsigned int temp_b;

	while (a != 0) {
		if (a & 1) {
			if (b >= characteristic - element_)
				element_ -= characteristic;
			element_ += b;
		}
		a >>= 1;

		temp_b = b;
		if (b >= characteristic - b)
			temp_b -= characteristic;
		b += temp_b;
	}

	return *this;
}

template<unsigned  characteristic>
inline Zp_field_element<characteristic> &Zp_field_element<characteristic>::operator=(Zp_field_element other)
{
	std::swap(element_, other.element_);
	return *this;
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
Zp_field_element<characteristic> operator+(Zp_field_element<characteristic> f1, Zp_field_element<characteristic> const& f2)
{
	f1 += f2;
	return f1;
}

template<unsigned int characteristic>
Zp_field_element<characteristic> operator-(Zp_field_element<characteristic> f1, Zp_field_element<characteristic> const& f2)
{
	f1 -= f2;
	return f1;
}

template<unsigned int characteristic>
Zp_field_element<characteristic> operator*(Zp_field_element<characteristic> f1, Zp_field_element<characteristic> const& f2)
{
	f1 *= f2;
	return f1;
}

template<unsigned int friendCharacteristic>
void swap(Zp_field_element<friendCharacteristic>& f1, Zp_field_element<friendCharacteristic>& f2)
{
	std::swap(f1.element_, f2.element_);
}

#endif  // MATRIX_FIELD_ZP_H_
