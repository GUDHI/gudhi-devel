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
	Zp_field_element(const Zp_field_element& toCopy);
	Zp_field_element(Zp_field_element&& toMove) noexcept;

	Zp_field_element& operator+=(Zp_field_element const &f);
	friend Zp_field_element operator+(Zp_field_element f1, Zp_field_element const& f2){
		f1 += f2;
		return f1;
	}
	Zp_field_element& operator+=(unsigned int const v);
	friend Zp_field_element operator+(Zp_field_element f, unsigned int const v){
		f += v;
		return f;
	}
	friend unsigned int operator+(unsigned int v, Zp_field_element const& f){
		v += f.element_;
		v %= characteristic;
		return v;
	}

	Zp_field_element& operator-=(Zp_field_element const &f);
	friend Zp_field_element operator-(Zp_field_element f1, Zp_field_element const& f2){
		f1 -= f2;
		return f1;
	}
	Zp_field_element& operator-=(unsigned int const v);
	friend Zp_field_element operator-(Zp_field_element f, unsigned int const v){
		f -= v;
		return f;
	}
	friend unsigned int operator-(unsigned int v, Zp_field_element const& f){
		if (v >= characteristic) v %= characteristic;
		if (f.element_ > v) v += characteristic;
		v -= f.element_;
		return v;
	}

	Zp_field_element& operator*=(Zp_field_element const &f);
	friend Zp_field_element operator*(Zp_field_element f1, Zp_field_element const& f2){
		f1 *= f2;
		return f1;
	}
	Zp_field_element& operator*=(unsigned int const v);
	friend Zp_field_element operator*(Zp_field_element f, unsigned int const v){
		f *= v;
		return f;
	}
	friend unsigned int operator*(unsigned int v, Zp_field_element const& f){
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
	friend bool operator==(const unsigned int v, const Zp_field_element& f){
		if (v < characteristic) return v == f.element_;
		return (v % characteristic) == f.element_;
	}
	friend bool operator==(const Zp_field_element& f, const unsigned int v){
		if (v < characteristic) return v == f.element_;
		return (v % characteristic) == f.element_;
	}

	Zp_field_element& operator=(Zp_field_element other);
	Zp_field_element& operator=(const unsigned int value);
	operator unsigned int() const;
	friend void swap(Zp_field_element& f1, Zp_field_element& f2){
		std::swap(f1.element_, f2.element_);
	}

	Zp_field_element get_inverse() const;

	static Zp_field_element get_additive_identity();
	static Zp_field_element get_multiplicative_identity();
	static constexpr int get_characteristic();

	unsigned int get_value() const;

private:
	unsigned int element_;
	static inline std::array<unsigned int,characteristic> inverse_;

	void _add(unsigned int v);
	void _substract(unsigned int v);
	void _multiply(unsigned int v);

	static constexpr bool _is_prime(const int p);
};

template<unsigned int characteristic>
inline Zp_field_element<characteristic>::Zp_field_element()
	: element_(0)
{
	static_assert(_is_prime(characteristic), "Characteristic has to be a prime number.");
}

template<unsigned int characteristic>
inline Zp_field_element<characteristic>::Zp_field_element(unsigned int element)
	: element_(element % characteristic)
{
	static_assert(_is_prime(characteristic), "Characteristic has to be a prime number.");
}

template<unsigned int characteristic>
inline Zp_field_element<characteristic>::Zp_field_element(const Zp_field_element<characteristic> &toCopy)
	: element_(toCopy.element_)
{
	static_assert(_is_prime(characteristic), "Characteristic has to be a prime number.");
}

template<unsigned int characteristic>
inline Zp_field_element<characteristic>::Zp_field_element(Zp_field_element<characteristic> &&toMove) noexcept
	: element_(std::exchange(toMove.element_, 0))
{
	static_assert(_is_prime(characteristic), "Characteristic has to be a prime number.");
}

template<unsigned int characteristic>
inline Zp_field_element<characteristic> &Zp_field_element<characteristic>::operator+=(Zp_field_element<characteristic> const &f)
{
	_add(f.element_);
	return *this;
}

template<unsigned int characteristic>
inline Zp_field_element<characteristic> &Zp_field_element<characteristic>::operator+=(unsigned int const v)
{
	_add(v % characteristic);
	return *this;
}

template<unsigned int characteristic>
inline Zp_field_element<characteristic> &Zp_field_element<characteristic>::operator-=(Zp_field_element<characteristic> const &f)
{
	_substract(f.element_);
	return *this;
}

template<unsigned int characteristic>
inline Zp_field_element<characteristic> &Zp_field_element<characteristic>::operator-=(unsigned int const v)
{
	_substract(v % characteristic);
	return *this;
}

template<unsigned int characteristic>
inline Zp_field_element<characteristic> &Zp_field_element<characteristic>::operator*=(Zp_field_element<characteristic> const &f)
{
	_multiply(f.element_);
	return *this;
}

template<unsigned int characteristic>
inline Zp_field_element<characteristic> &Zp_field_element<characteristic>::operator*=(unsigned int const v)
{
	_multiply(v % characteristic);
	return *this;
}

template<unsigned  characteristic>
inline Zp_field_element<characteristic> &Zp_field_element<characteristic>::operator=(Zp_field_element other)
{
	std::swap(element_, other.element_);
	return *this;
}

template<unsigned  characteristic>
inline Zp_field_element<characteristic> &Zp_field_element<characteristic>::operator=(unsigned int const value)
{
	element_ = value % characteristic;
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
inline void Zp_field_element<characteristic>::_add(unsigned int v)
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
inline void Zp_field_element<characteristic>::_substract(unsigned int v)
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
inline void Zp_field_element<characteristic>::_multiply(unsigned int v)
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

template<unsigned int characteristic>
inline constexpr bool Zp_field_element<characteristic>::_is_prime(const int p)
{
	if (p <= 1) return false;
	if (p <= 3) return true;
	if (p % 2 == 0 || p % 3 == 0) return false;

	for (long i = 5; i * i <= p; i = i + 6)
		if (p % i == 0 || p % (i + 2) == 0)
			return false;

	return true;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif  // MATRIX_FIELD_ZP_H_
