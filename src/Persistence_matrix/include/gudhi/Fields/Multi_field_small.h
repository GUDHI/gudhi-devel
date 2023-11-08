/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef MATRIX_FIELD_MULTI_SMALL_H_
#define MATRIX_FIELD_MULTI_SMALL_H_

#include <utility>
#include <vector>
#include <limits.h>
#include <iostream>
#include <gmpxx.h>
#include <cmath>
#include <numeric>

namespace Gudhi {
namespace persistence_matrix {

//productOfAllCharacteristics_ has to fit in an unsigned int, ie, productOfAllCharacteristics_ < 4294967295
template<unsigned int minimum, unsigned int maximum>
class Multi_field_element_with_small_characteristics {
public:
	using element_type = unsigned int;
	template <class T>
	using isInteger = std::enable_if_t<std::is_integral_v<T> >;

	Multi_field_element_with_small_characteristics();
	Multi_field_element_with_small_characteristics(unsigned int element);
	Multi_field_element_with_small_characteristics(int element);
	Multi_field_element_with_small_characteristics(const Multi_field_element_with_small_characteristics& toCopy);
	Multi_field_element_with_small_characteristics(Multi_field_element_with_small_characteristics&& toMove) noexcept;

	friend void operator+=(Multi_field_element_with_small_characteristics& f1, Multi_field_element_with_small_characteristics const& f2){
		f1.element_ = _add(f1.element_, f2.element_);
	}
	friend Multi_field_element_with_small_characteristics operator+(Multi_field_element_with_small_characteristics f1, Multi_field_element_with_small_characteristics const& f2){
		f1 += f2;
		return f1;
	}
	friend void operator+=(Multi_field_element_with_small_characteristics& f, unsigned int const v){
		f.element_ = _add(f.element_, v < productOfAllCharacteristics_ ? v : v % productOfAllCharacteristics_);
	}
	//v is assumed to be positive and will be casted into an unsigned int
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend Multi_field_element_with_small_characteristics operator+(Multi_field_element_with_small_characteristics f, const Integer_type v){
		f += v;
		return f;
	}
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend Integer_type operator+(Integer_type v, Multi_field_element_with_small_characteristics const& f){
		v += f.element_;
		v %= productOfAllCharacteristics_;
		return v;
	}

	friend void operator-=(Multi_field_element_with_small_characteristics& f1, Multi_field_element_with_small_characteristics const& f2){
		f1.element_ = _substract(f1.element_, f2.element_);
	}
	friend Multi_field_element_with_small_characteristics operator-(Multi_field_element_with_small_characteristics f1, Multi_field_element_with_small_characteristics const& f2){
		f1 -= f2;
		return f1;
	}
	friend void operator-=(Multi_field_element_with_small_characteristics& f, unsigned int const v){
		f.element_ = _substract(f.element_, v < productOfAllCharacteristics_ ? v : v % productOfAllCharacteristics_);
	}
	//v is assumed to be positive and will be casted into an unsigned int
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend Multi_field_element_with_small_characteristics operator-(Multi_field_element_with_small_characteristics f, const Integer_type v){
		f -= v;
		return f;
	}
	//v is assumed to be positive
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend Integer_type operator-(Integer_type v, Multi_field_element_with_small_characteristics const& f){
		if (v >= productOfAllCharacteristics_) v %= productOfAllCharacteristics_;
		if (f.element_ > v) v += productOfAllCharacteristics_;
		v -= f.element_;
		return v;
	}

	friend void operator*=(Multi_field_element_with_small_characteristics& f1, Multi_field_element_with_small_characteristics const& f2){
		f1.element_ = _multiply(f1.element_, f2.element_);
	}
	friend Multi_field_element_with_small_characteristics operator*(Multi_field_element_with_small_characteristics f1, Multi_field_element_with_small_characteristics const& f2){
		f1 *= f2;
		return f1;
	}
	friend void operator*=(Multi_field_element_with_small_characteristics& f, unsigned int const v){
		f.element_ = _multiply(f.element_, v < productOfAllCharacteristics_ ? v : v % productOfAllCharacteristics_);
	}
	//v is assumed to be positive and will be casted into an unsigned int
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend Multi_field_element_with_small_characteristics operator*(Multi_field_element_with_small_characteristics f, const Integer_type v){
		f *= v;
		return f;
	}
	//uses bitwise operations on v, so be carefull with signed integers
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend Integer_type operator*(Integer_type v, Multi_field_element_with_small_characteristics const& f){
		unsigned int b = f.element_;
		unsigned int res = 0;
		unsigned int temp_b;

		while (v != 0) {
			if (v & 1) {
				if (b >= productOfAllCharacteristics_ - res)
					res -= productOfAllCharacteristics_;
				res += b;
			}
			v >>= 1;

			temp_b = b;
			if (b >= productOfAllCharacteristics_ - b)
				temp_b -= productOfAllCharacteristics_;
			b += temp_b;
		}

		return res;
	}

	friend bool operator==(const Multi_field_element_with_small_characteristics& f1, const Multi_field_element_with_small_characteristics& f2){
		return f1.element_ == f2.element_;
	}
	//v is assumed to be positive
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend bool operator==(const Integer_type v, const Multi_field_element_with_small_characteristics& f){
		if (v < productOfAllCharacteristics_) return v == f.element_;
		return (v % productOfAllCharacteristics_) == f.element_;
	}
	//v is assumed to be positive
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend bool operator==(const Multi_field_element_with_small_characteristics& f, const Integer_type v){
		if (v < productOfAllCharacteristics_) return v == f.element_;
		return (v % productOfAllCharacteristics_) == f.element_;
	}
	friend bool operator!=(const Multi_field_element_with_small_characteristics& f1, const Multi_field_element_with_small_characteristics& f2){
		return !(f1 == f2);
	}
	//v is assumed to be positive
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend bool operator!=(const Integer_type v, const Multi_field_element_with_small_characteristics& f){
		return !(v == f);
	}
	//v is assumed to be positive
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend bool operator!=(const Multi_field_element_with_small_characteristics& f, const Integer_type v){
		return !(v == f);
	}

	Multi_field_element_with_small_characteristics& operator=(Multi_field_element_with_small_characteristics other);
	Multi_field_element_with_small_characteristics& operator=(const unsigned int value);
	operator unsigned int() const;
	friend void swap(Multi_field_element_with_small_characteristics& f1, Multi_field_element_with_small_characteristics& f2){
		std::swap(f1.element_, f2.element_);
	}

	Multi_field_element_with_small_characteristics get_inverse() const;
	std::pair<Multi_field_element_with_small_characteristics, unsigned int> get_partial_inverse(unsigned int productOfCharacteristics) const;

	static Multi_field_element_with_small_characteristics get_additive_identity();
	static Multi_field_element_with_small_characteristics get_multiplicative_identity();
	Multi_field_element_with_small_characteristics get_partial_multiplicative_identity();
	static constexpr unsigned int get_characteristic();

	unsigned int get_value() const;

	static constexpr bool handles_only_z2(){
		return false;
	}

private:
	static constexpr bool _is_prime(const int p);
	static constexpr unsigned int _multiply(unsigned int a, unsigned int b);
	static constexpr unsigned int _add(unsigned int element, unsigned int v);
	static constexpr unsigned int _substract(unsigned int element, unsigned int v);
	static constexpr int _get_inverse(unsigned int element, const unsigned int mod);

	unsigned int element_;
	static inline const std::vector<unsigned int> primes_ = [](){
		std::vector<unsigned int> res;
		for (unsigned int i = minimum; i <= maximum; ++i){
			if (_is_prime(i)){
				res.push_back(i);
			}
		}
		return res;
	}();
	static inline constexpr unsigned int productOfAllCharacteristics_ = [](){
		unsigned int res = 1;
		for (unsigned int i = minimum; i <= maximum; ++i){
			if (_is_prime(i)){
				res *= i;
			}
		}
		return res;
	}();
	static inline const std::vector<unsigned int> partials_ = [](){
		std::vector<unsigned int> res;

		if (productOfAllCharacteristics_ == 1) return res;

		for (unsigned int i = 0; i < primes_.size(); ++i){
			unsigned int p = primes_[i];
			unsigned int base = productOfAllCharacteristics_ / p;
			unsigned int exp = p - 1;
			res.push_back(1);

			while (exp > 0) {
				// If exp is odd, multiply with result
				if (exp & 1)
					res.back() = _multiply(res.back(), base);
				// y must be even now
				exp = exp >> 1;
				base = _multiply(base, base);
			}
		}

		return res;
	}();
	//If I understood the paper well, multiplicativeID_ always equals to 1. But in Clement's code,
	//multiplicativeID_ is computed (see commented lambda function below). TODO: verify with Clement.
	static inline constexpr unsigned int multiplicativeID_ = 1; /*= [](){
		unsigned int res = 0;
		for (unsigned int i = 0; i < partials_.size(); ++i){
			res = (res + partials_[i]) % productOfAllCharacteristics_;
		}

		return res;
	}();*/
};

template<unsigned int minimum, unsigned int maximum>
inline Multi_field_element_with_small_characteristics<minimum,maximum>::Multi_field_element_with_small_characteristics()
	: element_(0)
{
	static_assert(maximum >= 2, "Characteristics have to be positive.");
	static_assert(minimum <= maximum, "The given interval is not valid.");
	static_assert(minimum != maximum || _is_prime(minimum), "The given interval does not contain a prime number.");
	static_assert(productOfAllCharacteristics_ != 1, "The given interval does not contain a prime number.");
}

template<unsigned int minimum, unsigned int maximum>
inline Multi_field_element_with_small_characteristics<minimum,maximum>::Multi_field_element_with_small_characteristics(unsigned int element)
	: element_(element % productOfAllCharacteristics_)
{
	static_assert(maximum >= 2, "Characteristics has to be positive.");
	static_assert(minimum <= maximum, "The given interval is not valid.");
	static_assert(minimum != maximum || _is_prime(minimum), "The given interval does not contain a prime number.");
	static_assert(productOfAllCharacteristics_ != 1, "The given interval does not contain a prime number.");
}

template<unsigned int minimum, unsigned int maximum>
inline Multi_field_element_with_small_characteristics<minimum,maximum>::Multi_field_element_with_small_characteristics(int element)
	: element_(element % productOfAllCharacteristics_)
{
	static_assert(maximum >= 2, "Characteristics has to be positive.");
	static_assert(minimum <= maximum, "The given interval is not valid.");
	static_assert(minimum != maximum || _is_prime(minimum), "The given interval does not contain a prime number.");
	static_assert(productOfAllCharacteristics_ != 1, "The given interval does not contain a prime number.");
}

template<unsigned int minimum, unsigned int maximum>
inline Multi_field_element_with_small_characteristics<minimum,maximum>::Multi_field_element_with_small_characteristics(const Multi_field_element_with_small_characteristics<minimum,maximum> &toCopy)
	: element_(toCopy.element_)
{}

template<unsigned int minimum, unsigned int maximum>
inline Multi_field_element_with_small_characteristics<minimum,maximum>::Multi_field_element_with_small_characteristics(Multi_field_element_with_small_characteristics<minimum,maximum> &&toMove) noexcept
	: element_(std::exchange(toMove.element_, 0))
{}

//template<unsigned int minimum, unsigned int maximum>
//inline Multi_field_element_with_small_characteristics<minimum,maximum> &Multi_field_element_with_small_characteristics<minimum,maximum>::operator+=(Multi_field_element_with_small_characteristics<minimum,maximum> const &f)
//{
//	_add(f.element_);
//	return *this;
//}

//template<unsigned int minimum, unsigned int maximum>
//inline Multi_field_element_with_small_characteristics<minimum,maximum> &Multi_field_element_with_small_characteristics<minimum,maximum>::operator+=(unsigned int const v)
//{
//	_add(v % productOfAllCharacteristics_);
//	return *this;
//}

//template<unsigned int minimum, unsigned int maximum>
//inline Multi_field_element_with_small_characteristics<minimum,maximum> &Multi_field_element_with_small_characteristics<minimum,maximum>::operator-=(Multi_field_element_with_small_characteristics<minimum,maximum> const &f)
//{
//	_substract(f.element_);
//	return *this;
//}

//template<unsigned int minimum, unsigned int maximum>
//inline Multi_field_element_with_small_characteristics<minimum,maximum> &Multi_field_element_with_small_characteristics<minimum,maximum>::operator-=(unsigned int const v)
//{
//	_substract(v % productOfAllCharacteristics_);
//	return *this;
//}

//template<unsigned int minimum, unsigned int maximum>
//inline Multi_field_element_with_small_characteristics<minimum,maximum> &Multi_field_element_with_small_characteristics<minimum,maximum>::operator*=(Multi_field_element_with_small_characteristics<minimum,maximum> const &f)
//{
//	element_ = _multiply(element_, f.element_);
//	return *this;
//}

//template<unsigned int minimum, unsigned int maximum>
//inline Multi_field_element_with_small_characteristics<minimum,maximum> &Multi_field_element_with_small_characteristics<minimum,maximum>::operator*=(unsigned int const v)
//{
//	element_ = _multiply(element_, v % productOfAllCharacteristics_);
//	return *this;
//}

template<unsigned int minimum, unsigned int maximum>
inline Multi_field_element_with_small_characteristics<minimum,maximum> &Multi_field_element_with_small_characteristics<minimum,maximum>::operator=(Multi_field_element_with_small_characteristics other)
{
	std::swap(element_, other.element_);
	return *this;
}

template<unsigned int minimum, unsigned int maximum>
inline Multi_field_element_with_small_characteristics<minimum,maximum> &Multi_field_element_with_small_characteristics<minimum,maximum>::operator=(unsigned int const value)
{
	element_ = value % productOfAllCharacteristics_;
	return *this;
}

template<unsigned int minimum, unsigned int maximum>
inline Multi_field_element_with_small_characteristics<minimum,maximum>::operator unsigned int() const
{
	return element_;
}

template<unsigned int minimum, unsigned int maximum>
inline Multi_field_element_with_small_characteristics<minimum,maximum> Multi_field_element_with_small_characteristics<minimum,maximum>::get_inverse() const
{
	return get_partial_inverse(productOfAllCharacteristics_).first;
}

template<unsigned int minimum, unsigned int maximum>
inline std::pair<Multi_field_element_with_small_characteristics<minimum,maximum>, unsigned int> Multi_field_element_with_small_characteristics<minimum,maximum>::get_partial_inverse(unsigned int productOfCharacteristics) const
{
	unsigned int gcd = std::gcd(element_, productOfAllCharacteristics_);

	if (gcd == productOfCharacteristics)
		return {Multi_field_element_with_small_characteristics(), multiplicativeID_};  // partial inverse is 0

	unsigned int QT = productOfCharacteristics / gcd;
	Multi_field_element_with_small_characteristics res(QT);

	const unsigned int inv_qt = _get_inverse(element_, QT);

	res = res.get_partial_multiplicative_identity();
	res *= inv_qt;

	return {res, QT};
}

template<unsigned int minimum, unsigned int maximum>
inline Multi_field_element_with_small_characteristics<minimum,maximum> Multi_field_element_with_small_characteristics<minimum,maximum>::get_additive_identity()
{
	return Multi_field_element_with_small_characteristics<minimum,maximum>();
}

template<unsigned int minimum, unsigned int maximum>
inline Multi_field_element_with_small_characteristics<minimum,maximum> Multi_field_element_with_small_characteristics<minimum,maximum>::get_multiplicative_identity()
{
	return Multi_field_element_with_small_characteristics<minimum,maximum>(multiplicativeID_);
}

template<unsigned int minimum, unsigned int maximum>
inline Multi_field_element_with_small_characteristics<minimum,maximum> Multi_field_element_with_small_characteristics<minimum,maximum>::get_partial_multiplicative_identity()
{
	if (element_ == 0) {
		return Multi_field_element_with_small_characteristics<minimum,maximum>(multiplicativeID_);
	}
	Multi_field_element_with_small_characteristics<minimum,maximum> mult;
	for (unsigned int idx = 0; idx < primes_.size(); ++idx) {
		if ((element_ % primes_[idx]) == 0) {
			mult += partials_[idx];
		}
	}
	return mult;
}

template<unsigned int minimum, unsigned int maximum>
inline constexpr unsigned int Multi_field_element_with_small_characteristics<minimum, maximum>::get_characteristic()
{
	return productOfAllCharacteristics_;
}

template<unsigned int minimum, unsigned int maximum>
inline unsigned int Multi_field_element_with_small_characteristics<minimum,maximum>::get_value() const
{
	return element_;
}

template<unsigned int minimum, unsigned int maximum>
inline constexpr unsigned int Multi_field_element_with_small_characteristics<minimum,maximum>::_add(unsigned int element, unsigned int v)
{
	if (UINT_MAX - element < v) {
		//automatic unsigned integer overflow behaviour will make it work
		element += v;
		element -= productOfAllCharacteristics_;
		return element;
	}

	element += v;
	if (element >= productOfAllCharacteristics_) element -= productOfAllCharacteristics_;

	return element;
}

template<unsigned int minimum, unsigned int maximum>
inline constexpr unsigned int Multi_field_element_with_small_characteristics<minimum,maximum>::_substract(unsigned int element, unsigned int v)
{
	if (element < v){
		element += productOfAllCharacteristics_;
	}
	element -= v;

	return element;
}

template<unsigned int minimum, unsigned int maximum>
inline constexpr unsigned int Multi_field_element_with_small_characteristics<minimum,maximum>::_multiply(unsigned int a, unsigned int b)
{
	unsigned int res = 0;
	unsigned int temp_b = 0;

	if (b < a) std::swap(a, b);

	while (a != 0) {
		if (a & 1) {
			/* Add b to res, modulo m, without overflow */
			if (b >= productOfAllCharacteristics_ - res)
				res -= productOfAllCharacteristics_;
			res += b;
		}
		a >>= 1;

		/* Double b, modulo m */
		temp_b = b;
		if (b >= productOfAllCharacteristics_ - b)
			temp_b -= productOfAllCharacteristics_;
		b += temp_b;
	}
	return res;
}

template<unsigned int minimum, unsigned int maximum>
inline constexpr int Multi_field_element_with_small_characteristics<minimum,maximum>::_get_inverse(unsigned int element, const unsigned int mod)
{
	//to solve: Ax + My = 1
	int M = mod;
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
		x += mod;

	return x;
}

template<unsigned int minimum, unsigned int maximum>
inline constexpr bool Multi_field_element_with_small_characteristics<minimum,maximum>::_is_prime(const int p)
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

#endif  // MATRIX_FIELD_MULTI_SMALL_H_
