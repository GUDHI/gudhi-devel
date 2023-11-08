/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef MATRIX_FIELD_MULTI_SMALL_SHARED_H_
#define MATRIX_FIELD_MULTI_SMALL_SHARED_H_

#include <utility>
#include <vector>
#include <limits.h>
#include <iostream>
#include <gmpxx.h>
#include <cmath>
#include <numeric>

namespace Gudhi {
namespace persistence_matrix {

//productOfAllCharacteristics_ ^ 2 has to fit in an unsigned int
class Shared_multi_field_element_with_small_characteristics {
public:
	using element_type = unsigned int;
	template <class T>
	using isInteger = std::enable_if_t<std::is_integral_v<T> >;

	Shared_multi_field_element_with_small_characteristics();
	Shared_multi_field_element_with_small_characteristics(unsigned int element);
	Shared_multi_field_element_with_small_characteristics(int element);
	Shared_multi_field_element_with_small_characteristics(const Shared_multi_field_element_with_small_characteristics& toCopy);
	Shared_multi_field_element_with_small_characteristics(Shared_multi_field_element_with_small_characteristics&& toMove) noexcept;

	static void initialize(unsigned int minimum, unsigned int maximum);

	friend void operator+=(Shared_multi_field_element_with_small_characteristics& f1, Shared_multi_field_element_with_small_characteristics const& f2){
		f1.element_ = _add(f1.element_, f2.element_);
	}
	friend Shared_multi_field_element_with_small_characteristics operator+(Shared_multi_field_element_with_small_characteristics f1, Shared_multi_field_element_with_small_characteristics const& f2){
		f1 += f2;
		return f1;
	}
	friend void operator+=(Shared_multi_field_element_with_small_characteristics& f, unsigned int const v){
		f.element_ = _add(f.element_, v < productOfAllCharacteristics_ ? v : v % productOfAllCharacteristics_);
	}
	//v is assumed to be positive and will be casted into an unsigned int
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend Shared_multi_field_element_with_small_characteristics operator+(Shared_multi_field_element_with_small_characteristics f, const Integer_type v){
		f += v;
		return f;
	}
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend Integer_type operator+(Integer_type v, Shared_multi_field_element_with_small_characteristics const& f){
		v += f.element_;
		v %= productOfAllCharacteristics_;
		return v;
	}

	friend void operator-=(Shared_multi_field_element_with_small_characteristics& f1, Shared_multi_field_element_with_small_characteristics const& f2){
		f1.element_ = _substract(f1.element_, f2.element_);
	}
	friend Shared_multi_field_element_with_small_characteristics operator-(Shared_multi_field_element_with_small_characteristics f1, Shared_multi_field_element_with_small_characteristics const& f2){
		f1 -= f2;
		return f1;
	}
	friend void operator-=(Shared_multi_field_element_with_small_characteristics& f, unsigned int const v){
		f.element_ = _substract(f.element_, v < productOfAllCharacteristics_ ? v : v % productOfAllCharacteristics_);
	}
	//v is assumed to be positive and will be casted into an unsigned int
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend Shared_multi_field_element_with_small_characteristics operator-(Shared_multi_field_element_with_small_characteristics f, const Integer_type v){
		f -= v;
		return f;
	}
	//v is assumed to be positive
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend Integer_type operator-(Integer_type v, Shared_multi_field_element_with_small_characteristics const& f){
		if (v >= productOfAllCharacteristics_) v %= productOfAllCharacteristics_;
		if (f.element_ > v) v += productOfAllCharacteristics_;
		v -= f.element_;
		return v;
	}

	friend void operator*=(Shared_multi_field_element_with_small_characteristics& f1, Shared_multi_field_element_with_small_characteristics const& f2){
		f1.element_ = _multiply(f1.element_, f2.element_);
	}
	friend Shared_multi_field_element_with_small_characteristics operator*(Shared_multi_field_element_with_small_characteristics f1, Shared_multi_field_element_with_small_characteristics const& f2){
		f1 *= f2;
		return f1;
	}
	friend void operator*=(Shared_multi_field_element_with_small_characteristics& f, unsigned int const v){
		f.element_ = _multiply(f.element_, v < productOfAllCharacteristics_ ? v : v % productOfAllCharacteristics_);
	}
	//v is assumed to be positive and will be casted into an unsigned int
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend Shared_multi_field_element_with_small_characteristics operator*(Shared_multi_field_element_with_small_characteristics f, const Integer_type v){
		f *= v;
		return f;
	}
	//uses bitwise operations on v, so be carefull with signed integers
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend Integer_type operator*(Integer_type v, Shared_multi_field_element_with_small_characteristics const& f){
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

	friend bool operator==(const Shared_multi_field_element_with_small_characteristics& f1, const Shared_multi_field_element_with_small_characteristics& f2){
		return f1.element_ == f2.element_;
	}
	//v is assumed to be positive
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend bool operator==(const Integer_type v, const Shared_multi_field_element_with_small_characteristics& f){
		if (v < productOfAllCharacteristics_) return v == f.element_;
		return (v % productOfAllCharacteristics_) == f.element_;
	}
	//v is assumed to be positive
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend bool operator==(const Shared_multi_field_element_with_small_characteristics& f, const Integer_type v){
		if (v < productOfAllCharacteristics_) return v == f.element_;
		return (v % productOfAllCharacteristics_) == f.element_;
	}
	friend bool operator!=(const Shared_multi_field_element_with_small_characteristics& f1, const Shared_multi_field_element_with_small_characteristics& f2){
		return !(f1 == f2);
	}
	//v is assumed to be positive
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend bool operator!=(const Integer_type v, const Shared_multi_field_element_with_small_characteristics& f){
		return !(v == f);
	}
	//v is assumed to be positive
	template<typename Integer_type, class = isInteger<Integer_type> >
	friend bool operator!=(const Shared_multi_field_element_with_small_characteristics& f, const Integer_type v){
		return !(v == f);
	}

	Shared_multi_field_element_with_small_characteristics& operator=(Shared_multi_field_element_with_small_characteristics other);
	Shared_multi_field_element_with_small_characteristics& operator=(const unsigned int value);
	operator unsigned int() const;
	friend void swap(Shared_multi_field_element_with_small_characteristics& f1, Shared_multi_field_element_with_small_characteristics& f2){
		std::swap(f1.element_, f2.element_);
	}

	Shared_multi_field_element_with_small_characteristics get_inverse() const;
	std::pair<Shared_multi_field_element_with_small_characteristics, unsigned int> get_partial_inverse(unsigned int productOfCharacteristics) const;

	static Shared_multi_field_element_with_small_characteristics get_additive_identity();
	static Shared_multi_field_element_with_small_characteristics get_multiplicative_identity();
	Shared_multi_field_element_with_small_characteristics get_partial_multiplicative_identity();
	static unsigned int get_characteristic();

	unsigned int get_value() const;

	static constexpr bool handles_only_z2(){
		return false;
	}

private:
	static constexpr bool _is_prime(const int p);
	static unsigned int _multiply(unsigned int a, unsigned int b);
	static unsigned int _add(unsigned int element, unsigned int v);
	static unsigned int _substract(unsigned int element, unsigned int v);
	static constexpr int _get_inverse(unsigned int element, const unsigned int mod);

	unsigned int element_;
	static inline std::vector<unsigned int> primes_;
	static inline unsigned int productOfAllCharacteristics_;
	static inline std::vector<unsigned int> partials_;
	static inline constexpr unsigned int multiplicativeID_ = 1;
};

inline Shared_multi_field_element_with_small_characteristics::Shared_multi_field_element_with_small_characteristics()
	: element_(0)
{}

inline Shared_multi_field_element_with_small_characteristics::Shared_multi_field_element_with_small_characteristics(unsigned int element)
	: element_(element % productOfAllCharacteristics_)
{}

inline Shared_multi_field_element_with_small_characteristics::Shared_multi_field_element_with_small_characteristics(int element)
	: element_(element % productOfAllCharacteristics_)
{}

inline Shared_multi_field_element_with_small_characteristics::Shared_multi_field_element_with_small_characteristics(const Shared_multi_field_element_with_small_characteristics &toCopy)
	: element_(toCopy.element_)
{}

inline Shared_multi_field_element_with_small_characteristics::Shared_multi_field_element_with_small_characteristics(Shared_multi_field_element_with_small_characteristics &&toMove) noexcept
	: element_(std::exchange(toMove.element_, 0))
{}

inline void Shared_multi_field_element_with_small_characteristics::initialize(unsigned int minimum, unsigned int maximum)
{
	if (maximum < 2)
		throw std::invalid_argument("Characteristic must be strictly positive");
	if (minimum > maximum)
		throw std::invalid_argument("The given interval is not valid.");
	if (minimum == maximum && !_is_prime(minimum))
		throw std::invalid_argument("The given interval does not contain a prime number.");

	productOfAllCharacteristics_ = 1;
	primes_.clear();
	for (unsigned int i = minimum; i <= maximum; ++i){
		if (_is_prime(i)){
			primes_.push_back(i);
			productOfAllCharacteristics_ *= i;
		}
	}

	if (primes_.empty())
		throw std::invalid_argument("The given interval does not contain a prime number.");

	partials_.resize(primes_.size());
	for (unsigned int i = 0; i < primes_.size(); ++i){
		unsigned int p = primes_[i];
		unsigned int base = productOfAllCharacteristics_ / p;
		unsigned int exp = p - 1;
		partials_[i] = 1;

		while (exp > 0) {
			// If exp is odd, multiply with result
			if (exp & 1)
				partials_[i] = _multiply(partials_[i], base);
			// y must be even now
			exp = exp >> 1; // y = y/2
			base = _multiply(base, base);
		}
	}

	//If I understood the paper well, multiplicativeID_ always equals to 1. But in Clement's code,
	//multiplicativeID_ is computed (see commented loop below). TODO: verify with Clement.
//	for (unsigned int i = 0; i < partials_.size(); ++i){
//		multiplicativeID_ = (multiplicativeID_ + partials_[i]) % productOfAllCharacteristics_;
//	}
}

//inline Shared_multi_field_element_with_small_characteristics &Shared_multi_field_element_with_small_characteristics::operator+=(Shared_multi_field_element_with_small_characteristics const &f)
//{
//	_add(f.element_);
//	return *this;
//}

//inline Shared_multi_field_element_with_small_characteristics &Shared_multi_field_element_with_small_characteristics::operator+=(unsigned int const v)
//{
//	_add(v % productOfAllCharacteristics_);
//	return *this;
//}

//inline Shared_multi_field_element_with_small_characteristics &Shared_multi_field_element_with_small_characteristics::operator-=(Shared_multi_field_element_with_small_characteristics const &f)
//{
//	_substract(f.element_);
//	return *this;
//}

//inline Shared_multi_field_element_with_small_characteristics &Shared_multi_field_element_with_small_characteristics::operator-=(unsigned int const v)
//{
//	_substract(v % productOfAllCharacteristics_);
//	return *this;
//}

//inline Shared_multi_field_element_with_small_characteristics &Shared_multi_field_element_with_small_characteristics::operator*=(Shared_multi_field_element_with_small_characteristics const &f)
//{
//	element_ = _multiply(element_, f.element_);
//	return *this;
//}

//inline Shared_multi_field_element_with_small_characteristics &Shared_multi_field_element_with_small_characteristics::operator*=(unsigned int const v)
//{
//	element_ = _multiply(element_, v % productOfAllCharacteristics_);
//	return *this;
//}

inline Shared_multi_field_element_with_small_characteristics &Shared_multi_field_element_with_small_characteristics::operator=(Shared_multi_field_element_with_small_characteristics other)
{
	std::swap(element_, other.element_);
	return *this;
}

inline Shared_multi_field_element_with_small_characteristics &Shared_multi_field_element_with_small_characteristics::operator=(unsigned int const value)
{
	element_ = value % productOfAllCharacteristics_;
	return *this;
}

inline Shared_multi_field_element_with_small_characteristics::operator unsigned int() const
{
	return element_;
}

inline Shared_multi_field_element_with_small_characteristics Shared_multi_field_element_with_small_characteristics::get_inverse() const
{
	return get_partial_inverse(productOfAllCharacteristics_).first;
}

inline std::pair<Shared_multi_field_element_with_small_characteristics, unsigned int> Shared_multi_field_element_with_small_characteristics::get_partial_inverse(unsigned int productOfCharacteristics) const
{
	unsigned int gcd = std::gcd(element_, productOfAllCharacteristics_);

	if (gcd == productOfCharacteristics)
		return {Shared_multi_field_element_with_small_characteristics(), multiplicativeID_};  // partial inverse is 0

	unsigned int QT = productOfCharacteristics / gcd;
	Shared_multi_field_element_with_small_characteristics res(QT);

	const unsigned int inv_qt = _get_inverse(element_, QT);

	res = res.get_partial_multiplicative_identity();
	res *= inv_qt;

	return {res, QT};
}

inline Shared_multi_field_element_with_small_characteristics Shared_multi_field_element_with_small_characteristics::get_additive_identity()
{
	return Shared_multi_field_element_with_small_characteristics();
}

inline Shared_multi_field_element_with_small_characteristics Shared_multi_field_element_with_small_characteristics::get_multiplicative_identity()
{
	return Shared_multi_field_element_with_small_characteristics(multiplicativeID_);
}

inline Shared_multi_field_element_with_small_characteristics Shared_multi_field_element_with_small_characteristics::get_partial_multiplicative_identity()
{
	if (element_ == 0) {
		return Shared_multi_field_element_with_small_characteristics(multiplicativeID_);
	}
	Shared_multi_field_element_with_small_characteristics mult;
	for (unsigned int idx = 0; idx < primes_.size(); ++idx) {
		if ((element_ % primes_[idx]) == 0) {
			mult += partials_[idx];
		}
	}
	return mult;
}

inline unsigned int Shared_multi_field_element_with_small_characteristics::get_characteristic()
{
	return productOfAllCharacteristics_;
}

inline unsigned int Shared_multi_field_element_with_small_characteristics::get_value() const
{
	return element_;
}

inline unsigned int Shared_multi_field_element_with_small_characteristics::_add(unsigned int element, unsigned int v)
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

inline unsigned int Shared_multi_field_element_with_small_characteristics::_substract(unsigned int element, unsigned int v)
{
	if (element < v){
		element += productOfAllCharacteristics_;
	}
	element -= v;

	return element;
}

inline unsigned int Shared_multi_field_element_with_small_characteristics::_multiply(unsigned int a, unsigned int b)
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

inline constexpr int Shared_multi_field_element_with_small_characteristics::_get_inverse(unsigned int element, const unsigned int mod)
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

inline constexpr bool Shared_multi_field_element_with_small_characteristics::_is_prime(const int p)
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

#endif  // MATRIX_FIELD_MULTI_SMALL_SHARED_H_
