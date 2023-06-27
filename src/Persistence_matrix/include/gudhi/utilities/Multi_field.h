/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef MATRIX_FIELD_MULTI_H_
#define MATRIX_FIELD_MULTI_H_

#include <utility>
#include <vector>
#include <limits.h>
#include <iostream>
#include <gmpxx.h>
#include <stdexcept>

namespace Gudhi {
namespace persistence_matrix {

template<unsigned int minimum, unsigned int maximum>
class Multi_field_element {
public:
	using element_type = mpz_class;

	Multi_field_element();
	Multi_field_element(mpz_class element);
	Multi_field_element(const Multi_field_element& toCopy);
	Multi_field_element(Multi_field_element&& toMove) noexcept;

	friend void operator+=(Multi_field_element& f1, Multi_field_element const& f2){
		mpz_add(f1.element_.get_mpz_t(), f1.element_.get_mpz_t(), f2.element_.get_mpz_t());
		mpz_mod(f1.element_.get_mpz_t(), f1.element_.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
	}
	friend Multi_field_element operator+(Multi_field_element f1, Multi_field_element const& f2){
		f1 += f2;
		return f1;
	}
	friend void operator+=(Multi_field_element& f, mpz_class const v){
		mpz_add(f.element_.get_mpz_t(), f.element_.get_mpz_t(), v.get_mpz_t());
		mpz_mod(f.element_.get_mpz_t(), f.element_.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
	}
	friend Multi_field_element operator+(Multi_field_element f, mpz_class const v){
		f += v;
		return f;
	}
	friend mpz_class operator+(mpz_class v, Multi_field_element const& f){
		mpz_class e(v);
		mpz_add(e.get_mpz_t(), e.get_mpz_t(), f.element_.get_mpz_t());
		mpz_mod(e.get_mpz_t(), e.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
		return e;
	}

	friend void operator-=(Multi_field_element& f1, Multi_field_element const& f2){
		mpz_sub(f1.element_.get_mpz_t(), f1.element_.get_mpz_t(), f2.element_.get_mpz_t());
		mpz_mod(f1.element_.get_mpz_t(), f1.element_.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
	}
	friend Multi_field_element operator-(Multi_field_element f1, Multi_field_element const& f2){
		f1 -= f2;
		return f1;
	}
	friend void operator-=(Multi_field_element& f, mpz_class const v){
		mpz_sub(f.element_.get_mpz_t(), f.element_.get_mpz_t(), v.get_mpz_t());
		mpz_mod(f.element_.get_mpz_t(), f.element_.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
	}
	friend Multi_field_element operator-(Multi_field_element f, mpz_class const v){
		f -= v;
		return f;
	}
	friend mpz_class operator-(mpz_class v, Multi_field_element const& f){
		mpz_class e(v);
		if (e >= productOfAllCharacteristics_)
			mpz_mod(e.get_mpz_t(), e.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
		if (f.element_ > e)
			mpz_add(e.get_mpz_t(), e.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
		mpz_sub(e.get_mpz_t(), e.get_mpz_t(), f.element_.get_mpz_t());
		return e;
	}

	friend void operator*=(Multi_field_element& f1, Multi_field_element const& f2){
		mpz_mul(f1.element_.get_mpz_t(), f1.element_.get_mpz_t(), f2.element_.get_mpz_t());
		mpz_mod(f1.element_.get_mpz_t(), f1.element_.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
	}
	friend Multi_field_element operator*(Multi_field_element f1, Multi_field_element const& f2){
		f1 *= f2;
		return f1;
	}
	friend void operator*=(Multi_field_element& f, mpz_class const v){
		mpz_mul(f.element_.get_mpz_t(), f.element_.get_mpz_t(), v.get_mpz_t());
		mpz_mod(f.element_.get_mpz_t(), f.element_.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
	}
	friend Multi_field_element operator*(Multi_field_element f, mpz_class const v){
		f *= v;
		return f;
	}
	friend mpz_class operator*(mpz_class v, Multi_field_element const& f){
		mpz_class e(v);
		mpz_mul(e.get_mpz_t(), e.get_mpz_t(), f.element_.get_mpz_t());
		mpz_mod(e.get_mpz_t(), e.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
		return e;
	}

	friend bool operator==(const Multi_field_element& f1, const Multi_field_element& f2){
		return f1.element_ == f2.element_;
	}
	friend bool operator==(const mpz_class v, const Multi_field_element& f){
		if (v < productOfAllCharacteristics_) return v == f.element_;
		mpz_class e(v);
		mpz_mod(e.get_mpz_t(), e.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
		return e == f.element_;
	}
	friend bool operator==(const Multi_field_element& f, const mpz_class v){
		if (v < productOfAllCharacteristics_) return v == f.element_;
		mpz_class e(v);
		mpz_mod(e.get_mpz_t(), e.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
		return e == f.element_;
	}
	friend bool operator==(const unsigned int v, const Multi_field_element& f){
		mpz_class e(v);
		if (e < productOfAllCharacteristics_) return e == f.element_;
		mpz_mod(e.get_mpz_t(), e.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
		return e == f.element_;
	}
	friend bool operator==(const Multi_field_element& f, const unsigned int v){
		mpz_class e(v);
		if (e < productOfAllCharacteristics_) return e == f.element_;
		mpz_mod(e.get_mpz_t(), e.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
		return e == f.element_;
	}
	friend bool operator!=(const Multi_field_element& f1, const Multi_field_element& f2){
		return !(f1 == f2);
	}
	friend bool operator!=(const mpz_class v, const Multi_field_element& f){
		return !(v == f);
	}
	friend bool operator!=(const Multi_field_element& f, const mpz_class v){
		return !(v == f);
	}
	friend bool operator!=(const unsigned int v, const Multi_field_element& f){
		return !(v == f);
	}
	friend bool operator!=(const Multi_field_element& f, const unsigned int v){
		return !(v == f);
	}

	Multi_field_element& operator=(Multi_field_element other);
	Multi_field_element& operator=(const mpz_class value);
	operator unsigned int() const;
	operator mpz_class() const;
	friend void swap(Multi_field_element& f1, Multi_field_element& f2){
		std::swap(f1.element_, f2.element_);
	}

	Multi_field_element get_inverse() const;
	std::pair<Multi_field_element, mpz_class> get_partial_inverse(const mpz_class& productOfCharacteristics) const;

	static Multi_field_element get_additive_identity();
	static Multi_field_element get_multiplicative_identity();
	Multi_field_element get_partial_multiplicative_identity();
	static mpz_class get_characteristic();

	mpz_class get_value() const;

	static constexpr bool handles_only_z2(){
		return false;
	}

private:
	mpz_class element_;
	static inline const std::vector<unsigned int> primes_ = [](){
		std::vector<unsigned int> res;

		unsigned int curr_prime = minimum;
		mpz_t tmp_prime;
		mpz_init_set_ui(tmp_prime, minimum);
		// test if min_prime is prime
		int is_prime = mpz_probab_prime_p(tmp_prime, 25);  // probabilistic primality test

		if (is_prime == 0) {  // min_prime is composite
		  mpz_nextprime(tmp_prime, tmp_prime);
		  curr_prime = mpz_get_ui(tmp_prime);
		}

		while (curr_prime <= maximum) {
		  res.push_back(curr_prime);
		  mpz_nextprime(tmp_prime, tmp_prime);
		  curr_prime = mpz_get_ui(tmp_prime);
		}
		mpz_clear(tmp_prime);

		return res;
	}();
	static inline const mpz_class productOfAllCharacteristics_ = [](){
		mpz_class res = 1;
		for (const auto p : primes_) {
			res *= p;
		}

		return res;
	}();
	static inline const std::vector<mpz_class> partials_ = [](){
		std::vector<mpz_class> res;

		if (productOfAllCharacteristics_ == 1) return res;

		for (unsigned int i = 0; i < primes_.size(); ++i){
			unsigned int p = primes_[i];
			res.push_back(productOfAllCharacteristics_ / p);
			mpz_powm_ui(res.back().get_mpz_t(), res.back().get_mpz_t(), p - 1,
						productOfAllCharacteristics_.get_mpz_t());
		}

		return res;
	}();
	//If I understood the paper well, multiplicativeID_ always equals to 1. But in Clement's code,
	//multiplicativeID_ is computed (see commented lambda function below). TODO: verify with Clement.
	static inline const mpz_class multiplicativeID_ = 1;/*[](){
		mpz_class res = 0;
		for (unsigned int i = 0; i < partials_.size(); ++i){
			res = (res + partials_[i]) % productOfAllCharacteristics_;
		}

		return res;
	}();*/

	static constexpr bool _is_prime(const int p);
};

template<unsigned int minimum, unsigned int maximum>
inline Multi_field_element<minimum,maximum>::Multi_field_element()
	: element_(0)
{
	static_assert(maximum >= 2, "Characteristics have to be positive.");
	static_assert(minimum <= maximum, "The given interval is not valid.");
	static_assert(minimum != maximum || _is_prime(minimum), "The given interval does not contain a prime number.");

	if (productOfAllCharacteristics_ == 1)
		throw std::runtime_error("The given interval does not contain a prime number.");
}

template<unsigned int minimum, unsigned int maximum>
inline Multi_field_element<minimum,maximum>::Multi_field_element(mpz_class element)
	: element_(element)
{
	static_assert(maximum >= 2, "Characteristics has to be positive.");
	static_assert(minimum <= maximum, "The given interval is not valid.");
	static_assert(minimum != maximum || _is_prime(minimum), "The given interval does not contain a prime number.");

	if (productOfAllCharacteristics_ == 1)
		throw std::runtime_error("The given interval does not contain a prime number.");

	mpz_mod(element_.get_mpz_t(), element_.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
}

template<unsigned int minimum, unsigned int maximum>
inline Multi_field_element<minimum,maximum>::Multi_field_element(const Multi_field_element<minimum,maximum> &toCopy)
	: element_(toCopy.element_)
{}

template<unsigned int minimum, unsigned int maximum>
inline Multi_field_element<minimum,maximum>::Multi_field_element(Multi_field_element<minimum,maximum> &&toMove) noexcept
	: element_(std::move(toMove.element_))
{}

//template<unsigned int minimum, unsigned int maximum>
//inline Multi_field_element<minimum,maximum> &Multi_field_element<minimum,maximum>::operator+=(Multi_field_element<minimum,maximum> const &f)
//{
//	mpz_add(element_.get_mpz_t(), element_.get_mpz_t(), f.element_.get_mpz_t());
//	mpz_mod(element_.get_mpz_t(), element_.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
//	return *this;
//}

//template<unsigned int minimum, unsigned int maximum>
//inline Multi_field_element<minimum,maximum> &Multi_field_element<minimum,maximum>::operator+=(mpz_class const v)
//{
//	mpz_add(element_.get_mpz_t(), element_.get_mpz_t(), v.get_mpz_t());
//	mpz_mod(element_.get_mpz_t(), element_.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
//	return *this;
//}

//template<unsigned int minimum, unsigned int maximum>
//inline Multi_field_element<minimum,maximum> &Multi_field_element<minimum,maximum>::operator-=(Multi_field_element<minimum,maximum> const &f)
//{
//	mpz_sub(element_.get_mpz_t(), element_.get_mpz_t(), f.element_.get_mpz_t());
//	mpz_mod(element_.get_mpz_t(), element_.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
//	return *this;
//}

//template<unsigned int minimum, unsigned int maximum>
//inline Multi_field_element<minimum,maximum> &Multi_field_element<minimum,maximum>::operator-=(mpz_class const v)
//{
//	mpz_sub(element_.get_mpz_t(), element_.get_mpz_t(), v.get_mpz_t());
//	mpz_mod(element_.get_mpz_t(), element_.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
//	return *this;
//}

//template<unsigned int minimum, unsigned int maximum>
//inline Multi_field_element<minimum,maximum> &Multi_field_element<minimum,maximum>::operator*=(Multi_field_element<minimum,maximum> const &f)
//{
//	mpz_mul(element_.get_mpz_t(), element_.get_mpz_t(), f.element_.get_mpz_t());
//	mpz_mod(element_.get_mpz_t(), element_.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
//	return *this;
//}

//template<unsigned int minimum, unsigned int maximum>
//inline Multi_field_element<minimum,maximum> &Multi_field_element<minimum,maximum>::operator*=(mpz_class const v)
//{
//	mpz_mul(element_.get_mpz_t(), element_.get_mpz_t(), v.get_mpz_t());
//	mpz_mod(element_.get_mpz_t(), element_.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
//	return *this;
//}

template<unsigned int minimum, unsigned int maximum>
inline Multi_field_element<minimum,maximum> &Multi_field_element<minimum,maximum>::operator=(Multi_field_element other)
{
	std::swap(element_, other.element_);
	return *this;
}

template<unsigned int minimum, unsigned int maximum>
inline Multi_field_element<minimum,maximum> &Multi_field_element<minimum,maximum>::operator=(mpz_class const value)
{
	mpz_mod(element_.get_mpz_t(), value.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
	return *this;
}

template<unsigned int minimum, unsigned int maximum>
inline Multi_field_element<minimum,maximum>::operator unsigned int() const
{
	return element_.get_ui();
}

template<unsigned int minimum, unsigned int maximum>
inline Multi_field_element<minimum,maximum>::operator mpz_class() const
{
	return element_;
}

template<unsigned int minimum, unsigned int maximum>
inline Multi_field_element<minimum,maximum> Multi_field_element<minimum,maximum>::get_inverse() const
{
	return get_partial_inverse(productOfAllCharacteristics_).first;
}

template<unsigned int minimum, unsigned int maximum>
inline std::pair<Multi_field_element<minimum,maximum>,mpz_class> Multi_field_element<minimum,maximum>::get_partial_inverse(const mpz_class& productOfCharacteristics) const
{
	mpz_class QR;
	mpz_gcd(QR.get_mpz_t(), element_.get_mpz_t(), productOfCharacteristics.get_mpz_t());  // QR <- gcd(x,QS)

	if (QR == productOfCharacteristics)
		return {Multi_field_element(), multiplicativeID_};  // partial inverse is 0

	mpz_class QT = productOfCharacteristics / QR;
	Multi_field_element res(QT);

	mpz_class inv_qt;
	mpz_invert(inv_qt.get_mpz_t(), element_.get_mpz_t(), QT.get_mpz_t());

	res = res.get_partial_multiplicative_identity();
	res *= inv_qt;

	return {res, QT};
}

template<unsigned int minimum, unsigned int maximum>
inline Multi_field_element<minimum,maximum> Multi_field_element<minimum,maximum>::get_additive_identity()
{
	return Multi_field_element<minimum,maximum>();
}

template<unsigned int minimum, unsigned int maximum>
inline Multi_field_element<minimum,maximum> Multi_field_element<minimum,maximum>::get_multiplicative_identity()
{
	return Multi_field_element<minimum,maximum>(multiplicativeID_);
}

template<unsigned int minimum, unsigned int maximum>
inline Multi_field_element<minimum,maximum> Multi_field_element<minimum,maximum>::get_partial_multiplicative_identity()
{
	if (element_ == 0) {
		return Multi_field_element<minimum,maximum>(multiplicativeID_);
	}
	Multi_field_element<minimum,maximum> mult;
	for (unsigned int idx = 0; idx < primes_.size(); ++idx) {
		if ((element_ % primes_[idx]) == 0) {
			mult += partials_[idx];
		}
	}
	return mult;
}

template<unsigned int minimum, unsigned int maximum>
inline mpz_class Multi_field_element<minimum, maximum>::get_characteristic()
{
	return productOfAllCharacteristics_;
}

template<unsigned int minimum, unsigned int maximum>
inline mpz_class Multi_field_element<minimum,maximum>::get_value() const
{
	return element_;
}

template<unsigned int minimum, unsigned int maximum>
inline constexpr bool Multi_field_element<minimum,maximum>::_is_prime(const int p)
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

#endif  // MATRIX_FIELD_MULTI_H_
