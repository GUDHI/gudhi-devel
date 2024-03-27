/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef MATRIX_FIELD_MULTI_SMALL_OPERATORS_H_
#define MATRIX_FIELD_MULTI_SMALL_OPERATORS_H_

#include <utility>
#include <vector>
#include <limits.h>
#include <gmpxx.h>
#include <numeric>

namespace Gudhi {
namespace persistence_matrix {

//productOfAllCharacteristics_ ^ 2 has to fit in an unsigned int
class Multi_field_operators_with_small_characteristics {
public:
	using element_type = unsigned int;
	using characteristic_type = element_type;

	Multi_field_operators_with_small_characteristics() : productOfAllCharacteristics_(0)/* , multiplicativeID_(1) */ {}
	Multi_field_operators_with_small_characteristics(int minCharacteristic, int maxCharacteristic) 
		: productOfAllCharacteristics_(0)//, multiplicativeID_(1) 
	{
		set_characteristic(minCharacteristic, maxCharacteristic);
	}
	Multi_field_operators_with_small_characteristics(const Multi_field_operators_with_small_characteristics& toCopy) 
		: primes_(toCopy.primes_), 
		  productOfAllCharacteristics_(toCopy.productOfAllCharacteristics_),
		  partials_(toCopy.partials_)/* ,
		  multiplicativeID_(toCopy.multiplicativeID_) */
	{}
	Multi_field_operators_with_small_characteristics(Multi_field_operators_with_small_characteristics&& toMove) noexcept
		: primes_(std::move(toMove.primes_)),
		  productOfAllCharacteristics_(std::move(toMove.productOfAllCharacteristics_)),
		  partials_(std::move(toMove.partials_))/* ,
		  multiplicativeID_(std::move(toMove.multiplicativeID_)) */
	{}

	void set_characteristic(int minimum, int maximum){
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
			characteristic_type base = productOfAllCharacteristics_ / p;
			unsigned int exp = p - 1;
			partials_[i] = 1;

			while (exp > 0) {
				// If exp is odd, multiply with result
				if (exp & 1)
					partials_[i] = _multiply(partials_[i], base, productOfAllCharacteristics_);
				// y must be even now
				exp = exp >> 1; // y = y/2
				base = _multiply(base, base, productOfAllCharacteristics_);
			}
		}

		//If I understood the paper well, multiplicativeID_ always equals to 1. But in Clement's code,
		//multiplicativeID_ is computed (see commented loop below). TODO: verify with Clement.
	//	for (unsigned int i = 0; i < partials_.size(); ++i){
	//		multiplicativeID_ = (multiplicativeID_ + partials_[i]) % productOfAllCharacteristics_;
	//	}
	}
	characteristic_type get_characteristic() const{
		return productOfAllCharacteristics_;
	}

	element_type get_value(element_type e) const{
		return e < productOfAllCharacteristics_ ? e : e % productOfAllCharacteristics_;
	}

	//r = e1 + e2
	element_type add(element_type e1, element_type e2) const{
		return _add(get_value(e1), get_value(e2), productOfAllCharacteristics_);
	}

	//r = e1 - e2
	element_type substract(element_type e1, element_type e2) const{
		return _substract(get_value(e1), get_value(e2), productOfAllCharacteristics_);
	}

	//r = e1 * e2
	element_type multiply(element_type e1, element_type e2) const{
		return _multiply(get_value(e1), get_value(e2), productOfAllCharacteristics_);
	}

	//r = e * m + a
	//WARNING: not overflow safe.
	element_type multiply_and_add(element_type e, element_type m, element_type a) const{
		return get_value(e * m + a);
	}

	//r = (e + a) * m
	//WARNING: not overflow safe.
	element_type add_and_multiply(element_type e, element_type a, element_type m) const{
		return get_value((e + a) * m);
	}

	bool are_equal(element_type e1, element_type e2) const{
		return get_value(e1) == get_value(e2);
	}

	element_type get_inverse(const element_type& e) const{
		return get_partial_inverse(e, productOfAllCharacteristics_).first;
	}
	std::pair<element_type,characteristic_type> get_partial_inverse(const element_type& e, const characteristic_type& productOfCharacteristics) const{
		characteristic_type gcd = std::gcd(e, productOfAllCharacteristics_);

		if (gcd == productOfCharacteristics)
			return {0, get_multiplicative_identity()};  // partial inverse is 0

		characteristic_type QT = productOfCharacteristics / gcd;

		const characteristic_type inv_qt = _get_inverse(e, QT);

		auto res = get_partial_multiplicative_identity(QT);
		res = _multiply(res, inv_qt, productOfAllCharacteristics_);

		return {res, QT};
	}

	static constexpr element_type get_additive_identity(){ return 0; }
	static constexpr element_type get_multiplicative_identity(){ return 1; }
	// static element_type get_multiplicative_identity(){ return multiplicativeID_; }
	element_type get_partial_multiplicative_identity(const characteristic_type& productOfCharacteristics) const
	{
		if (productOfCharacteristics == 0) {
			return get_multiplicative_identity();
		}
		element_type multIdentity = 0;
		for (unsigned int idx = 0; idx < primes_.size(); ++idx) {
			if ((productOfCharacteristics % primes_[idx]) == 0) {
				multIdentity = _add(multIdentity, partials_[idx], productOfAllCharacteristics_);
			}
		}
		return multIdentity;
	}

	static constexpr bool handles_only_z2(){
		return false;
	}

	Multi_field_operators_with_small_characteristics& operator=(Multi_field_operators_with_small_characteristics other){
		primes_.swap(other.primes_);
		productOfAllCharacteristics_ = other.productOfAllCharacteristics_;
		partials_.swap(other.partials_);

		return *this;
	}

	friend void swap(Multi_field_operators_with_small_characteristics& f1, Multi_field_operators_with_small_characteristics& f2){
		f1.primes_.swap(f2.primes_);
		std::swap(f1.productOfAllCharacteristics_, f2.productOfAllCharacteristics_);
		f1.partials_.swap(f2.partials_);
	}

private:
	std::vector<unsigned int> primes_;
	characteristic_type productOfAllCharacteristics_;
	std::vector<characteristic_type> partials_;
	// static inline constexpr unsigned int multiplicativeID_ = 1;

	static element_type _add(element_type element, element_type v, characteristic_type characteristic);
	static element_type _substract(element_type element, element_type v, characteristic_type characteristic);
	static element_type _multiply(element_type a, element_type b, characteristic_type characteristic);
	static constexpr long int _get_inverse(element_type element, characteristic_type mod);
	static constexpr bool _is_prime(const int p);
};

inline Multi_field_operators_with_small_characteristics::element_type Multi_field_operators_with_small_characteristics::_add(element_type element, element_type v, characteristic_type characteristic)
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

inline Multi_field_operators_with_small_characteristics::element_type Multi_field_operators_with_small_characteristics::_substract(element_type element, element_type v, characteristic_type characteristic)
{
	if (element < v){
		element += characteristic;
	}
	element -= v;

	return element;
}

inline Multi_field_operators_with_small_characteristics::element_type Multi_field_operators_with_small_characteristics::_multiply(element_type a, element_type b, characteristic_type characteristic)
{
	element_type res = 0;
	element_type temp_b = 0;

	if (b < a) std::swap(a, b);

	while (a != 0) {
		if (a & 1) {
			/* Add b to res, modulo m, without overflow */
			if (b >= characteristic - res)
				res -= characteristic;
			res += b;
		}
		a >>= 1;

		/* Double b, modulo m */
		temp_b = b;
		if (b >= characteristic - b)
			temp_b -= characteristic;
		b += temp_b;
	}
	return res;
}

inline constexpr long int Multi_field_operators_with_small_characteristics::_get_inverse(element_type element, characteristic_type mod)
{
	//to solve: Ax + My = 1
	element_type M = mod;
	element_type A = element;
	long int y = 0, x = 1;
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

inline constexpr bool Multi_field_operators_with_small_characteristics::_is_prime(const int p)
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

#endif  // MATRIX_FIELD_MULTI_SMALL_OPERATORS_H_
