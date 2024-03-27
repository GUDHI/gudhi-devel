/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef MATRIX_FIELD_ZP_OPERATOR_H_
#define MATRIX_FIELD_ZP_OPERATOR_H_

#include <stdexcept>
#include <utility>
#include <vector>
#include <limits.h>

namespace Gudhi {
namespace persistence_matrix {

template<typename Unsigned_integer_type = unsigned int, class = std::enable_if_t<std::is_unsigned_v<Unsigned_integer_type> > >
class Zp_field_operators {
public:
	using element_type = Unsigned_integer_type;
	using characteristic_type = element_type;
	template <class T>
	using isSignedInteger = std::enable_if_t<std::is_signed_v<T> >;

	Zp_field_operators(characteristic_type characteristic = 0): characteristic_(0) {
		if (characteristic != 0) set_characteristic(characteristic);
	}
	Zp_field_operators(const Zp_field_operators& toCopy)
		: characteristic_(toCopy.characteristic_), inverse_(toCopy.inverse_)
	{}
	Zp_field_operators(Zp_field_operators&& toMove) noexcept
		: characteristic_(std::exchange(toMove.characteristic_, 0)), inverse_(std::move(toMove.inverse_))
	{}

	void set_characteristic(characteristic_type characteristic){
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
	}
	characteristic_type get_characteristic() const{
		return characteristic_;
	}

	element_type get_value(element_type e) const{
		return e < characteristic_ ? e : e % characteristic_;
	}

	template<typename Signed_integer_type, class = isSignedInteger<Signed_integer_type> >
	element_type get_value(Signed_integer_type e) const{
		if (e < -characteristic_) e = e % characteristic_;
		if (e < 0) return e += characteristic_;
		return e < characteristic_ ? e : e % characteristic_;
	}

	//r = e1 + e2
	element_type add(element_type e1, element_type e2) const{
		return _add(get_value(e1), get_value(e2), characteristic_);
	}

	//r = e1 - e2
	element_type substract(element_type e1, element_type e2) const{
		return _substract(get_value(e1), get_value(e2), characteristic_);
	}

	//r = e1 * e2
	element_type multiply(element_type e1, element_type e2) const{
		return _multiply(get_value(e1), get_value(e2), characteristic_);
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

	element_type get_inverse(element_type e) const{
		return inverse_[get_value(e)];
	}
	std::pair<element_type,characteristic_type> get_partial_inverse(element_type e, characteristic_type productOfCharacteristics) const{
		return {get_inverse(e), productOfCharacteristics};
	}

	static constexpr element_type get_additive_identity(){ return 0; }
	static constexpr element_type get_multiplicative_identity(){ return 1; }
	static constexpr element_type get_partial_multiplicative_identity([[maybe_unused]] characteristic_type productOfCharacteristics){ return 1; }

	static constexpr bool handles_only_z2(){
		return false;
	}

	Zp_field_operators& operator=(Zp_field_operators other){
		std::swap(characteristic_, other.characteristic_);
		inverse_.swap(other.inverse_);
		return *this;
	}
	friend void swap(Zp_field_operators& f1, Zp_field_operators& f2){
		std::swap(f1.characteristic_, f2.characteristic_);
		f1.inverse_.swap(f2.inverse_);
	}

private:
	characteristic_type characteristic_;
	std::vector<element_type> inverse_;

	static element_type _add(element_type e1, element_type e2, characteristic_type characteristic){
		if (UINT_MAX - e1 < e2) {
			//automatic unsigned integer overflow behaviour will make it work
			e1 += e2;
			e1 -= characteristic;
			return e1;
		}

		e1 += e2;
		if (e1 >= characteristic) e1 -= characteristic;

		return e1;
	}
	static element_type _substract(element_type e1, element_type e2, characteristic_type characteristic){
		if (e1 < e2){
			e1 += characteristic;
		}
		e1 -= e2;

		return e1;
	}
	static element_type _multiply(element_type e1, element_type e2, characteristic_type characteristic){
		unsigned int a = e1;
		e1 = 0;
		unsigned int temp_b;

		while (a != 0) {
			if (a & 1) {
				if (e2 >= characteristic - e1)
					e1 -= characteristic;
				e1 += e2;
			}
			a >>= 1;

			temp_b = e2;
			if (e2 >= characteristic - e2)
				temp_b -= characteristic;
			e2 += temp_b;
		}

		return e1;
	}
};

} //namespace persistence_matrix
} //namespace Gudhi

#endif  // MATRIX_FIELD_ZP_OPERATOR_H_
