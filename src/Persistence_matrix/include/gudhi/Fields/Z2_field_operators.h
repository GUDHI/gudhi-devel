/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef MATRIX_FIELD_Z2_OPERATORS_H_
#define MATRIX_FIELD_Z2_OPERATORS_H_

#include <utility>

namespace Gudhi {
namespace persistence_matrix {

class Z2_field_operators {
public:
	using element_type = bool;
	using characteristic_type = unsigned int;
	template <class T>
	using isUnsignedInteger = std::enable_if_t<std::is_unsigned_v<T> >;
	template <class T>
	using isInteger = std::enable_if_t<std::is_integral_v<T> >;

	Z2_field_operators(){};

	static constexpr characteristic_type get_characteristic(){
		return 2;
	}

	template<typename Integer_type, class = isInteger<Integer_type> >
	static element_type get_value(Integer_type e){
		if constexpr (std::is_same_v<Integer_type, bool>){
			return e;
		} else {
			return e < 2 && e >= 0 ? e : e % 2;		//returns bool, so %2 won't be negative and is optimized to &
		}
	}

	//r = e1 + e2
	template<typename Unsigned_integer_type, class = isUnsignedInteger<Unsigned_integer_type> >
	static element_type add(Unsigned_integer_type e1, Unsigned_integer_type e2){
		if constexpr (std::is_same_v<Unsigned_integer_type, bool>){
			return e1 != e2;
		} else {
			return get_value(e1) != get_value(e2);
		}
	}

	//r = e1 - e2
	template<typename Unsigned_integer_type, class = isUnsignedInteger<Unsigned_integer_type> >
	static element_type substract(Unsigned_integer_type e1, Unsigned_integer_type e2){
		if constexpr (std::is_same_v<Unsigned_integer_type, bool>){
			return e1 != e2;
		} else {
			return get_value(e1) != get_value(e2);
		}
	}

	//r = e1 * e2
	template<typename Unsigned_integer_type, class = isUnsignedInteger<Unsigned_integer_type> >
	static element_type multiply(Unsigned_integer_type e1, Unsigned_integer_type e2){
		if constexpr (std::is_same_v<Unsigned_integer_type, bool>){
			return e1 && e2;
		} else {
			return get_value(e1) ? get_value(e2) : false;
		}
	}

	//r = e * m + a
	template<typename Unsigned_integer_type, class = isUnsignedInteger<Unsigned_integer_type> >
	static element_type multiply_and_add(Unsigned_integer_type e, Unsigned_integer_type m, Unsigned_integer_type a){
		if constexpr (std::is_same_v<Unsigned_integer_type, bool>){
			return (e && m) != a;
		} else {
			return multiply(e, m) != get_value(a);
		}
	}

	//r = (e + a) * m
	template<typename Unsigned_integer_type, class = isUnsignedInteger<Unsigned_integer_type> >
	static element_type add_and_multiply(Unsigned_integer_type e, Unsigned_integer_type a, Unsigned_integer_type m){
		if constexpr (std::is_same_v<Unsigned_integer_type, bool>){
			return (e != a) && m;
		} else {
			return add(e, a) ? get_value(m) : false;
		}
	}

	template<typename Unsigned_integer_type, class = isUnsignedInteger<Unsigned_integer_type> >
	static bool are_equal(Unsigned_integer_type e1, Unsigned_integer_type e2){
		if constexpr (std::is_same_v<Unsigned_integer_type, bool>){
			return e1 == e2;
		} else {
			return get_value(e1) == get_value(e2);
		}
	}

	template<typename Unsigned_integer_type, class = isUnsignedInteger<Unsigned_integer_type> >
	static element_type get_inverse(Unsigned_integer_type e){
		if constexpr (std::is_same_v<Unsigned_integer_type, bool>){
			return e;
		} else {
			return get_value(e);
		}
	}
	template<typename Unsigned_integer_type, class = isUnsignedInteger<Unsigned_integer_type> >
	static std::pair<element_type,characteristic_type> get_partial_inverse(Unsigned_integer_type e, characteristic_type productOfCharacteristics){
		return {get_inverse(e), productOfCharacteristics};
	}

	static constexpr element_type get_additive_identity(){ return false; }
	static constexpr element_type get_multiplicative_identity(){ return true; }
	static constexpr element_type get_partial_multiplicative_identity([[maybe_unused]] characteristic_type productOfCharacteristics){ return true; }

	static constexpr bool handles_only_z2(){
		return true;
	}
};

} //namespace persistence_matrix
} //namespace Gudhi

#endif  // MATRIX_FIELD_Z2_OPERATORS_H_
