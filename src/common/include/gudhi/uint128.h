/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Marc Glisse
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef GUDHI_UINT128_H_
#define GUDHI_UINT128_H_
#include <gudhi/Debug_utils.h>
#include <cstdint>
#include <utility>
#include <boost/container_hash/hash.hpp>

// GUDHI_FORCE_FAKE_UINT128 is only used for tests
#if !defined __SIZEOF_INT128__ || defined GUDHI_FORCE_FAKE_UINT128
namespace Gudhi::numbers {
class Fake_uint128 {
  // Debug
 #ifdef __SIZEOF_INT128__
  unsigned __int128 native() const { return ((unsigned __int128)high << 64) + low; }
  #define GUDHI_VERIF(X) GUDHI_CHECK(res.native() == (X), "")
 #else
  #define GUDHI_VERIF(X)
 #endif
  public:
  constexpr Fake_uint128(): high(0), low(0) {}
  constexpr Fake_uint128(std::uint64_t a): high(0), low(a) {}
  // Arithmetic
  // (multiplication and division are not needed for now)
  friend Fake_uint128 operator+(Fake_uint128 a, Fake_uint128 b){
    Fake_uint128 res;
    res.low = a.low + b.low;
    res.high = a.high + b.high + (res.low < a.low);
    GUDHI_VERIF (a.native() + b.native());
    return res;
  }
  friend Fake_uint128 operator-(Fake_uint128 a, Fake_uint128 b){
    Fake_uint128 res;
    res.low = a.low - b.low;
    res.high = a.high - b.high - (res.low > a.low);
    GUDHI_VERIF (a.native() - b.native());
    return res;
  }
  friend Fake_uint128 operator<<(Fake_uint128 a, uint8_t b){
    Fake_uint128 res;
    GUDHI_CHECK(b < 128, "");
    if (b >= 64) { res.low = 0; res.high = a.low << (b-64); }
    else if (b == 0) { res = a; }
    else { res.low = a.low << b; res.high = a.high << b | a.low >> (64-b); }
    GUDHI_VERIF (a.native() << b);
    return res;
  }
  friend Fake_uint128 operator>>(Fake_uint128 a, uint8_t b){
    Fake_uint128 res;
    GUDHI_CHECK(b < 128, "");
    if (b >= 64) { res.high = 0; res.low = a.high >> (b-64); }
    else if (b == 0) { res = a; }
    else { res.high = a.high >> b; res.low = a.low >> b | a.high << (64-b); }
    GUDHI_VERIF (a.native() >> b);
    return res;
  }
  friend Fake_uint128 operator&(Fake_uint128 a, Fake_uint128 b){
    Fake_uint128 res;
    res.low = a.low & b.low;
    res.high = a.high & b.high;
    GUDHI_VERIF (a.native() & b.native());
    return res;
  }
  friend Fake_uint128 operator|(Fake_uint128 a, Fake_uint128 b){
    Fake_uint128 res;
    res.low = a.low | b.low;
    res.high = a.high | b.high;
    GUDHI_VERIF (a.native() | b.native());
    return res;
  }
  friend Fake_uint128 operator~(Fake_uint128 a){
    Fake_uint128 res;
    res.low = ~a.low;
    res.high = ~a.high;
    return res;
  }
  // In-place arithmetic
  Fake_uint128& operator+=(Fake_uint128 a) { *this = *this + a; return *this; }
  Fake_uint128& operator-=(Fake_uint128 a) { *this = *this - a; return *this; }
  Fake_uint128& operator++() { if (++low == 0) ++high; return *this; }
  Fake_uint128& operator--() { if (low-- == 0) --high; return *this; }
  Fake_uint128& operator<<=(uint8_t a) { *this = *this << a; return *this; }
  Fake_uint128& operator>>=(uint8_t a) { *this = *this >> a; return *this; }
  Fake_uint128& operator&=(Fake_uint128 a) { *this = *this & a; return *this; }
  Fake_uint128& operator|=(Fake_uint128 a) { *this = *this | a; return *this; }
  // Comparisons
  friend bool operator==(Fake_uint128 a, Fake_uint128 b){
    return a.low == b.low && a.high == b.high;
  }
  friend bool operator!=(Fake_uint128 a, Fake_uint128 b){
    return a.low != b.low || a.high != b.high;
  }
  friend bool operator<(Fake_uint128 a, Fake_uint128 b){
    return a.high < b.high || (a.high == b.high && a.low < b.low);
  }
  friend bool operator>(Fake_uint128 a, Fake_uint128 b){
    return a.high > b.high || (a.high == b.high && a.low > b.low);
  }
  friend bool operator<=(Fake_uint128 a, Fake_uint128 b){
    return a.high < b.high || (a.high == b.high && a.low <= b.low);
  }
  friend bool operator>=(Fake_uint128 a, Fake_uint128 b){
    return a.high > b.high || (a.high == b.high && a.low >= b.low);
  }
  // Misc
  friend std::size_t hash_value(Fake_uint128 a) {
    typedef std::pair<std::uint64_t, std::uint64_t> P;
    return boost::hash_value(P(a.high, a.low));
  }
  template <class T, class=std::enable_if_t<std::is_integral_v<T>>>
  explicit operator T() const {
    GUDHI_CHECK(high == 0 && low <= std::numeric_limits<T>::max(), "");
    return static_cast<T>(low);
  }
  private:
  std::uint64_t high, low; // does the order matter?
 #undef GUDHI_VERIF
};
typedef Fake_uint128 uint128_t;
} // namespace Gudhi::numbers
template<> class std::numeric_limits<Gudhi::numbers::Fake_uint128> {
  public:
  static constexpr bool is_specialized = true;
  static constexpr bool is_signed = false;
  static constexpr bool is_integer = true;
  static constexpr bool is_exact = true;
  static constexpr bool has_infinity = false;
  static constexpr bool is_modulo = true;
  static constexpr int digits = 128;
  static constexpr int radix = 2;
  // etc
};
#else
namespace Gudhi::numbers {
typedef unsigned __int128 uint128_t;
} // namespace Gudhi::numbers
#endif
#endif // GUDHI_UINT128_H_
