/* Based on Ripser commit 140670f2c76997404601e43d8054151f46be9fd7
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 *      - 2024 Marc Glisse: Heavy refactoring
*/

/*

 Ripser: a lean C++ code for computation of Vietoris-Rips persistence barcodes

 MIT License

 Copyright (c) 2015–2021 Ulrich Bauer

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.

 You are under no obligation whatsoever to provide any bug fixes, patches, or
 upgrades to the features, functionality or performance of the source code
 ("Enhancements") to anyone; however, if you choose to make your Enhancements
 available either publicly, or directly to the author of this software, without
 imposing a separate written license agreement for such Enhancements, then you
 hereby grant the following license: a non-exclusive, royalty-free perpetual
 license to install, use, modify, prepare derivative works, incorporate into
 other computer software, distribute, and sublicense such enhancements or
 derivative works thereof, in binary and source code form.

*/

//#define GUDHI_INDICATE_PROGRESS

// #define GUDHI_RIPSER_USE_BOOST_HEAP // only useful when heap::push dominates?
// #define GUDHI_RIPSER_USE_HASHMAP_FOR_SPARSE_DIST_MAT // only useful in cases where edge-collapse is more important?

/* TODO:
 * from branch representative-cocycles
   - for vertices: brute force check all vertices that are in the same connected component (could do better, but there is already quadraticness elsewhere in most cases, even the output for dim 0 may be quadratic)
   - for others: print_chain just calls iteratively get_pivot on working_reduction_column
 * from branch representative-cycles
   - dim 0: trivial
   - parametrize a number of functions (like add_coboundary) by the (co)boundary iterator so they can be also used for homology
   - once cohomology is computed for some dim, assemble the relevant simplices (don't forget the essential ones) and reduce them homology-style. I think we have 2 choices: reduce the birth and check which columns we add (only possibility for essential classes, but do we really need to do cohomology first if we are going to do that?), or reduce the death and look at the column after reduction (Ripser(er) chose this).
 * check out the persistence image branch

 * allow non-0 filtration value on vertices, so we can handle all flag-type filtrations, not just plain Rips. cf scikit-tda
*/

#ifndef GUDHI_RIPSER_H
#define GUDHI_RIPSER_H

#include <algorithm>
#include <cmath> // sqrt
#include <numeric>
#include <queue>
#include <optional>
#include <stdexcept>
#ifdef GUDHI_INDICATE_PROGRESS
#include <chrono>
#endif

#include <boost/range/iterator_range_core.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION >= 108100
#include <boost/unordered/unordered_flat_map.hpp>
#else
#include <boost/unordered_map.hpp>
#endif

#ifdef GUDHI_RIPSER_USE_BOOST_HEAP
#include <boost/heap/d_ary_heap.hpp>
#endif

#include <gudhi/uint128.h>
#include <gudhi/Debug_utils.h>


namespace Gudhi::ripser {

#define GUDHI_assert(X) GUDHI_CHECK(X, std::logic_error(""))

#ifdef GUDHI_RIPSER_USE_BOOST_HEAP
template <class T, class, class C>
using Heap = boost::heap::d_ary_heap<T, boost::heap::arity<8>, boost::heap::compare<C>>;
#else
template <class T, class V, class C>
struct Heap : std::priority_queue<T, V, C> {
  typedef std::priority_queue<T, V, C> Base;
  using Base::Base;
  void clear() { this->c.clear(); }
};
#endif

template<class vertex_t_>
class Union_find {
  public:
    typedef vertex_t_ vertex_t;
  private:
    std::vector<vertex_t> parent;
    std::vector<uint8_t> rank;
  public:
    Union_find(const vertex_t n) : parent(n), rank(n, 0) {
      for (vertex_t i = 0; i < n; ++i) parent[i] = i;
    }

    // TODO: check which one is faster (probably doesn't matter)
    vertex_t find(vertex_t x) {
#if 0
      vertex_t y = x, z;
      while ((z = parent[y]) != y) y = z;
      while ((z = parent[x]) != y) {
        parent[x] = y;
        x = z;
      }
      return z;
#else
         // Path halving
      vertex_t y = parent[x];
      vertex_t z = parent[y];
      while (y != z)
      {
        parent[x] = z;
        x = z;
        y = parent[x];
        z = parent[y];
      }
      return y;
#endif
    }

    void link(vertex_t x, vertex_t y) {
      // this line is redundant, the caller already does it
      if ((x = find(x)) == (y = find(y))) return;
      if (rank[x] > rank[y])
        parent[y] = x;
      else {
        parent[x] = y;
        if (rank[x] == rank[y]) ++rank[y];
      }
    }
};

template<class coefficient_t>
bool is_prime(const coefficient_t n) {
  if (!(n & 1) || n < 2) return n == 2;
  for (coefficient_t p = 3; p * p <= n; p += 2)
    if (!(n % p)) return false;
  return true;
}

// For pretty printing, modulo 11, we prefer -1 to 10.
template<class coefficient_t>
coefficient_t normalize(const coefficient_t n, const coefficient_t modulus) {
  return n > modulus / 2 ? n - modulus : n;
}

template<class coefficient_storage_t, class coefficient_t>
std::vector<coefficient_storage_t> multiplicative_inverse_vector(const coefficient_t m) {
  std::vector<coefficient_storage_t> inverse(m);
  if (!is_prime(m))
    throw std::domain_error("Modulus must be a prime number");
  if ((m - 1) != (coefficient_storage_t)(m - 1))
    throw std::overflow_error("Modulus is too large");
  inverse[1] = 1;
  // m = a * (m / a) + m % a
  // Multipying with inverse(a) * inverse(m % a):
  // 0 = inverse(m % a) * (m / a) + inverse(a)  (mod m)
  for (coefficient_t a = 2; a < m; ++a) inverse[a] = m - (inverse[m % a] * (m / a)) % m;
  return inverse;
}

#ifdef GUDHI_INDICATE_PROGRESS
constexpr std::chrono::milliseconds time_step(40);
constexpr const char* clear_line="\r\033[K";
#endif

/* concept DistanceMatrix
struct DistanceMatrix {
  typedef Category;
  typedef vertex_t;
  typedef value_t;
  value_t operator()(vertex_t i, vertex_t j) const;
  vertex_t size() const;
  // Optional:
  static const bool is_sparse;
  vector<vector<vertex_diameter_t>> neighbors;
};
*/

class Tag_dense {}; // Use directly, iterate on vertices
class Tag_sparse {}; // Use directly, iterate on edges
class Tag_other {}; // Could use directly as dense, but prefer converting it.

template <class Params>
struct Full_distance_matrix {
  public:
    typedef Tag_dense Category;
    typedef typename Params::vertex_t vertex_t;
    typedef typename Params::value_t value_t;
    std::vector<value_t> distances; // TODO: private
  private:
    vertex_t n;
  public:
    //Full_distance_matrix(std::vector<value_t>&& _distances)
    //  : distances(std::move(_distances)), n(std::sqrt(_distances.size())) {}

    template <typename DistanceMatrix>
      Full_distance_matrix(const DistanceMatrix& mat)
      : distances(static_cast<std::size_t>(mat.size()) * mat.size()), n(mat.size()) { // vertex_t is meant for vertices. Using it for edges could be unsafe, so we cast to size_t.
        for (vertex_t i = 0; i < size(); ++i)
          for (vertex_t j = 0; j < size(); ++j)
            distances[i * n + j] = mat(i, j); // If mat::operator() involves computations, should we try to take advantage of the symmetry?
      }

    // Confusing: why is the order j,i significantly faster than i,j?
    value_t operator()(const vertex_t j, const vertex_t i) const {
      return distances[i * n + j];
    }
    vertex_t size() const { return n; }
};

enum Compressed_matrix_layout { LOWER_TRIANGULAR, UPPER_TRIANGULAR };

template <class Params, Compressed_matrix_layout Layout>
struct Compressed_distance_matrix {
  public:
    typedef Tag_dense Category;
    typedef typename Params::vertex_t vertex_t;
    typedef typename Params::value_t value_t;
    std::vector<value_t> distances; // TODO: private
  private:
    std::vector<value_t*> rows; // Surprisingly, this is more efficient than computing i*(i-1)/2

  public:
    Compressed_distance_matrix(std::vector<value_t>&& _distances)
      : distances(std::move(_distances)), rows((1 + std::sqrt(1 + 8 * distances.size())) / 2) {
        GUDHI_assert(distances.size() == (std::size_t)size() * (size() - 1) / 2);
        init_rows();
      }

    template <typename DistanceMatrix>
      Compressed_distance_matrix(const DistanceMatrix& mat)
      : distances(static_cast<std::size_t>(mat.size()) * (mat.size() - 1) / 2), rows(mat.size()) { // vertex_t is meant for vertices. Using it for edges could be unsafe, so we cast to size_t.
        init_rows();

        for (vertex_t i = 1; i < size(); ++i)
          for (vertex_t j = 0; j < i; ++j) rows[i][j] = mat(i, j);
      }

    value_t operator()(const vertex_t i, const vertex_t j) const {
      if (i == j) return 0;
      if ((Layout == LOWER_TRIANGULAR) ? (i < j) : (i > j))
        return rows[j][i];
      else
        return rows[i][j];
    }
    vertex_t size() const { return rows.size(); }
  private:
    void init_rows() {
      if constexpr (Layout == LOWER_TRIANGULAR) {
        value_t* pointer = &distances[0];
        for (vertex_t i = 1; i < size(); ++i) {
          rows[i] = pointer;
          pointer += i;
        }
      } else { // UPPER_TRIANGULAR
        value_t* pointer = &distances[0] - 1;
        for (vertex_t i = 0; i < size() - 1; ++i) {
          rows[i] = pointer;
          pointer += size() - i - 2;
        }
      }
    }
};

template <class Params>
struct Sparse_distance_matrix {
  public:
    typedef Tag_sparse Category;
    static constexpr bool is_sparse = true;
    typedef typename Params::vertex_t vertex_t;
    typedef typename Params::value_t value_t;
    struct vertex_diameter_t {
      vertex_diameter_t() =default;
      vertex_diameter_t(vertex_t i_, value_t d_) : i(i_), d(d_) {}
#if 1
      // gcc is faster with those 2 lines (or with a std::pair) :-(
      vertex_diameter_t(vertex_diameter_t const&o)noexcept :i(o.i),d(o.d){}
      vertex_diameter_t&operator=(vertex_diameter_t const&o) noexcept{i=o.i;d=o.d;return*this;}
#endif
      vertex_t i; value_t d;
      friend vertex_t get_vertex(const vertex_diameter_t& i) { return i.i; }
      friend value_t get_diameter(const vertex_diameter_t& i) { return i.d; }
      friend bool operator<(vertex_diameter_t const& a, vertex_diameter_t const& b) {
        if (a.i < b.i) return true;
        if (a.i > b.i) return false;
        return a.d < b.d;
      }
    };

    std::vector<std::vector<vertex_diameter_t>> neighbors;
    std::size_t num_edges; // TODO: useless, remove

  private:
#ifdef GUDHI_RIPSER_USE_HASHMAP_FOR_SPARSE_DIST_MAT
#if BOOST_VERSION >= 108100
    boost::unordered_flat_map<std::pair<vertex_t,vertex_t>,value_t> m;
#else
    boost::unordered_map<std::pair<vertex_t,vertex_t>,value_t> m;
#endif
    // Would a vector<unordered_map> (one map per vertex) have any advantage?
#endif

  public:
    Sparse_distance_matrix(std::vector<std::vector<vertex_diameter_t>>&& _neighbors,
        std::size_t _num_edges = 0)
      : neighbors(std::move(_neighbors)), num_edges(_num_edges) {init();}

    template <typename DistanceMatrix>
      Sparse_distance_matrix(const DistanceMatrix& mat, const value_t threshold)
      : neighbors(mat.size()), num_edges(0) {

        for (vertex_t i = 0; i < size(); ++i)
          for (vertex_t j = 0; j < size(); ++j)
            if (i != j) {
              auto d = mat(i, j);
              if (d <= threshold) {
                ++num_edges;
                neighbors[i].emplace_back(j, d);
              }
            }
        init();
      }

    value_t operator()(const vertex_t i, const vertex_t j) const {
#ifdef GUDHI_RIPSER_USE_HASHMAP_FOR_SPARSE_DIST_MAT
      // We could insert in both orders to save the minmax for each query.
      return m.at(std::minmax(i,j));
      //// We never hit the infinity case?
      // auto it = m.find(std::minmax(i,j));
      // return (it != m.end()) ? *it : std::numeric_limits<value_t>::infinity();
#else
      auto neighbor =
        std::lower_bound(neighbors[i].begin(), neighbors[i].end(), vertex_diameter_t{j, 0});
      return (neighbor != neighbors[i].end() && get_vertex(*neighbor) == j)
        ? get_diameter(*neighbor)
        : std::numeric_limits<value_t>::infinity();
#endif
    }
    vertex_t size() const { return neighbors.size(); }
  private:
    void init() {
#ifdef GUDHI_RIPSER_USE_HASHMAP_FOR_SPARSE_DIST_MAT
      for(vertex_t i=0; i<size(); ++i){
        for(auto n:neighbors[i]){
          m[std::minmax(i,get_vertex(n))]=get_diameter(n);
        }
      }
#endif
    }
};

// Do not feed this directly to ripser (slow), first convert to another matrix type
template <class Params>
struct Euclidean_distance_matrix {
  public:
    typedef Tag_other Category;
    typedef typename Params::vertex_t vertex_t;
    typedef typename Params::value_t value_t;

    Euclidean_distance_matrix(std::vector<std::vector<value_t>>&& _points)
      : points(std::move(_points)) {
        for (auto p : points) { GUDHI_assert(p.size() == points.front().size()); }
      }

    value_t operator()(const vertex_t i, const vertex_t j) const {
      GUDHI_assert((std::size_t)i < points.size());
      GUDHI_assert((std::size_t)j < points.size());
      return std::sqrt(std::inner_product(
            points[i].begin(), points[i].end(), points[j].begin(), value_t(), std::plus<value_t>(),
            [](value_t u, value_t v) { return (u - v) * (u - v); }));
    }

    vertex_t size() const { return points.size(); }
  private:
    std::vector<std::vector<value_t>> points;
};

// The gratuitous restrictions on what can be specialized in C++ are annoying.
template <class DistanceMatrix, class=std::bool_constant<true>> struct Is_sparse_impl : std::bool_constant<false> {};
template <class DistanceMatrix> struct Is_sparse_impl<DistanceMatrix, std::bool_constant<DistanceMatrix::is_sparse>> : std::bool_constant<true> {};
template <class DistanceMatrix> constexpr bool is_sparse () { return Is_sparse_impl<DistanceMatrix>::value; }

template <typename ValueType> class Compressed_sparse_matrix_ {
  std::vector<size_t> bounds;
  std::vector<ValueType> entries;

  typedef typename std::vector<ValueType>::const_iterator iterator;
  typedef boost::iterator_range<iterator> iterator_pair;

  public:
  iterator_pair subrange(const size_t index) const {
    return {entries.begin() + (index == 0 ? 0 : bounds[index - 1]),
      entries.begin() + bounds[index]};
  }

  // Close the current column, the next push_back will be for a new column
  void append_column() { bounds.push_back(entries.size()); }

  void push_back(const ValueType e) {
    GUDHI_assert(0 < bounds.size());
    entries.push_back(e);
    ++bounds.back();
  }
};

/* concept SimplexEncoding
struct SimplexEncoding {
  typedef dimension_t;
  typedef vertex_t;
  typedef simplex_t; // has to be an integer
  SimplexEncoding(vertex_t n_vertices, dimension_t max_dim_plus_one); // max_vertices_per_simplex
  simplex_t operator()(vertex_t v, dimension_t position) const; // position starts at 1?
  vertex_t get_max_vertex(simplex_t s, dimension_t position, vertex_t n_vertices) const;
  int num_extra_bits() const;
};
*/

// number of bits necessary to store x with 0 <= x < n
template <class vertex_t>
constexpr int log2up(vertex_t n) {
  --n;
  int k = 0;
  while(n>0) { n>>=1; ++k; }
  return k;
}

template <class Params>
class Cns_encoding {
  public:
    typedef typename Params::dimension_t dimension_t;
    typedef typename Params::vertex_t vertex_t;
    typedef typename Params::simplex_t simplex_t;

  private:
    std::vector<std::vector<simplex_t>> B; // table of binomial coefficients
    int extra_bits;

  public:
    Cns_encoding(vertex_t n, dimension_t k) : B(k + 1, std::vector<simplex_t>(n + 1, 0)) {
      static_assert(std::numeric_limits<simplex_t>::radix == 2);
      const int available_bits = std::numeric_limits<simplex_t>::digits;
      simplex_t max_simplex_index = 0;
      for (vertex_t i = 0; i <= n; ++i) {
        B[0][i] = 1;
        for (dimension_t j = 1; (vertex_t)j < std::min<vertex_t>(i, k + 1); ++j)
          B[j][i] = B[j - 1][i - 1] + B[j][i - 1];
        if (i <= k) B[i][i] = 1;
        vertex_t mi = std::min<vertex_t>(i >> 1, (vertex_t)k); // max
        max_simplex_index = B[mi][i];
        if (i > 1 && max_simplex_index < B[mi][i-1]) { // overflow
          throw std::overflow_error("cannot encode all simplices of dimension " + std::to_string(k) + " with " + std::to_string(n) + " vertices using only " + std::to_string(available_bits) + " bits");
        }
      }
      extra_bits = available_bits - log2up(max_simplex_index + 1);
    }

    simplex_t operator()(vertex_t n, dimension_t k) const {
      GUDHI_assert(n >= k - 1);
      return B[k][n];
    }

    // We could get `n` from B and avoid passing it as argument
    vertex_t get_max_vertex(const simplex_t idx, const dimension_t k, const vertex_t n) const {
      return get_max(n, k - 1, [&](vertex_t w) -> bool { return (*this)(w, k) <= idx; });
    }

    int num_extra_bits() const { return extra_bits; }

  private:
    template <class Predicate>
      static vertex_t get_max(vertex_t top, const vertex_t bottom, const Predicate pred) {
        if (!pred(top)) {
          vertex_t count = top - bottom;
          while (count > 0) {
            vertex_t step = count >> 1, mid = top - step;
            if (!pred(mid)) {
              top = mid - 1;
              count -= step + 1;
            } else
              count = step;
          }
        }
        return top;
      }
};

template <class Params>
class Bitfield_encoding {
  public:
    typedef typename Params::dimension_t dimension_t;
    typedef typename Params::vertex_t vertex_t;
    typedef typename Params::simplex_t simplex_t;

  private:
    int bits_per_vertex;
    int extra_bits;

  public:
    Bitfield_encoding(vertex_t n, dimension_t k) : bits_per_vertex(log2up(n)) {
      static_assert(std::numeric_limits<simplex_t>::radix == 2);
      const int available_bits = std::numeric_limits<simplex_t>::digits;
      extra_bits = available_bits - bits_per_vertex * k;
      if (extra_bits < 0)
        throw std::overflow_error("cannot encode all simplices of dimension " + std::to_string(k - 1) + " with " + std::to_string(n) + " vertices using only " + std::to_string(available_bits) + " bits");
      // The message is a bit misleading, it is tuples that we cannot encode, and just with this representation.
    }

    simplex_t operator()(vertex_t n, dimension_t k) const {
      if(k==0) return 1; // because of odd use in (co)boundary...
      --k;
      // USE_N_MINUS_K only useful if it somehow helps remove the test k==0 above
#ifdef USE_N_MINUS_K
      return (simplex_t)(n - k) << (bits_per_vertex * k);
#else
      return (simplex_t)n << (bits_per_vertex * k);
#endif
    }

    vertex_t get_max_vertex(const simplex_t idx, dimension_t k, const vertex_t) const {
      GUDHI_assert(k > 0);
      --k;
#ifdef USE_N_MINUS_K
      return static_cast<vertex_t>(idx >> (bits_per_vertex * k)) + k;
#else
      return static_cast<vertex_t>(idx >> (bits_per_vertex * k));
#endif
    }

    int num_extra_bits() const { return extra_bits; }
};

template <typename DistanceMatrix, typename SimplexEncoding, typename Params> struct Rips_filtration {
  using size_t = typename Params::size_t;
  using vertex_t = typename SimplexEncoding::vertex_t;
  // static_assert(std::is_same_v<vertex_t, typename DistanceMatrix::vertex_t>); // too strict
  using simplex_t = typename SimplexEncoding::simplex_t;
  using dimension_t = typename SimplexEncoding::dimension_t;
  using value_t = typename DistanceMatrix::value_t;
  using coefficient_storage_t = typename Params::coefficient_storage_t;
  using coefficient_t = typename Params::coefficient_t;
  static constexpr bool use_coefficients = Params::use_coefficients;

  // The definition of entry_t could be added in some intermediate layer between SimplexEncoding and here
  struct entry_with_coeff_t {
    simplex_t content;
    entry_with_coeff_t(simplex_t _index, coefficient_t _coefficient, int bits_for_coeff)
      : content((_index << bits_for_coeff) | (_coefficient - 1)) { GUDHI_assert(_coefficient != 0); }
    entry_with_coeff_t() {}
    // We never store a coefficient of 0, so we can store coef-1 so %2 requires 0 bit and %3 only 1 bit
    friend const entry_with_coeff_t& get_entry(const entry_with_coeff_t& e) { return e; }
  };
  simplex_t get_index(const entry_with_coeff_t& e) const { return e.content >> num_bits_for_coeff(); }
  coefficient_t get_coefficient(const entry_with_coeff_t& e) const { return static_cast<coefficient_t>(e.content & (((simplex_t)1 << num_bits_for_coeff()) - 1)) + 1; }
  void set_coefficient(entry_with_coeff_t& e, const coefficient_t c) const { GUDHI_assert(c!=0); e.content = (e.content & ((simplex_t)(-1) << num_bits_for_coeff())) | (c - 1); }
  // Should we cache the masks derived from num_bits_for_coeff?

  struct entry_plain_t {
    simplex_t index;
    entry_plain_t(simplex_t _index, coefficient_t, int) : index(_index) {}
    entry_plain_t() {}
    friend const entry_plain_t& get_entry(const entry_plain_t& e) { return e; }
  };
  static simplex_t get_index(const entry_plain_t& i) { return i.index; }
  static coefficient_t get_coefficient(const entry_plain_t& i) { return 1; }
  static void set_coefficient(entry_plain_t& e, const coefficient_t c) { GUDHI_assert(c==1); }

  typedef std::conditional_t<use_coefficients, entry_with_coeff_t, entry_plain_t> entry_t;
  entry_t make_entry(simplex_t i, coefficient_t c) const { return entry_t(i, c, num_bits_for_coeff()); }

  static_assert(sizeof(entry_t) == sizeof(simplex_t), "size of entry_t is not the same as simplex_t");

  // TODO: avoid storing filtp when !use_coefficients (same for Equal_index and greater_*)
  struct Entry_hash {
    Entry_hash(Rips_filtration const& filt) : filtp(&filt) {}
    Rips_filtration const* filtp;
    std::size_t operator()(const entry_t& e) const { return boost::hash<simplex_t>()(filtp->get_index(e)); }
  };

  struct Equal_index {
    Equal_index(Rips_filtration const& filt) : filtp(&filt) {}
    Rips_filtration const* filtp;
    bool operator()(const entry_t& e, const entry_t& f) const {
      return filtp->get_index(e) == filtp->get_index(f);
    }
  };

  struct diameter_simplex_t {
    value_t diameter;
    simplex_t index;
    friend value_t get_diameter(const diameter_simplex_t& i) { return i.diameter; }
  };
  static simplex_t get_index(const diameter_simplex_t& i) { return i.index; }

  struct diameter_entry_t : std::pair<value_t, entry_t> {
    using std::pair<value_t, entry_t>::pair;
    friend const entry_t& get_entry(const diameter_entry_t& p) { return p.second; }
    friend entry_t& get_entry(diameter_entry_t& p) { return p.second; }
    friend value_t get_diameter(const diameter_entry_t& p) { return p.first; }
  };
  simplex_t get_index(const diameter_entry_t& p) const { return get_index(get_entry(p)); }
  coefficient_t get_coefficient(const diameter_entry_t& p) const {
    return get_coefficient(get_entry(p));
  }
  void set_coefficient(diameter_entry_t& p, const coefficient_t c) const {
    set_coefficient(get_entry(p), c);
  }
  diameter_entry_t make_diameter_entry(value_t diameter, simplex_t index, coefficient_t coefficient) const {
    return diameter_entry_t(diameter, make_entry(index, coefficient));
  }
  diameter_entry_t make_diameter_entry(const diameter_simplex_t& diameter_index, coefficient_t coefficient) const {
    return diameter_entry_t(get_diameter(diameter_index), make_entry(get_index(diameter_index), coefficient));
  }

  template <typename Entry> struct Greater_diameter_or_smaller_index {
    Greater_diameter_or_smaller_index(Rips_filtration const& filt) : filtp(&filt) {}
    Rips_filtration const* filtp;
    bool operator()(const Entry& a, const Entry& b) const {
      return (get_diameter(a) > get_diameter(b)) ||
        ((get_diameter(a) == get_diameter(b)) && (filtp->get_index(a) < filtp->get_index(b)));
    }
  };

  const DistanceMatrix dist; // only store a reference instead?
  const vertex_t n; // redundant with dist?
  const dimension_t dim_max;
  const value_t threshold; // It would be nice if this was only in DistanceMatrix, but inconvenient. Only used to list the edges of dense distance matrices.
  const coefficient_t modulus;
  const SimplexEncoding simplex_encoding; // only store a reference instead?
  mutable std::vector<vertex_t> vertices; // we must not have several threads looking at the same complex
  int bits_for_coeff;

  Rips_filtration(DistanceMatrix&& _dist, dimension_t _dim_max, value_t _threshold, coefficient_t _modulus)
    : dist(std::move(_dist)), n(dist.size()),
    dim_max(std::min<vertex_t>(_dim_max, dist.size() - 2)), threshold(_threshold),
    modulus(_modulus), simplex_encoding(n, dim_max + 2), bits_for_coeff(log2up(modulus-1)) {
      // See entry_with_coeff_t for the logic for log2up(modulus-1) (storing coeff-1)
      if (use_coefficients && simplex_encoding.num_extra_bits() < num_bits_for_coeff())
        // TODO: include relevant numbers in the message
        throw std::overflow_error("Not enough spare bits in the simplex encoding to store a coefficient");
    }

  vertex_t num_vertices() const { return n; }
  int num_bits_for_coeff() const { return bits_for_coeff; }

  simplex_t get_edge_index(const vertex_t i, const vertex_t j) const {
    return simplex_encoding(i, 2) + j;
  }

  template <typename OutputIterator>
    OutputIterator get_simplex_vertices(simplex_t idx, const dimension_t dim, vertex_t n,
        OutputIterator out) const {
      --n;
      for (dimension_t k = dim + 1; k > 1; --k) {
        n = simplex_encoding.get_max_vertex(idx, k, n);
        *out++ = n;
        idx -= simplex_encoding(n, k);
      }
      *out = static_cast<vertex_t>(idx);
      return out;
    }

  value_t compute_diameter(const simplex_t index, const dimension_t dim) const {
    value_t diam = -std::numeric_limits<value_t>::infinity();

    vertices.resize(dim + 1);
    get_simplex_vertices(index, dim, dist.size(), vertices.rbegin());

    for (dimension_t i = 0; i <= dim; ++i)
      for (dimension_t j = 0; j < i; ++j) {
        diam = std::max(diam, dist(vertices[i], vertices[j]));
      }
    return diam;
  }

  std::vector<diameter_simplex_t> get_edges() {
    if constexpr (!std::is_same_v<typename DistanceMatrix::Category, Tag_sparse>) { // Compressed_lower_distance_matrix
      // TODO: it would be convenient to have DistanceMatrix provide a range of neighbors at dist<=threshold even in the dense case (as a filtered_range)
      std::vector<diameter_simplex_t> edges;
      for (vertex_t i = 0; i < n; ++i) {
        for (vertex_t j = 0; j < i; ++j) {
          value_t length = dist(i, j);
          if (length <= threshold) edges.push_back({length, get_edge_index(i, j)});
        }
      }
      return edges;
    } else { // Sparse_distance_matrix
      std::vector<diameter_simplex_t> edges;
      for (vertex_t i = 0; i < n; ++i)
        for (auto n : dist.neighbors[i]) {
          vertex_t j = get_vertex(n);
          if (i > j) edges.push_back({get_diameter(n), get_edge_index(i, j)});
        }
      return edges;
    }
  }

  // TODO: document in what way (if any) the order matters
  template<class DistanceMatrix2, class=typename DistanceMatrix2::Category> class Simplex_coboundary_enumerator_ { // Compressed_lower_distance_matrix
    simplex_t idx_below, idx_above;
    vertex_t j;
    dimension_t k;
    std::vector<vertex_t> vertices;
    diameter_entry_t simplex;
    const DistanceMatrix2& dist;
    const SimplexEncoding& simplex_encoding;
    const Rips_filtration& parent; // for n and get_simplex_vertices
    // at least dist and simplex_encoding are redundant with parent, but using parent.dist and parent.simplex_encoding seems to have a bad impact on performance.

    public:
    Simplex_coboundary_enumerator_(const Rips_filtration& _parent) : dist(_parent.dist),
    simplex_encoding(_parent.simplex_encoding), parent(_parent) {}

    void set_simplex(const diameter_entry_t _simplex, const dimension_t _dim) {
      idx_below = parent.get_index(_simplex);
      idx_above = 0;
      j = dist.size() - 1;
      k = _dim + 1;
      simplex = _simplex;
      vertices.resize(_dim + 1);
      parent.get_simplex_vertices(parent.get_index(_simplex), _dim, dist.size(), vertices.rbegin());
    }

    bool has_next(bool all_cofacets = true) {
      return (j >= k && (all_cofacets || simplex_encoding(j, k) > idx_below));
    }

    std::optional<diameter_entry_t> next_raw(bool all_cofacets = true) {
      if (!has_next(all_cofacets)) return std::nullopt;
      // this requires simplex_encoding(x,0)>0
      while (simplex_encoding(j, k) <= idx_below) {
        idx_below -= simplex_encoding(j, k);
        idx_above += simplex_encoding(j, k + 1);
        --j;
        --k;
        GUDHI_assert(k != -1);
      }
      value_t cofacet_diameter = get_diameter(simplex);
      // The order of j and i matters for performance
      for (vertex_t i : vertices) cofacet_diameter = std::max(cofacet_diameter, dist(j, i));
      simplex_t cofacet_index = idx_above + simplex_encoding(j--, k + 1) + idx_below;
      coefficient_t cofacet_coefficient = parent.get_coefficient(simplex);
      if (k & 1) cofacet_coefficient = parent.modulus - cofacet_coefficient;
      return parent.make_diameter_entry(cofacet_diameter, cofacet_index, cofacet_coefficient);
    }
    std::optional<diameter_entry_t> next(bool all_cofacets = true) {
      while(true) {
        std::optional<diameter_entry_t> res = next_raw(all_cofacets);
        if(!res || get_diameter(*res) <= parent.threshold) return res;
      }
    }
  };

  template <class DistanceMatrix2> class Simplex_coboundary_enumerator_<DistanceMatrix2, Tag_sparse> {
    typedef typename DistanceMatrix2::vertex_diameter_t vertex_diameter_t;
    simplex_t idx_below, idx_above;
    dimension_t k;
    std::vector<vertex_t> vertices;
    diameter_entry_t simplex;
    const DistanceMatrix2& dist;
    const SimplexEncoding& simplex_encoding;
    std::vector<typename std::vector<vertex_diameter_t>::const_reverse_iterator> neighbor_it;
    std::vector<typename std::vector<vertex_diameter_t>::const_reverse_iterator> neighbor_end;
    vertex_diameter_t neighbor;
    const Rips_filtration& parent; // for n and get_simplex_vertices

    public:
    Simplex_coboundary_enumerator_(const Rips_filtration& _parent)
      : dist(_parent.dist),
      simplex_encoding(_parent.simplex_encoding), parent(_parent) {}

    void set_simplex(const diameter_entry_t _simplex, const dimension_t _dim) {
      idx_below = parent.get_index(_simplex);
      idx_above = 0;
      k = _dim + 1;
      simplex = _simplex;
      vertices.resize(_dim + 1);
      parent.get_simplex_vertices(idx_below, _dim, dist.size(), vertices.rbegin());

      neighbor_it.resize(_dim + 1);
      neighbor_end.resize(_dim + 1);
      for (dimension_t i = 0; i <= _dim; ++i) {
        auto v = vertices[i];
        neighbor_it[i] = dist.neighbors[v].rbegin();
        neighbor_end[i] = dist.neighbors[v].rend();
      }
    }

    bool has_next(bool all_cofacets = true) {
      for (auto &it0 = neighbor_it[0], &end0 = neighbor_end[0]; it0 != end0; ++it0) {
        neighbor = *it0;
        for (size_t idx = 1; idx < neighbor_it.size(); ++idx) {
          auto &it = neighbor_it[idx], end = neighbor_end[idx];
          while (get_vertex(*it) > get_vertex(neighbor))
            if (++it == end) return false;
          if (get_vertex(*it) != get_vertex(neighbor))
            goto continue_outer;
          else
            neighbor = std::max(neighbor, *it);
        }
        while (k > 0 && vertices[k - 1] > get_vertex(neighbor)) {
          if (!all_cofacets) return false;
          idx_below -= simplex_encoding(vertices[k - 1], k);
          idx_above += simplex_encoding(vertices[k - 1], k + 1);
          --k;
        }
        return true;
continue_outer:;
      }
      return false;
    }

    std::optional<diameter_entry_t> next_raw(bool all_cofacets = true) {
      if (!has_next(all_cofacets)) return std::nullopt;
      ++neighbor_it[0];
      value_t cofacet_diameter = std::max(get_diameter(simplex), get_diameter(neighbor));
      simplex_t cofacet_index = idx_above + simplex_encoding(get_vertex(neighbor), k + 1) + idx_below;
      const coefficient_t modulus = parent.modulus;
      coefficient_t cofacet_coefficient =
        (k & 1 ? modulus - 1 : 1) * parent.get_coefficient(simplex) % modulus;
      return parent.make_diameter_entry(cofacet_diameter, cofacet_index, cofacet_coefficient);
    }
    std::optional<diameter_entry_t> next(bool all_cofacets = true) {
      return next_raw(all_cofacets);
    }
  };

  typedef Simplex_coboundary_enumerator_<DistanceMatrix> Simplex_coboundary_enumerator;

  class Simplex_boundary_enumerator {
    private:
      simplex_t idx_below, idx_above;
      vertex_t j;
      dimension_t k;
      diameter_entry_t simplex;
      dimension_t dim;
      const SimplexEncoding& simplex_encoding;
      const Rips_filtration& parent; // for n, get_max_vertex, compute_diameter

    public:
      Simplex_boundary_enumerator(const Rips_filtration& _parent)
        : simplex_encoding(_parent.simplex_encoding), parent(_parent) {}

      void set_simplex(const diameter_entry_t _simplex, const dimension_t _dim) {
        idx_below = parent.get_index(_simplex);
        idx_above = 0;
        j = parent.n - 1;
        k = _dim;
        simplex = _simplex;
        dim = _dim;
      }

      bool has_next() { return (k >= 0); }

      std::optional<diameter_entry_t> next() {
        if (!has_next()) return std::nullopt;
        j = simplex_encoding.get_max_vertex(idx_below, k + 1, j);

        simplex_t face_index = idx_above - simplex_encoding(j, k + 1) + idx_below;

        // It would make sense to extract the vertices once in set_simplex
        // and pass the proper subset to compute_diameter, but even in cases
        // where this dominates it does not seem to help (probably because we
        // stop at the first coface).
        value_t face_diameter = parent.compute_diameter(face_index, dim - 1);

        const coefficient_t modulus = parent.modulus;
        coefficient_t face_coefficient =
          (k & 1 ? -1 + modulus : 1) * parent.get_coefficient(simplex) % modulus;

        idx_below -= simplex_encoding(j, k + 1);
        idx_above += simplex_encoding(j, k);

        --k;

        return parent.make_diameter_entry(face_diameter, face_index, face_coefficient);
      }
  };
};

#if BOOST_VERSION >= 108100
template <class Key, class T, class H, class E> using Hash_map = boost::unordered_flat_map<Key, T, H, E>;
#else
template <class Key, class T, class H, class E> using Hash_map = boost::unordered_map<Key, T, H, E>;
#endif

template <typename Filtration> class Persistent_cohomology {
  using size_t = typename Filtration::size_t;
  using coefficient_t = typename Filtration::coefficient_t;
  using coefficient_storage_t = typename Filtration::coefficient_storage_t;
  using dimension_t = typename Filtration::dimension_t;
  using value_t = typename Filtration::value_t;
  using vertex_t = typename Filtration::vertex_t;
  using simplex_t = typename Filtration::simplex_t;
  using diameter_simplex_t = typename Filtration::diameter_simplex_t;
  using entry_t = typename Filtration::entry_t;
  using diameter_entry_t = typename Filtration::diameter_entry_t;
  using Entry_hash = typename Filtration::Entry_hash;
  using Equal_index = typename Filtration::Equal_index;
  using Simplex_boundary_enumerator = typename Filtration::Simplex_boundary_enumerator;
  using Simplex_coboundary_enumerator = typename Filtration::Simplex_coboundary_enumerator;
  template<class T>using Greater_diameter_or_smaller_index = typename Filtration::template Greater_diameter_or_smaller_index<T>;

  typedef Compressed_sparse_matrix_<diameter_entry_t> Compressed_sparse_matrix;
  typedef Hash_map<entry_t, size_t, Entry_hash, Equal_index> entry_hash_map;

  Filtration filt;
  const vertex_t n;
  const dimension_t dim_max;
  const coefficient_t modulus;
  const std::vector<coefficient_storage_t> multiplicative_inverse_;
  mutable std::vector<diameter_entry_t> cofacet_entries;
  mutable std::vector<vertex_t> vertices;
  Simplex_boundary_enumerator facets;
  // Creating a new one in each function that needs it wastes a bit of time, but we may need 2 at the same time.
  Simplex_coboundary_enumerator cofacets1, cofacets2;

  coefficient_t multiplicative_inverse(coefficient_t c) const {
    return multiplicative_inverse_[c];
  }
  public:
  Persistent_cohomology(Filtration&& _filt, dimension_t _dim_max, coefficient_t _modulus)
    : filt(std::move(_filt)), n(filt.num_vertices()),
    dim_max(std::min<vertex_t>(_dim_max, n - 2)),
    modulus(_modulus),
    multiplicative_inverse_(multiplicative_inverse_vector<coefficient_storage_t>(_modulus)),
    facets(filt), cofacets1(filt), cofacets2(filt)
  {
    }


  std::optional<diameter_entry_t> get_zero_pivot_facet(const diameter_entry_t simplex, const dimension_t dim) {
    facets.set_simplex(simplex, dim);
    while(true) {
      std::optional<diameter_entry_t> facet = facets.next();
      if (!facet) break;
      if (get_diameter(*facet) == get_diameter(simplex)) return *facet;
    }
    return std::nullopt;
  }

  std::optional<diameter_entry_t> get_zero_pivot_cofacet(const diameter_entry_t simplex, const dimension_t dim) {
    cofacets1.set_simplex(simplex, dim);
    while(true) {
      std::optional<diameter_entry_t> cofacet = cofacets1.next_raw();
      if (!cofacet) break;
      if (get_diameter(*cofacet) == get_diameter(simplex)) return *cofacet;
    }
    return std::nullopt;
  }

  // Apparent pairs are implicit in Ripser.
  // pro: we don't need to store them
  // con: we may have to recompute them many times, and each test is more expensive than emergent pairs
  std::optional<diameter_entry_t> get_zero_apparent_facet(const diameter_entry_t simplex, const dimension_t dim) {
    std::optional<diameter_entry_t> facet = get_zero_pivot_facet(simplex, dim);
    if (!facet) return std::nullopt;
    std::optional<diameter_entry_t> cofacet = get_zero_pivot_cofacet(*facet, dim - 1);
    if (!cofacet || filt.get_index(*cofacet) != filt.get_index(simplex)) return std::nullopt;
    return *facet;
  }

  std::optional<diameter_entry_t> get_zero_apparent_cofacet(const diameter_entry_t simplex, const dimension_t dim) {
    std::optional<diameter_entry_t> cofacet = get_zero_pivot_cofacet(simplex, dim);
    if (!cofacet) return std::nullopt;
    std::optional<diameter_entry_t> facet = get_zero_pivot_facet(*cofacet, dim + 1);
    if (!facet || filt.get_index(*facet) != filt.get_index(simplex)) return std::nullopt;
    return *cofacet;
  }

  bool is_in_zero_apparent_pair(const diameter_entry_t simplex, const dimension_t dim) {
    return get_zero_apparent_cofacet(simplex, dim) || get_zero_apparent_facet(simplex, dim);
  }

  void assemble_columns_to_reduce(std::vector<diameter_simplex_t>& simplices,
      std::vector<diameter_simplex_t>& columns_to_reduce,
      entry_hash_map& pivot_column_index, dimension_t dim) {

#ifdef GUDHI_INDICATE_PROGRESS
    std::cerr << clear_line << "assembling columns" << std::flush;
    std::chrono::steady_clock::time_point next = std::chrono::steady_clock::now() + time_step;
#endif

    columns_to_reduce.clear();
    std::vector<diameter_simplex_t> next_simplices;

    for (diameter_simplex_t& simplex : simplices) {
      cofacets2.set_simplex(filt.make_diameter_entry(simplex, 1), dim - 1);

      while(true) {
        std::optional<diameter_entry_t> cofacet = cofacets2.next(false);
        if (!cofacet) break;
#ifdef GUDHI_INDICATE_PROGRESS
        if (std::chrono::steady_clock::now() > next) {
          std::cerr << clear_line << "assembling " << next_simplices.size()
            << " columns (processing " << std::distance(&simplices[0], &simplex)
            << "/" << simplices.size() << " simplices)" << std::flush;
          next = std::chrono::steady_clock::now() + time_step;
        }
#endif
        if (dim < dim_max) next_simplices.push_back({get_diameter(*cofacet), filt.get_index(*cofacet)});
        // Wouldn't it be cheaper in the reverse order? Seems negligible
        if (!is_in_zero_apparent_pair(*cofacet, dim) &&
            (pivot_column_index.find(get_entry(*cofacet)) == pivot_column_index.end()))
          columns_to_reduce.push_back({get_diameter(*cofacet), filt.get_index(*cofacet)});
      }
    }

    if (dim < dim_max) simplices.swap(next_simplices);

#ifdef GUDHI_INDICATE_PROGRESS
    std::cerr << clear_line << "sorting " << columns_to_reduce.size() << " columns"
      << std::flush;
#endif

    std::sort(columns_to_reduce.begin(), columns_to_reduce.end(),
        Greater_diameter_or_smaller_index<diameter_simplex_t>(filt));
#ifdef GUDHI_INDICATE_PROGRESS
    std::cerr << clear_line << std::flush;
#endif
  }

  template<class OutPair>
    void compute_dim_0_pairs(std::vector<diameter_simplex_t>& edges,
        std::vector<diameter_simplex_t>& columns_to_reduce, OutPair& output_pair) {
      Union_find<vertex_t> dset(n);

      edges = filt.get_edges();
      std::sort(edges.rbegin(), edges.rend(),
          Greater_diameter_or_smaller_index<diameter_simplex_t>(filt));
      std::vector<vertex_t> vertices_of_edge(2);
      for (auto e : edges) {
        // Should we work with pairs of vertices instead of edges, to skip get_simplex_vertices?
        filt.get_simplex_vertices(filt.get_index(e), 1, n, vertices_of_edge.rbegin());
        vertex_t u = dset.find(vertices_of_edge[0]), v = dset.find(vertices_of_edge[1]);

        if (u != v) {
          if (get_diameter(e) != 0)
            output_pair(0, get_diameter(e));
          dset.link(u, v);
        } else if ((dim_max > 0) && !get_zero_apparent_cofacet(filt.make_diameter_entry(e, 1), 1))
          columns_to_reduce.push_back(e);
      }
      if (dim_max > 0) std::reverse(columns_to_reduce.begin(), columns_to_reduce.end());

      for (vertex_t i = 0; i < n; ++i)
        if (dset.find(i) == i) output_pair(0, std::numeric_limits<value_t>::infinity());
    }

  template <typename Column> std::optional<diameter_entry_t> pop_pivot(Column& column) {
    while(!column.empty()) { // At this stage the partial sum is 0
      diameter_entry_t pivot = column.top();
      column.pop();
      while(true) { // At this stage the partial sum is led by pivot
        if (column.empty() || filt.get_index(column.top()) != filt.get_index(pivot)) return pivot;
        coefficient_t sum = (filt.get_coefficient(pivot) + filt.get_coefficient(column.top())) % modulus;
        column.pop();
        if (sum == 0) {
          break;
        }
        filt.set_coefficient(pivot, sum);
      }
    }
    return std::nullopt;
  }

  template <typename Column> std::optional<diameter_entry_t> get_pivot(Column& column) {
    // We could look for the pivot without needing to push it back, but it does not seem to help in benchmarks.
    std::optional<diameter_entry_t> result = pop_pivot(column);
    if (result) column.push(*result);
    return result;
  }

  // Either return the pivot in an emergent pair, or fill working_coboundary with the full coboundary (and return its pivot)
  template <typename Column>
    std::optional<diameter_entry_t> init_coboundary_and_get_pivot(const diameter_entry_t simplex,
        Column& working_coboundary, const dimension_t dim,
        entry_hash_map& pivot_column_index) {
      bool check_for_emergent_pair = true;
      cofacet_entries.clear();
      cofacets2.set_simplex(simplex, dim);
      while(true) {
        std::optional<diameter_entry_t> cofacet = cofacets2.next();
        if (!cofacet) break;
        cofacet_entries.push_back(*cofacet);
        if (check_for_emergent_pair && (get_diameter(simplex) == get_diameter(*cofacet))) {
          if ((pivot_column_index.find(get_entry(*cofacet)) == pivot_column_index.end()) &&
              !get_zero_apparent_facet(*cofacet, dim + 1))
            return *cofacet;
          check_for_emergent_pair = false;
        }
      }
      for (auto cofacet : cofacet_entries) working_coboundary.push(cofacet);
      return get_pivot(working_coboundary);
    }

  // Beyond apparent/emergent pairs, Ripser does the column reduction eagerly.
  // To keep with the lazy paradigm, it would be possible to represent a column as a heap of coboundary_iterator, using the order of the current simplices pointed to by the iterators. The equivalent of `pop` would be ++ on the top iterator, which decreases its priority (or removes it if it reached the end). If we do it right, it could also help us notice if we have twice the same iterator and avoid computing its full coboundary twice (although we may still end up computing the first simplex of the coboundary for both, since it may be the easiest way to detect duplicates without maintaining yet another structure on the side). Another advantage is that its size would be bounded by the number of simplices of dimension d instead of d+1 currently. (Dory seems to do something related but too complicated for me.)
  template <typename Column>
    void add_simplex_coboundary(const diameter_entry_t simplex, const dimension_t dim,
        Column& working_reduction_column, Column& working_coboundary) {
      working_reduction_column.push(simplex);
      cofacets1.set_simplex(simplex, dim);
      while(true) { // stupid C++
        std::optional<diameter_entry_t> cofacet = cofacets1.next();
        if (!cofacet) break;
        working_coboundary.push(*cofacet);
      }
    }

  // add an already reduced column, i.e. add all the simplex coboundaries that were involved in that reduction
  template <typename Column>
    void add_coboundary(Compressed_sparse_matrix& reduction_matrix,
        const std::vector<diameter_simplex_t>& columns_to_reduce,
        const size_t index_column_to_add, const coefficient_t factor,
        const dimension_t dim, Column& working_reduction_column,
        Column& working_coboundary) {
      diameter_entry_t column_to_add = filt.make_diameter_entry(columns_to_reduce[index_column_to_add], factor);
      add_simplex_coboundary(column_to_add, dim, working_reduction_column, working_coboundary);

      for (diameter_entry_t simplex : reduction_matrix.subrange(index_column_to_add)) {
        filt.set_coefficient(simplex, filt.get_coefficient(simplex) * factor % modulus);
        add_simplex_coboundary(simplex, dim, working_reduction_column, working_coboundary);
      }
    }

  template<class OutPair>
    void compute_pairs(const std::vector<diameter_simplex_t>& columns_to_reduce,
        entry_hash_map& pivot_column_index, const dimension_t dim, OutPair& output_pair) {
      Compressed_sparse_matrix reduction_matrix;
      Greater_diameter_or_smaller_index<diameter_entry_t> cmp(filt);
      Heap<diameter_entry_t, std::vector<diameter_entry_t>,
        Greater_diameter_or_smaller_index<diameter_entry_t>>
          working_reduction_column(cmp), working_coboundary(cmp);

#ifdef GUDHI_INDICATE_PROGRESS
      std::chrono::steady_clock::time_point next = std::chrono::steady_clock::now() + time_step;
#endif
      for (size_t index_column_to_reduce = 0; index_column_to_reduce < columns_to_reduce.size();
          ++index_column_to_reduce) {

        diameter_entry_t column_to_reduce = filt.make_diameter_entry(columns_to_reduce[index_column_to_reduce], 1);
        value_t diameter = get_diameter(column_to_reduce);

        reduction_matrix.append_column();

        working_reduction_column.clear(); working_coboundary.clear();

        std::optional<diameter_entry_t> pivot = init_coboundary_and_get_pivot(
            column_to_reduce, working_coboundary, dim, pivot_column_index);
        // When we found an emergent pair, we could avoid checking again below, but it does not seem to gain anything in practice.

        while (true) {
#ifdef GUDHI_INDICATE_PROGRESS
          if (std::chrono::steady_clock::now() > next) {
            std::cerr << clear_line << "reducing column " << index_column_to_reduce + 1
              << "/" << columns_to_reduce.size() << " (diameter " << diameter << ")"
              << std::flush;
            next = std::chrono::steady_clock::now() + time_step;
          }
#endif
          if (pivot) {
            auto pair = pivot_column_index.find(get_entry(*pivot));
            if (pair != pivot_column_index.end()) {
              entry_t other_pivot = pair->first;
              size_t index_column_to_add = pair->second;
              coefficient_t factor =
                modulus - filt.get_coefficient(*pivot) *
                multiplicative_inverse(filt.get_coefficient(other_pivot)) %
                modulus;

              // It saves a little bit (3% on an example, 0% on another) if we pass pivot to add_coboundary and avoid pushing entries smaller than pivot in working_coboundary
              add_coboundary(reduction_matrix, columns_to_reduce, index_column_to_add,
                  factor, dim, working_reduction_column, working_coboundary);

              pivot = get_pivot(working_coboundary);
            } else if (std::optional<diameter_entry_t> e = get_zero_apparent_facet(*pivot, dim + 1); e) {
              filt.set_coefficient(*e, modulus - filt.get_coefficient(*e));

              add_simplex_coboundary(*e, dim, working_reduction_column, working_coboundary);

              pivot = get_pivot(working_coboundary);
            } else {
              value_t death = get_diameter(*pivot);
              output_pair(diameter, death);
              pivot_column_index.insert({get_entry(*pivot), index_column_to_reduce});
              // CubicalRipser suggests caching the column here, at least if it took many operations to reduce it.

              while (true) {
                std::optional<diameter_entry_t> e = pop_pivot(working_reduction_column);
                if (!e) break;
                GUDHI_assert(filt.get_coefficient(*e) > 0);
                reduction_matrix.push_back(*e);
              }
              break;
            }
          } else {
            output_pair(diameter, std::numeric_limits<value_t>::infinity());
            break;
          }
        }
      }
#ifdef GUDHI_INDICATE_PROGRESS
      std::cerr << clear_line << std::flush;
#endif
    }

  // Add a separate output_essential?
  // TODO: Should output_pair also take a simplex_t argument?
  template<class OutDim, class OutPair>
    void compute_barcodes(OutDim&& output_dim, OutPair&& output_pair) {
      std::vector<diameter_simplex_t> simplices, columns_to_reduce;

      output_dim(0);
      compute_dim_0_pairs(simplices, columns_to_reduce, output_pair);

      for (dimension_t dim = 1; dim <= dim_max; ++dim) {
        entry_hash_map pivot_column_index(0, filt, filt);
        pivot_column_index.reserve(columns_to_reduce.size());

        output_dim(dim);
        compute_pairs(columns_to_reduce, pivot_column_index, dim, output_pair);

        if (dim < dim_max)
          assemble_columns_to_reduce(simplices, columns_to_reduce, pivot_column_index,
              dim + 1);
        // If for some odd reason one wanted all the infinite intervals in the last dimension, assemble_columns_to_reduce should give us that.
      }
    }
};

// An example of what Params can be
#if 0
struct Params1 {
  // size_t (not always from Params...) is used for counting (ok) and storage for the index of columns in a Hash_map.
  // To gain on a pair<entry_t,size_t> by reducing size_t, simplex_t has to be smaller than size_t I guess, which is very small, not worth it.
  typedef std::size_t size_t;
  typedef float value_t;
  typedef int8_t dimension_t; // Does it need to be signed? Experimentally no.
  typedef int vertex_t; // Currently needs to be signed for Simplex_coboundary_enumerator_<Compressed_lower_distance_matrix>::has_next. Reducing to int16_t helps perf a bit.
  typedef Gudhi::numbers::uint128_t simplex_t;
  // We could introduce a smaller edge_t, but it is not trivial to separate and probably not worth it
  typedef uint16_t coefficient_storage_t; // For the table of multiplicative inverses
  typedef uint_least32_t coefficient_t; // We need x * y % z to work, but promotion from uint16_t to int is not enough
  static const bool use_coefficients = false;

  // Assumptions used in the code:
  // * dimension_t smaller than vertex_t
  // Check which ones need to be signed, and which ones could be unsigned instead
};
#endif

// Trying to write a magic function
template<class Params, class SimplexEncoding, class DistanceMatrix, class OutDim, class OutPair>
void help2(DistanceMatrix&& dist, int dim_max, typename DistanceMatrix::value_t threshold, unsigned modulus, OutDim&& output_dim, OutPair&& output_pair) {
  typedef Rips_filtration<DistanceMatrix, SimplexEncoding, Params> Filt;
  typedef Persistent_cohomology<Filt> Pcoh;
  Filt filt(std::move(dist), dim_max, threshold, modulus);
  Pcoh(std::move(filt), dim_max, modulus).compute_barcodes(output_dim, output_pair);
}
template<bool use_coefficients_, class simplex_t_, class value_t_> struct TParams {
  // hardcode most options
  typedef std::size_t size_t;
  typedef value_t_ value_t;
  typedef int8_t dimension_t;
  typedef int vertex_t;
  typedef simplex_t_ simplex_t;
  typedef uint16_t coefficient_storage_t;
  typedef uint_least32_t coefficient_t;
  static const bool use_coefficients = use_coefficients_;
};
template<class value_t_> struct TParams2 {
  typedef int vertex_t;
  typedef value_t_ value_t;
};
// Choose simplex encoding
template<bool use_coefficients, class DistanceMatrix, class OutDim, class OutPair>
void help1(DistanceMatrix&& dist, int dim_max, typename DistanceMatrix::value_t threshold, unsigned modulus, OutDim&& output_dim, OutPair&& output_pair) {
  auto n = dist.size();
  // TODO: if n>=3 we could even go to n-3, because the Rips cannot have homology in dimension n-2
  // (e.g. for 3 points, there is no 1-homology)
  if (dim_max > n - 2) dim_max = n - 2; // FIXME: duplicated. problem if n unsigned and 0 or 1.
  int bits_per_vertex = log2up(n);
  int bits_for_coeff = log2up(modulus - 1); // duplicating logic :-( Also, if modulus is something absurd, the diagnostic is too late
  int bitfield_size = bits_per_vertex * (dim_max + 2) + bits_for_coeff;
  if (bitfield_size <= 64) { // use bitfield-64
    typedef TParams<use_coefficients, uint64_t, typename DistanceMatrix::value_t> P;
    help2<P, Bitfield_encoding<P>>(std::move(dist), dim_max, threshold, modulus, output_dim, output_pair);
  } else if (bitfield_size <= 128) { // use bitfield-128
    typedef TParams<use_coefficients, Gudhi::numbers::uint128_t, typename DistanceMatrix::value_t> P;
    help2<P, Bitfield_encoding<P>>(std::move(dist), dim_max, threshold, modulus, output_dim, output_pair);
  } else { // use cns-128
    typedef TParams<use_coefficients, Gudhi::numbers::uint128_t, typename DistanceMatrix::value_t> P;
    help2<P, Cns_encoding<P>>(std::move(dist), dim_max, threshold, modulus, output_dim, output_pair);
  }
  // Does cns-64 have its place?
}
// Select hardcoded Z/2Z or runtime Z/pZ
template<class DistanceMatrix, class OutDim, class OutPair>
void ripser(DistanceMatrix dist, int dim_max, typename DistanceMatrix::value_t threshold, unsigned modulus, OutDim&& output_dim, OutPair&& output_pair) {
  if (modulus == 2)
    help1<false>(std::move(dist), dim_max, threshold, modulus, output_dim, output_pair);
  else
    help1<true >(std::move(dist), dim_max, threshold, modulus, output_dim, output_pair);
}
#if 0
template<class DMParams, class OutDim, class OutPair>
void ripser_auto(Sparse_distance_matrix<DMParams> dist, int dim_max, typename DMParams::value_t threshold, unsigned modulus, OutDim&& output_dim, OutPair&& output_pair) {
  ripser(std::move(dist), dim_max, threshold, modulus, output_dim, output_pair);
}
template<class DMParams, Compressed_matrix_layout Layout, class OutDim, class OutPair>
void ripser_auto(Compressed_distance_matrix<DMParams, Layout> dist, int dim_max, typename DMParams::value_t threshold, unsigned modulus, OutDim&& output_dim, OutPair&& output_pair) {
  typedef typename DMParams::value_t value_t;
  typedef typename DMParams::vertex_t vertex_t;
  if (threshold < std::numeric_limits<value_t>::max()) { // or infinity()
    Sparse_distance_matrix<DMParams> new_dist(dist, threshold);
    ripser(std::move(new_dist), dim_max, threshold, modulus, output_dim, output_pair);
  } else {
    for (vertex_t i = 0; i < dist.size(); ++i) {
      value_t r_i = -std::numeric_limits<value_t>::infinity();
      for (vertex_t j = 0; j < dist.size(); ++j)
        r_i = std::max(r_i, dist(i, j));
      threshold = std::min(threshold, r_i);
      // Should we also compute this when an explicit threshold is passed, in case the user passed an unnecessarily large threshold?
    }
    ripser(std::move(dist), dim_max, threshold, modulus, output_dim, output_pair);
  }
}
// Other = Euclidean. TODO: more robust dispatching... (how?)
template<class DistanceMatrix, class OutDim, class OutPair>
void ripser_auto(DistanceMatrix dist, int dim_max, typename DistanceMatrix::value_t threshold, unsigned modulus, OutDim&& output_dim, OutPair&& output_pair) {
  typedef typename DistanceMatrix::value_t value_t;
  typedef TParams2<value_t> P;
  if (threshold < std::numeric_limits<value_t>::max()) { // or infinity()
    Sparse_distance_matrix<P> new_dist(dist, threshold);
    ripser(std::move(new_dist), dim_max, threshold, modulus, output_dim, output_pair);
  } else {
    Compressed_distance_matrix<P, LOWER_TRIANGULAR> new_dist(dist);
    ripser(std::move(new_dist), dim_max, threshold, modulus, output_dim, output_pair);
  }
}
#endif
template<class DistanceMatrix, class OutDim, class OutPair>
void ripser_auto(DistanceMatrix dist, int dim_max, typename DistanceMatrix::value_t threshold, unsigned modulus, OutDim&& output_dim, OutPair&& output_pair) {
  typedef typename DistanceMatrix::vertex_t vertex_t;
  typedef typename DistanceMatrix::value_t value_t;
  typedef TParams2<value_t> P;
  if constexpr (std::is_same_v<typename DistanceMatrix::Category, Tag_sparse>) {
    ripser(std::move(dist), dim_max, threshold, modulus, output_dim, output_pair);
  } else if (threshold < std::numeric_limits<value_t>::max()) { // or infinity()
    Sparse_distance_matrix<P> new_dist(dist, threshold);
    ripser(std::move(new_dist), dim_max, threshold, modulus, output_dim, output_pair);
  } else if constexpr (std::is_same_v<typename DistanceMatrix::Category, Tag_dense>) {
    for (vertex_t i = 0; i < dist.size(); ++i) {
      value_t r_i = -std::numeric_limits<value_t>::infinity();
      for (vertex_t j = 0; j < dist.size(); ++j)
        r_i = std::max(r_i, dist(i, j));
      threshold = std::min(threshold, r_i);
      // Should we also compute this when an explicit threshold is passed, in case the user passed an unnecessarily large threshold?
    }
    ripser(std::move(dist), dim_max, threshold, modulus, output_dim, output_pair);
  } else {
    Compressed_distance_matrix<P, LOWER_TRIANGULAR> new_dist(dist);
    ripser_auto(std::move(new_dist), dim_max, threshold, modulus, output_dim, output_pair);
  }
}
// - sparse input -> sparse matrix
// - euclidean input & threshold -> sparse matrix (don't build dense matrix)
// - euclidean input & !threshold -> dense matrix
// - dense matrix & threshold -> sparse matrix
// - dense matrix & !threshold -> compute minmax, keep dense
}
#undef GUDHI_assert
#endif // GUDHI_RIPSER_H

/* Relevant benchmarks where different functions dominate the profile
# push
ripser --format point-cloud ripser/examples/o3_1024.txt
# pop
ripser --format point-cloud ripser/examples/o3_4096.txt --threshold 1.6
# apparent_pair (get_simplex_vertices + dist) - edge_collapse would help
ripser --format point-cloud --dim 2 --threshold .7 tore
# apparent_pair (get_simplex_vertices) - no edge collapse possible
ripser --format point-cloud circle24 --dim 25

where equivalent datasets can be obtained through
  tore -> tore3D_1307.off
  circle24:
    t=np.linspace(0, 2*np.pi, 24, endpoint=False)
    np.stack((np.cos(t),np.sin(t))).T
    OR
    gudhi.datasets.generators.points.torus(24,1,'grid')
  o3 (?):
    scipy.stats.ortho_group.rvs(3,4096).reshape(-1,9)
*/
