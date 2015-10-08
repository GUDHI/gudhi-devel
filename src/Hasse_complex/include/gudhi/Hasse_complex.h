/*    This file is part of the Gudhi Library. The Gudhi library 
 *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
 *    library for computational topology.
 *
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014  INRIA Sophia Antipolis-Méditerranée (France)
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef HASSE_COMPLEX_H_
#define HASSE_COMPLEX_H_

#include <boost/iterator/counting_iterator.hpp>

#include <algorithm>

namespace Gudhi {

template < class HasseCpx >
struct Hasse_simplex {
  // Complex_ds must verify that cpx->key(sh) is the order of sh in the filtration

  template< class Complex_ds >
  Hasse_simplex(Complex_ds & cpx
                , typename Complex_ds::Simplex_handle sh)
      : key_(cpx.key(sh))
      , filtration_(cpx.filtration(sh))
      , boundary_() {
    boundary_.reserve(cpx.dimension(sh) + 1);
    for (auto b_sh : cpx.boundary_simplex_range(sh)) {
      boundary_.push_back(cpx.key(b_sh));
    }
  }

  Hasse_simplex(typename HasseCpx::Simplex_key key
                , typename HasseCpx::Filtration_value fil
                , std::vector<typename HasseCpx::Simplex_handle> boundary)
      : key_(key)
      , filtration_(fil)
      , boundary_(boundary) { }

  typename HasseCpx::Simplex_key key_;
  typename HasseCpx::Filtration_value filtration_;
  std::vector<typename HasseCpx::Simplex_handle> boundary_;
};

/** \private 
 * \brief Data structure representing a Hasse diagram, i.e. 
 * a complex where all codimension 1 incidence 
 * relations are explicitly encoded.
 *
 * \implements FilteredComplex
 * \ingroup simplex_tree
 */
template < typename FiltrationValue = double
, typename SimplexKey = int
, typename VertexHandle = int
>
class Hasse_complex {
 public:

  typedef Hasse_simplex<Hasse_complex> Hasse_simp;
  typedef FiltrationValue Filtration_value;
  typedef SimplexKey Simplex_key;
  typedef int Simplex_handle;  //index in vector complex_

  typedef boost::counting_iterator< Simplex_handle > Filtration_simplex_iterator;
  typedef boost::iterator_range<Filtration_simplex_iterator> Filtration_simplex_range;

  typedef typename std::vector< Simplex_handle >::iterator Boundary_simplex_iterator;
  typedef boost::iterator_range<Boundary_simplex_iterator> Boundary_simplex_range;

  typedef typename std::vector< Simplex_handle >::iterator Skeleton_simplex_iterator;
  typedef boost::iterator_range< Skeleton_simplex_iterator > Skeleton_simplex_range;

  /*  only dimension 0 skeleton_simplex_range(...) */
  Skeleton_simplex_range skeleton_simplex_range(int dim = 0) {
    if (dim != 0) {
      std::cerr << "Dimension must be 0 \n";
    }
    return Skeleton_simplex_range(vertices_.begin(), vertices_.end());
  }

  template < class Complex_ds >
  Hasse_complex(Complex_ds & cpx)
      : complex_()
      , vertices_()
      , threshold_(cpx.filtration())
      , num_vertices_()
      , dim_max_(cpx.dimension()) {
    complex_.reserve(cpx.num_simplices());
    int idx = 0;
    for (auto cpx_sh : cpx.filtration_simplex_range()) {
      complex_.push_back(Hasse_simp(cpx, cpx_sh));
      if (dimension(idx) == 0) {
        vertices_.push_back(idx);
      }
      ++idx;
    }
  }

  Hasse_complex()
      : complex_()
      , vertices_()
      , threshold_(0)
      , num_vertices_(0)
      , dim_max_(-1) { }

  size_t num_simplices() {
    return complex_.size();
  }

  Filtration_simplex_range filtration_simplex_range() {
    return Filtration_simplex_range(Filtration_simplex_iterator(0)
                                    , Filtration_simplex_iterator(complex_.size()));
  }

  Simplex_key key(Simplex_handle sh) {
    return complex_[sh].key_;
  }

  Simplex_key null_key() {
    return -1;
  }

  Simplex_handle simplex(Simplex_key key) {
    if (key == null_key()) return null_simplex();
    return key;
  }

  Simplex_handle null_simplex() {
    return -1;
  }

  Filtration_value filtration(Simplex_handle sh) {
    if (sh == null_simplex()) {
      return filtration();
    }
    return complex_[sh].filtration_;
  }

  Filtration_value filtration() {
    return threshold_;
  }

  int dimension(Simplex_handle sh) {
    if (complex_[sh].boundary_.empty()) return 0;
    return complex_[sh].boundary_.size() - 1;
  }

  int dimension() {
    return dim_max_;
  }

  std::pair<Simplex_handle, Simplex_handle> endpoints(Simplex_handle sh) {
    return std::pair<Simplex_handle, Simplex_handle>(complex_[sh].boundary_[0]
                                                     , complex_[sh].boundary_[1]);
  }

  void assign_key(Simplex_handle sh, Simplex_key key) {
    complex_[sh].key_ = key;
  }

  Boundary_simplex_range boundary_simplex_range(Simplex_handle sh) {
    return Boundary_simplex_range(complex_[sh].boundary_.begin()
                                  , complex_[sh].boundary_.end());
  }

  void display_simplex(Simplex_handle sh) {
    std::cout << dimension(sh) << "  ";
    for (auto sh_b : boundary_simplex_range(sh)) std::cout << sh_b << " ";
    std::cout << "  " << filtration(sh) << "         key=" << key(sh);
  }

  void initialize_filtration() {
    Simplex_key key = 0;
    for (auto & h_simp : complex_) {
      h_simp.key_ = key;
      ++key;
    }
  }

  std::vector< Hasse_simp > complex_;
  std::vector<Simplex_handle> vertices_;
  Filtration_value threshold_;
  size_t num_vertices_;
  int dim_max_;
};

template< typename T1, typename T2, typename T3 >
std::istream& operator>>(std::istream & is
                         , Hasse_complex< T1, T2, T3 > & hcpx) {
  assert(hcpx.num_simplices() == 0);

  size_t num_simp;
  is >> num_simp;
  hcpx.complex_.reserve(num_simp);

  std::vector< typename Hasse_complex<T1, T2, T3>::Simplex_key > boundary;
  typename Hasse_complex<T1, T2, T3>::Filtration_value fil;
  typename Hasse_complex<T1, T2, T3>::Filtration_value max_fil = 0;
  int max_dim = -1;
  int key = 0;
  while (read_hasse_simplex(is, boundary, fil))  //read all simplices in the file as a list of vertices
  {
    // insert every simplex in the simplex tree
    hcpx.complex_.push_back(Hasse_simplex< Hasse_complex<T1, T2, T3> >(key, fil, boundary));

    if (max_dim < hcpx.dimension(key)) {
      max_dim = hcpx.dimension(key);
    }
    if (hcpx.dimension(key) == 0) {
      hcpx.vertices_.push_back(key);
    }
    if (max_fil < fil) {
      max_fil = fil;
    }

    ++key;
    boundary.clear();
  }

  hcpx.dim_max_ = max_dim;
  hcpx.threshold_ = max_fil;

  return is;
}

}  // namespace Gudhi

#endif  // HASSE_COMPLEX_H_
