/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Marc Glisse
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <vector>
#include <array>
#include <cmath>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

#include <gudhi/ripser.h>

using namespace Gudhi::ripser;

namespace py = pybind11;
typedef std::vector<std::array< float, 2>> Vf;
typedef std::vector<std::array<double, 2>> Vd;
PYBIND11_MAKE_OPAQUE(Vf);
PYBIND11_MAKE_OPAQUE(Vd);

template<class T>struct Numpy_euclidean {
  typedef Tag_other category;
  typedef int vertex_t;
  typedef T value_t;

  decltype(std::declval<py::array_t<T>&>().template unchecked<2>()) data;

  int size() const { return data.shape(0); }

  T operator()(int i, int j) const {
    T dist = 0;
    for (int k=0; k<data.shape(1); ++k) {
      T diff = data(i, k) - data(j, k);
      dist += diff * diff;
    }
    return std::sqrt(dist);
  }
};

template<class T>struct Full {
  typedef Tag_dense category;
  typedef int vertex_t;
  typedef T value_t;
  decltype(std::declval<py::array_t<T>&>().template unchecked<2>()) data;
  int size() const { return data.shape(0); }
  T operator()(int i, int j) const {
    return data(i, j);
  }
};

template<class vertex_t_, class value_t_>struct DParams {
  typedef vertex_t_ vertex_t;
  typedef value_t_ value_t;
};

template<class DistanceMatrix>
py::list doit(DistanceMatrix&& dist, int max_dimension, typename DistanceMatrix::value_t max_edge_length, unsigned homology_coeff_field) {
  typedef typename DistanceMatrix::value_t T;
  std::vector<std::vector<std::array<T, 2>>> dgms;
  {
    py::gil_scoped_release release;
    auto output = [&](T birth, T death){ dgms.back().push_back({birth, death}); };
    auto switch_dim = [&](int new_dim){
      dgms.emplace_back();
    };
    ripser_auto(std::move(dist), max_dimension, max_edge_length, homology_coeff_field, switch_dim, output);
  }
  py::list ret;
  for (auto&& dgm : dgms)
    ret.append(py::array(py::cast(std::move(dgm))));
  return ret;
}

template<class T>
py::list euclidean(py::array_t<T> points, int max_dimension, T max_edge_length, unsigned homology_coeff_field) {
  Numpy_euclidean<T> dist{points.template unchecked<2>()};
  if(dist.data.ndim() != 2)
    throw std::runtime_error("points must be a 2-dimensional array");

  // ripser_auto should already do that
#if 0
  // optional as a trick to allow destruction where I want it, without hiding dist in some scope that ends too early.
  std::optional<py::gil_scoped_release> release_local(std::in_place);
  compressed_distance_matrix<DParams<int, T>, LOWER_TRIANGULAR> dist(dist_);
  release_local.reset();
#endif

  return doit(std::move(dist), max_dimension, max_edge_length, homology_coeff_field);
}

template<class T>
py::list full(py::array_t<T> matrix, int max_dimension, T max_edge_length, unsigned homology_coeff_field) {
  Full<T> dist{matrix.template unchecked<2>()};
  if(dist.data.ndim() != 2 || dist.data.shape(0) != dist.data.shape(1))
    throw std::runtime_error("Distance matrix must be a square 2-dimensional array");
  return doit(std::move(dist), max_dimension, max_edge_length, homology_coeff_field);
}

py::list lower(py::object low_mat, int max_dimension, double max_edge_length, unsigned homology_coeff_field) {
  std::vector<double> distances;
  int rowi = 0;
  for (auto&& row : low_mat) {
    if (rowi == 0) { ++rowi; continue; }
    int coli = 0;
    for (auto&& elem : row) {
      distances.push_back(elem.cast<double>()); // need a cast?
      if (++coli == rowi) break;
    }
    if (coli < rowi) throw std::invalid_argument("Not enough elements for a lower triangular matrix");
    ++rowi;
  };

  std::optional<py::gil_scoped_release> release_local(std::in_place);
  compressed_distance_matrix<DParams<int, double>, LOWER_TRIANGULAR> dist(std::move(distances));
  release_local.reset();

  return doit(std::move(dist), max_dimension, max_edge_length, homology_coeff_field);
}

template<class V, class T>
py::list sparse(py::array_t<V> is_, py::array_t<V> js_, py::array_t<T> fs_, int num_vertices, int max_dimension, T max_edge_length, unsigned homology_coeff_field) {
  auto is = is_.unchecked();
  auto js = js_.unchecked();
  auto fs = fs_.unchecked();
  if (is.ndim() != 1 || js.ndim() != 1 || fs.ndim() != 1)
    throw std::runtime_error("vertices and filtrations must be 1-dimensional arrays");
  if (is.shape(0) != js.shape(0) || is.shape(0) != js.shape(0))
    throw std::runtime_error("vertices and filtrations must have the same shape");

  typedef DParams<V, T> P;
  typedef sparse_distance_matrix_<P> Dist;
  typedef typename Dist::vertex_diameter_t vertex_diameter_t;

  std::optional<py::gil_scoped_release> release_local(std::in_place);
  std::vector<std::vector<vertex_diameter_t>> neighbors(num_vertices);
  for (py::ssize_t e = 0; e < is.shape(0); ++e) {
    neighbors[is(e)].emplace_back(js(e), fs(e));
    neighbors[js(e)].emplace_back(is(e), fs(e));
  }
  for (size_t i = 0; i < neighbors.size(); ++i)
    std::sort(neighbors[i].begin(), neighbors[i].end());
  Dist dist(std::move(neighbors));
  release_local.reset();

  return doit(std::move(dist), max_dimension, max_edge_length, homology_coeff_field);
}

PYBIND11_MODULE(_ripser, m) {
  py::bind_vector<Vf>(m, "VectorPairFloat" , py::buffer_protocol());
  py::bind_vector<Vd>(m, "VectorPairDouble", py::buffer_protocol());
  // Remove the default for max_dimension?
  m.def("_euclidean", euclidean<float>, py::arg("points").noconvert(), py::arg("max_dimension") = std::numeric_limits<int>::max(), py::arg("max_edge_length") = std::numeric_limits<float>::infinity(), py::arg("homology_coeff_field") = 2);
  m.def("_euclidean", euclidean<double>, py::arg("points"), py::arg("max_dimension") = std::numeric_limits<int>::max(), py::arg("max_edge_length") = std::numeric_limits<double>::infinity(), py::arg("homology_coeff_field") = 2);
  m.def("_full", full<float>, py::arg("matrix").noconvert(), py::arg("max_dimension") = std::numeric_limits<int>::max(), py::arg("max_edge_length") = std::numeric_limits<double>::infinity(), py::arg("homology_coeff_field") = 2);
  m.def("_full", full<double>, py::arg("matrix"), py::arg("max_dimension") = std::numeric_limits<int>::max(), py::arg("max_edge_length") = std::numeric_limits<double>::infinity(), py::arg("homology_coeff_field") = 2);
  m.def("_lower", lower, py::arg("matrix"), py::arg("max_dimension") = std::numeric_limits<int>::max(), py::arg("max_edge_length") = std::numeric_limits<double>::infinity(), py::arg("homology_coeff_field") = 2);
  // We could do a version with long, but copying the arrays of integers shouldn't be too costly
  // TODO doc: duplicate entries forbidden
  m.def("_sparse", sparse<int, float>, py::arg("row"), py::arg("col"), py::arg("data").noconvert(), py::arg("num_vertices"), py::arg("max_dimension") = std::numeric_limits<int>::max(), py::arg("max_edge_length") = std::numeric_limits<double>::infinity(), py::arg("homology_coeff_field") = 2);
  m.def("_sparse", sparse<int, double>, py::arg("row"), py::arg("col"), py::arg("data"), py::arg("num_vertices"), py::arg("max_dimension") = std::numeric_limits<int>::max(), py::arg("max_edge_length") = std::numeric_limits<double>::infinity(), py::arg("homology_coeff_field") = 2);
}

// TODO:
// * input matrice de distances "low" (et aussi "full"? ou que "full" et on convertit côté python?)
// * input matrice de distances sparse "coo"
//
// - sparse input -> sparse matrix
// - euclidean input & threshold -> sparse matrix (don't build dense matrix)
// - euclidean input & !threshold -> dense matrix
// - dense matrix & threshold -> sparse matrix
// - dense matrix & !threshold -> compute minmax, keep dense
