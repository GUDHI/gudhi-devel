/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Marc Glisse
 *
 *    Copyright (C) 2024 Inria
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
typedef std::vector<   int> Vi;
typedef std::vector<double> Vd;
typedef std::vector<std::array< float, 2>> V2f;
typedef std::vector<std::array<double, 2>> V2d;
PYBIND11_MAKE_OPAQUE(Vi);
PYBIND11_MAKE_OPAQUE(Vd);
PYBIND11_MAKE_OPAQUE(V2f);
PYBIND11_MAKE_OPAQUE(V2d);

template<class T>struct Full {
  typedef Tag_dense Category;
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
  //static_assert(!std::is_lvalue_reference_v<DistanceMatrix>);
  typedef typename DistanceMatrix::value_t T;
  // We could put everything in a single vector, and return slices of it, but I don't think there is much advantage.
  std::vector<std::vector<std::array<T, 2>>> dgms;
  {
    py::gil_scoped_release release;
    auto output = [&](T birth, T death){
      // Skip empty intervals
      if (birth < death)
        dgms.back().push_back({birth, death});
    };
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
py::list full(py::array_t<T> matrix, int max_dimension, T max_edge_length, unsigned homology_coeff_field) {
  Full<T> dist{matrix.template unchecked<2>()};
  if(dist.data.ndim() != 2 || dist.data.shape(0) != dist.data.shape(1))
    throw std::runtime_error("Distance matrix must be a square 2-dimensional array");
  return doit(std::move(dist), max_dimension, max_edge_length, homology_coeff_field);
}

py::list lower(py::object low_mat, int max_dimension, double max_edge_length, unsigned homology_coeff_field) {
  using Dist = Compressed_distance_matrix<DParams<int, double>, LOWER_TRIANGULAR>;
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

  // optional as a trick to allow destruction where I want it
  std::optional<py::gil_scoped_release> release_local(std::in_place);
  Dist dist(std::move(distances));
  release_local.reset();

  return doit(std::move(dist), max_dimension, max_edge_length, homology_coeff_field);
}

template<class V, class T>
py::list sparse(py::array_t<V> is_, py::array_t<V> js_, py::array_t<T> fs_, int num_vertices, int max_dimension, T max_edge_length, unsigned homology_coeff_field) {
  // Duplicate entries and self loops are forbidden
  auto is = is_.unchecked();
  auto js = js_.unchecked();
  auto fs = fs_.unchecked();
  if (is.ndim() != 1 || js.ndim() != 1 || fs.ndim() != 1)
    throw std::runtime_error("vertices and filtrations must be 1-dimensional arrays");
  if (is.shape(0) != js.shape(0) || is.shape(0) != js.shape(0))
    throw std::runtime_error("vertices and filtrations must have the same shape");

  typedef DParams<V, T> P;
  typedef Sparse_distance_matrix<P> Dist;
  typedef typename Dist::vertex_diameter_t vertex_diameter_t;

  std::optional<py::gil_scoped_release> release_local(std::in_place);
  std::vector<std::vector<vertex_diameter_t>> neighbors(num_vertices);
  for (py::ssize_t e = 0; e < is.shape(0); ++e) {
    neighbors[is(e)].emplace_back(js(e), fs(e));
    neighbors[js(e)].emplace_back(is(e), fs(e));
  }
  // We could easily parallelize this loop, but it is unlikely to be worth it.
  for (size_t i = 0; i < neighbors.size(); ++i)
    std::sort(neighbors[i].begin(), neighbors[i].end());
  Dist dist(std::move(neighbors));
  release_local.reset();

  return doit(std::move(dist), max_dimension, max_edge_length, homology_coeff_field);
}

py::list lower_to_coo(py::object low_mat, double max_edge_length) {
  // Cannot release the GIL since we keep accessing Python objects.
  // TODO: full_to_coo for numpy arrays?
  // Should we compute the cone radius at the same time?
  std::vector<int> is, js;
  std::vector<double> fs;
  int rowi = 0;
  for (auto&& row : low_mat) {
    if (rowi == 0) { ++rowi; continue; }
    int coli = 0;
    for (auto&& elem : row) {
      double d = elem.cast<double>(); // need a cast?
      if (d <= max_edge_length) {
        is.push_back(rowi);
        js.push_back(coli);
        fs.push_back(d);
      }
      if (++coli == rowi) break;
    }
    if (coli < rowi) throw std::invalid_argument("Not enough elements for a lower triangular matrix");
    ++rowi;
  };
  return py::make_tuple(
      py::array(py::cast(std::move(is))),
      py::array(py::cast(std::move(js))),
      py::array(py::cast(std::move(fs))));
}

double lower_cone_radius(py::object low_mat) {
  // It would be more efficient to read the matrix only once
  auto n = py::len(low_mat);
  std::vector<double> maxs(n, -std::numeric_limits<double>::infinity());
  int rowi = 0;
  for (auto&& row : low_mat) {
    if (rowi == 0) { ++rowi; continue; }
    int coli = 0;
    for (auto&& elem : row) {
      double d = elem.cast<double>();
      maxs[rowi] = std::max(maxs[rowi], d);
      maxs[coli] = std::max(maxs[coli], d);
      if (++coli == rowi) break;
    }
    if (coli < rowi) throw std::invalid_argument("Not enough elements for a lower triangular matrix");
    ++rowi;
  };
  return *std::max_element(maxs.begin(), maxs.end());
}

PYBIND11_MODULE(_ripser, m) {
  py::bind_vector<Vi >(m, "VectorInt"       , py::buffer_protocol());
  py::bind_vector<Vd >(m, "VectorDouble"    , py::buffer_protocol());
  py::bind_vector<V2f>(m, "VectorPairFloat" , py::buffer_protocol());
  py::bind_vector<V2d>(m, "VectorPairDouble", py::buffer_protocol());
  // Remove the default for max_dimension?
  m.def("_full", full<float>, py::arg("matrix").noconvert(), py::arg("max_dimension") = std::numeric_limits<int>::max(), py::arg("max_edge_length") = std::numeric_limits<double>::infinity(), py::arg("homology_coeff_field") = 2);
  m.def("_full", full<double>, py::arg("matrix"), py::arg("max_dimension") = std::numeric_limits<int>::max(), py::arg("max_edge_length") = std::numeric_limits<double>::infinity(), py::arg("homology_coeff_field") = 2);
  m.def("_lower", lower, py::arg("matrix"), py::arg("max_dimension") = std::numeric_limits<int>::max(), py::arg("max_edge_length") = std::numeric_limits<double>::infinity(), py::arg("homology_coeff_field") = 2);
  // We could do a version with long, but copying the arrays of integers shouldn't be too costly
  m.def("_sparse", sparse<int, float>, py::arg("row"), py::arg("col"), py::arg("data").noconvert(), py::arg("num_vertices"), py::arg("max_dimension") = std::numeric_limits<int>::max(), py::arg("max_edge_length") = std::numeric_limits<double>::infinity(), py::arg("homology_coeff_field") = 2);
  m.def("_sparse", sparse<int, double>, py::arg("row"), py::arg("col"), py::arg("data"), py::arg("num_vertices"), py::arg("max_dimension") = std::numeric_limits<int>::max(), py::arg("max_edge_length") = std::numeric_limits<double>::infinity(), py::arg("homology_coeff_field") = 2);
  // Not directly an interface to Ripser...
  m.def("_lower_to_coo", lower_to_coo, py::arg("matrix"), py::arg("max_edge_length"));
  m.def("_lower_cone_radius", lower_cone_radius, py::arg("matrix"));
}

// We could also create a RipsComplex class, that allows looking at a simplex, querying its (co)boundary, etc. But I am not convinced it is worth the effort.
