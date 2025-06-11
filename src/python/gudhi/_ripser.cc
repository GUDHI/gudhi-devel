/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Marc Glisse
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - 2025/03 Hannah Schreiber: Use nanobind instead of PyBind11 for python bindings.
 *      - YYYY/MM Author: Description of the modification
 */

#include <cassert>
#include <vector>
#include <cmath>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/vector.h>
// #include <nanobind/stl/bind_vector.h>

#include <gudhi/Debug_utils.h>
#include <gudhi/ripser.h>
#include <python_interfaces/numpy_utils.h>

using namespace Gudhi::ripser;

namespace nb = nanobind;

// typedef std::vector<int> Vi;
// typedef std::vector<double> Vd;
// typedef std::vector<std::array<float, 2>> V2f;
// typedef std::vector<std::array<double, 2>> V2d;
// NB_MAKE_OPAQUE(Vi);
// NB_MAKE_OPAQUE(Vd);
// NB_MAKE_OPAQUE(V2f);
// NB_MAKE_OPAQUE(V2d);

template <class T>
struct Full {
  using Data = decltype(std::declval<nb::ndarray<const T, nb::ndim<2>>>().view());
  typedef Tag_dense Category;
  typedef int vertex_t;
  typedef T value_t;
  Data data;

  int size() const { return data.shape(0); }

  T operator()(int i, int j) const { return data(i, j); }
};

template <class vertex_t_, class value_t_>
struct DParams {
  typedef vertex_t_ vertex_t;
  typedef value_t_ value_t;
};

template <class DistanceMatrix>
nb::list doit(DistanceMatrix&& dist,
              int max_dimension,
              typename DistanceMatrix::value_t max_edge_length,
              unsigned homology_coeff_field)
{
  // static_assert(!std::is_lvalue_reference_v<DistanceMatrix>);
  typedef typename DistanceMatrix::value_t T;
  std::vector<std::vector<T>*> dgms;
  {
    nb::gil_scoped_release release;
    auto output = [&](T birth, T death) {
      // Skip empty intervals
      if (birth < death) {
        dgms.back()->push_back(birth);
        dgms.back()->push_back(death);
      }
    };
    auto switch_dim = [&](int new_dim) { dgms.emplace_back(new std::vector<T>()); };
    ripser_auto(std::move(dist), max_dimension, max_edge_length, homology_coeff_field, switch_dim, output);
  }
  nb::list ret;
  for (auto&& dgm : dgms) ret.append(_wrap_as_numpy_array(dgm, dgm->size() / 2, 2));
  return ret;
}

template <class T>
nb::list full(nb::ndarray<const T, nb::ndim<2>> matrix, int max_dimension, T max_edge_length, unsigned homology_coeff_field)
{
  Full<T> dist{matrix.view()};
  if (dist.data.ndim() != 2 || dist.data.shape(0) != dist.data.shape(1))
    throw std::runtime_error("Distance matrix must be a square 2-dimensional array");
  return doit(std::move(dist), max_dimension, max_edge_length, homology_coeff_field);
}

nb::list lower(nb::object low_mat, int max_dimension, double max_edge_length, unsigned homology_coeff_field)
{
  using Dist = Compressed_distance_matrix<DParams<int, double>, LOWER_TRIANGULAR>;
  std::vector<double> distances;
  int rowi = 0;
  for (auto&& row : low_mat) {
    if (rowi == 0) {
      ++rowi;
      continue;
    }
    int coli = 0;
    for (auto&& elem : row) {
      distances.push_back(nb::cast<double>(elem));  // need a cast?
      if (++coli == rowi) break;
    }
    if (coli < rowi) throw std::invalid_argument("Not enough elements for a lower triangular matrix");
    ++rowi;
  };

  // optional as a trick to allow destruction where I want it
  std::optional<nb::gil_scoped_release> release_local(std::in_place);
  Dist dist(std::move(distances));
  release_local.reset();

  return doit(std::move(dist), max_dimension, max_edge_length, homology_coeff_field);
}

template <class V, class T>
nb::list sparse(nb::ndarray<const V, nb::ndim<1>> is,
                nb::ndarray<const V, nb::ndim<1>> js,
                nb::ndarray<const T, nb::ndim<1>> fs,
                int num_vertices,
                int max_dimension,
                T max_edge_length,
                unsigned homology_coeff_field)
{
  auto is_view = is.view();
  auto js_view = js.view();
  auto fs_view = fs.view();

  // Duplicate entries and self loops are forbidden
  if (is_view.ndim() != 1 || js_view.ndim() != 1 || fs_view.ndim() != 1)
    throw std::runtime_error("vertices and filtrations must be 1-dimensional arrays");
  if (is_view.shape(0) != js_view.shape(0) || is_view.shape(0) != js_view.shape(0))
    throw std::runtime_error("vertices and filtrations must have the same shape");

  typedef DParams<V, T> P;
  typedef Sparse_distance_matrix<P> Dist;
  typedef typename Dist::vertex_diameter_t vertex_diameter_t;

  std::optional<nb::gil_scoped_release> release_local(std::in_place);
  std::vector<std::vector<vertex_diameter_t>> neighbors(num_vertices);
  for (std::size_t e = 0; e < is_view.shape(0); ++e) {
    neighbors[is_view(e)].emplace_back(js_view(e), fs_view(e));
    neighbors[js_view(e)].emplace_back(is_view(e), fs_view(e));
  }
  // We could easily parallelize this loop, but it is unlikely to be worth it.
  for (std::size_t i = 0; i < neighbors.size(); ++i) std::sort(neighbors[i].begin(), neighbors[i].end());
  Dist dist(std::move(neighbors));
  release_local.reset();

  return doit(std::move(dist), max_dimension, max_edge_length, homology_coeff_field);
}

nb::tuple lower_to_coo(nb::object low_mat, double max_edge_length)
{
  // Cannot release the GIL since we keep accessing Python objects.
  // TODO: full_to_coo for numpy arrays?
  // Should we compute the cone radius at the same time?
  std::vector<int>*is = new std::vector<int>(), *js = new std::vector<int>();
  std::vector<double>* fs = new std::vector<double>();
  int rowi = 0;
  for (auto&& row : low_mat) {
    if (rowi == 0) {
      ++rowi;
      continue;
    }
    int coli = 0;
    for (auto&& elem : row) {
      double d = nb::cast<double>(elem);  // need a cast?
      if (d <= max_edge_length) {
        is->push_back(rowi);
        js->push_back(coli);
        fs->push_back(d);
      }
      if (++coli == rowi) break;
    }
    if (coli < rowi) throw std::invalid_argument("Not enough elements for a lower triangular matrix");
    ++rowi;
  };

  return nb::make_tuple(
      _wrap_as_numpy_array(is, is->size()), _wrap_as_numpy_array(js, js->size()), _wrap_as_numpy_array(fs, fs->size()));
}

double lower_cone_radius(nb::object low_mat)
{
  // It would be more efficient to read the matrix only once
  auto n = nb::len(low_mat);
  std::vector<double> maxs(n, -std::numeric_limits<double>::infinity());
  int rowi = 0;
  for (auto&& row : low_mat) {
    if (rowi == 0) {
      ++rowi;
      continue;
    }
    int coli = 0;
    for (auto&& elem : row) {
      double d = nb::cast<double>(elem);
      maxs[rowi] = std::max(maxs[rowi], d);
      maxs[coli] = std::max(maxs[coli], d);
      if (++coli == rowi) break;
    }
    if (coli < rowi) throw std::invalid_argument("Not enough elements for a lower triangular matrix");
    ++rowi;
  };
  return *std::min_element(maxs.begin(), maxs.end());
}

NB_MODULE(_ripser_ext, m)
{
  // nb::bind_vector<Vi>(m, "VectorInt");
  // nb::bind_vector<Vd>(m, "VectorDouble");
  // nb::bind_vector<V2f>(m, "VectorPairFloat");
  // nb::bind_vector<V2d>(m, "VectorPairDouble");

  m.attr("__license__") = "MIT";

  // Remove the default for max_dimension?
  m.def("_full",
        full<float>,
        nb::arg("matrix").noconvert(),
        nb::arg("max_dimension") = std::numeric_limits<int>::max(),
        nb::arg("max_edge_length") = std::numeric_limits<double>::infinity(),
        nb::arg("homology_coeff_field") = 2);
  m.def("_full",
        full<double>,
        nb::arg("matrix"),
        nb::arg("max_dimension") = std::numeric_limits<int>::max(),
        nb::arg("max_edge_length") = std::numeric_limits<double>::infinity(),
        nb::arg("homology_coeff_field") = 2);
  m.def("_lower",
        lower,
        nb::arg("matrix"),
        nb::arg("max_dimension") = std::numeric_limits<int>::max(),
        nb::arg("max_edge_length") = std::numeric_limits<double>::infinity(),
        nb::arg("homology_coeff_field") = 2);
  // We could do a version with long, but copying the arrays of integers shouldn't be too costly
  m.def("_sparse",
        sparse<int, float>,
        nb::arg("row"),
        nb::arg("col"),
        nb::arg("data").noconvert(),
        nb::arg("num_vertices"),
        nb::arg("max_dimension") = std::numeric_limits<int>::max(),
        nb::arg("max_edge_length") = std::numeric_limits<double>::infinity(),
        nb::arg("homology_coeff_field") = 2);
  m.def("_sparse",
        sparse<int, double>,
        nb::arg("row"),
        nb::arg("col"),
        nb::arg("data"),
        nb::arg("num_vertices"),
        nb::arg("max_dimension") = std::numeric_limits<int>::max(),
        nb::arg("max_edge_length") = std::numeric_limits<double>::infinity(),
        nb::arg("homology_coeff_field") = 2);
  // Not directly an interface to Ripser...
  m.def("_lower_to_coo", lower_to_coo, nb::arg("matrix"), nb::arg("max_edge_length"));
  m.def("_lower_cone_radius", lower_cone_radius, nb::arg("matrix"));
}

// We could also create a RipsComplex class, that allows looking at a simplex, querying its (co)boundary, etc. But I am
// not convinced it is worth the effort.
