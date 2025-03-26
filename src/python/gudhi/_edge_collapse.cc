/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Marc Glisse
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - 2025/03 Hannah Schreiber: Use nanobind instead of Pybind11 for python bindings.
 *      - YYYY/MM Author: Description of the modification
 */

#include <cstddef>    //std::size_t
#include <stdexcept>  //std::runtime_error
#include <vector>
#include <tuple>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include <gudhi/Flag_complex_edge_collapser.h>

namespace nb = nanobind;

template <class Index, class Filtr>
nb::object collapse(nb::ndarray<Index, nb::ndim<1>> is,
                    nb::ndarray<Index, nb::ndim<1>> js,
                    nb::ndarray<Filtr, nb::ndim<1>> fs,
                    int nb_iterations)
{
  typedef std::tuple<Index, Index, Filtr> Filtered_edge;
  typedef std::vector<Filtered_edge> Edges;

  if (is.ndim() != 1 || js.ndim() != 1 || fs.ndim() != 1)
    throw std::runtime_error("Input must be 1-dimensional arrays");
  if (is.shape(0) != js.shape(0) || is.shape(0) != fs.shape(0))
    throw std::runtime_error("Input arrays must have the same size");

  if (fs.shape(0) == 0) {
    return nb::make_tuple(nb::ndarray<Index, nanobind::numpy>(nullptr, {2, 0}),
                          nb::ndarray<Filtr, nanobind::numpy>(nullptr, {0}));
  }

  Edges edges;
  {
    nb::gil_scoped_release release;
    Index n_edges = static_cast<Index>(is.shape(0));
    edges.reserve(n_edges);
    auto strides_i = is.stride(0);
    auto strides_j = js.stride(0);
    auto strides_f = fs.stride(0);
    for (Index k = 0; k < n_edges; ++k) {
      Index i = *(is.data() + k * strides_i);
      Index j = *(js.data() + k * strides_j);
      Filtr f = *(fs.data() + k * strides_f);
      edges.emplace_back(i, j, f);
    }
    for (int k = 0; k < nb_iterations; ++k) {
      edges = Gudhi::collapse::flag_complex_collapse_edges(std::move(edges), [](auto const& d) { return d; });
    }
  }

  // nanobind needs the strides to be an element count and not a byte size, so every count needs to be of the same size
  // not sure this works with the vector `edges`.
  std::vector<Index>* indices = new std::vector<Index>(edges.size() * 2);
  std::vector<Filtr>* filtrs = new std::vector<Filtr>(edges.size());

  std::size_t i = 0;
  for (const Filtered_edge& e : edges) {
    (*indices)[i] = std::get<0>(e);
    (*indices)[edges.size() + i] = std::get<1>(e);
    (*filtrs)[i] = std::get<2>(e);
    ++i;
  }

  return nb::make_tuple(
      nb::ndarray<Index, nanobind::numpy>(
          indices->data(),
          {2, edges.size()},
          nb::capsule(indices, [](void* p) noexcept { delete reinterpret_cast<std::vector<Index>*>(p); })),
      nb::ndarray<Filtr, nanobind::numpy>(filtrs->data(), {edges.size()}, nb::capsule(filtrs, [](void* p) noexcept {
                                            delete reinterpret_cast<std::vector<Filtr>*>(p);
                                          })));
}

NB_MODULE(_edge_collapse_ext, m)
{
  m.def("_collapse_edges",
        collapse<int, float>,
        nb::arg("i").noconvert(),
        nb::arg("j").noconvert(),
        nb::arg("f").noconvert(),
        nb::arg("nb_iterations") = 1);
  m.def("_collapse_edges",
        collapse<nb::ssize_t, double>,
        nb::arg("i"),
        nb::arg("j"),
        nb::arg("f"),
        nb::arg("nb_iterations") = 1);
}
