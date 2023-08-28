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
#include <utility>
#include <stdexcept>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <gudhi/Flag_complex_edge_collapser.h>

namespace py = pybind11;

template<class Index, class Filtr>
py::object collapse(py::array_t<Index> is, py::array_t<Index> js, py::array_t<Filtr> fs, int nb_iterations) {
  typedef std::tuple<Index, Index, Filtr> Filtered_edge;
  typedef std::vector<Filtered_edge> Edges;
  py::buffer_info bufi = is.request();
  py::buffer_info bufj = js.request();
  py::buffer_info buff = fs.request();
  if (bufi.ndim != 1 || bufj.ndim != 1 || buff.ndim != 1)
    throw std::runtime_error("Input must be 1-dimensional arrays");
  if (bufi.shape[0] != bufj.shape[0] || bufi.shape[0] != buff.shape[0])
    throw std::runtime_error("Input arrays must have the same size");
  if (buff.shape[0] == 0) {
    py::array_t<Index> indices({{ 2, 0 }}, {{ 0, 0 }});
    py::array_t<Filtr> filtrs;
    return py::make_tuple(std::move(indices), std::move(filtrs));
  }
  auto& edges = *new Edges();
  {
    py::gil_scoped_release release;
    Index n_edges = static_cast<Index>(bufi.shape[0]);
    edges.reserve(n_edges);
    auto strides_i = bufi.strides[0];
    auto strides_j = bufj.strides[0];
    auto strides_f = buff.strides[0];
    for (Index k = 0; k < n_edges; ++k) {
      Index i = *reinterpret_cast<Index*>(static_cast<char*>(bufi.ptr) + k * strides_i);
      Index j = *reinterpret_cast<Index*>(static_cast<char*>(bufj.ptr) + k * strides_j);
      Filtr f = *reinterpret_cast<Filtr*>(static_cast<char*>(buff.ptr) + k * strides_f);
      edges.emplace_back(i, j, f);
    }
    for (int k = 0; k < nb_iterations; ++k) {
      edges = Gudhi::collapse::flag_complex_collapse_edges(std::move(edges), [](auto const&d){return d;});
    }
  }
  py::capsule owner(&edges, [](void*p){ delete reinterpret_cast<Edges*>(p); });
  const auto offset = reinterpret_cast<char*>(&std::get<1>(edges[0])) - reinterpret_cast<char*>(&std::get<0>(edges[0]));
  py::array_t<Index> indices({{ 2, static_cast<py::ssize_t>(edges.size()) }}, {{ offset, sizeof(Filtered_edge) }}, &std::get<0>(edges[0]), owner);
  py::array_t<Filtr> filtrs ({{    static_cast<py::ssize_t>(edges.size()) }}, {{         sizeof(Filtered_edge) }}, &std::get<2>(edges[0]), owner);
  return py::make_tuple(std::move(indices), std::move(filtrs));
}

PYBIND11_MODULE(_edge_collapse, m) {
  m.def("_collapse_edges", collapse<int, float>, py::arg("i").noconvert(), py::arg("j").noconvert(), py::arg("f").noconvert(), py::arg("nb_iterations")=1);
  m.def("_collapse_edges", collapse<py::ssize_t, double>, py::arg("i"), py::arg("j"), py::arg("f"), py::arg("nb_iterations")=1);
}
