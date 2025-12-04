/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <cstddef>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include <gudhi/vineyard_builder.h>
#include <gudhi/vineyard_helper.h>
#include <python_interfaces/Simplex_tree_interface.h>
#include <python_interfaces/numpy_utils.h>

namespace nb = nanobind;

namespace Gudhi {
namespace vineyard {

class Vineyard_interface
{
  template <typename T>
  struct Range {
    using value_type = T;
    using size_type = std::size_t;
    using const_reference = const value_type&;
    using const_pointer = value_type const*;
    using const_iterator = const_pointer;
    using iterator = const_iterator;
    using difference_type = std::ptrdiff_t;

    Range() {}

    Range(const nb::ndarray<const T, nb::ndim<1>>& range) : range_(range.view()) {}

    size_type size() const { return range_.shape(0); }

    const_reference operator[](size_type n) const { return range_(n); }

   private:
    using Data = decltype(std::declval<nb::ndarray<const T, nb::ndim<1>>>().view());

    Data range_;
  };

 public:
  using value_type = double;
  using Index = Default_vineyard_options::Index;
  using Dimension = Default_vineyard_options::Dimension;

  Vineyard_interface(bool storeRepCycles, Dimension repCyclesDim) : vineyard_(storeRepCycles, repCyclesDim) {}

  template <typename F, typename I>
  void initialize(nb::list boundaryMatrix,
                  nb::ndarray<const I, nb::ndim<1>>& dimensions,
                  nb::ndarray<const F, nb::ndim<1>>& filtrationValues,
                  int numberOfUpdates = 0)
  {
    std::vector<Range<I>> boundaries(boundaryMatrix.size());
    // views are faster than the direct `nb::ndarray`s, but their get destructed together with the array,
    // so everything needs to be store the time the initializing finishs. Normally, they should be no copy,
    // so it should be negligible enough compared to all the rest
    std::vector<nb::ndarray<const I, nb::ndim<1>>> numpy_b(boundaryMatrix.size());

    for (std::size_t i = 0; i < boundaries.size(); ++i) {
      numpy_b[i] = nb::cast<nb::ndarray<const I, nb::ndim<1>>>(boundaryMatrix[i]);
      boundaries[i] = Range<I>(numpy_b[i]);
    }

    vineyard_.initialize(boundaries, Range<I>(dimensions), Range<F>(filtrationValues), numberOfUpdates);
  }

  template <class FilteredComplex>
  void initialize_from_complex(FilteredComplex& complex, int numberOfUpdates = 0)
  {
    std::vector<std::vector<std::uint32_t>> boundaries;
    std::vector<int> dimensions;
    std::vector<double> filtrationValues;

    Gudhi::vineyard::build_boundary_matrix_from_complex(complex, boundaries, dimensions, filtrationValues);
    vineyard_.initialize(boundaries, dimensions, filtrationValues, numberOfUpdates);
  }

  void update(nb::ndarray<const value_type, nb::ndim<1>>& filtrationValues)
  {
    vineyard_.update(Range<value_type>(filtrationValues));
  }

  template <class FilteredComplex>
  void update_from_complex(FilteredComplex& complex)
  {
    std::vector<double> filtrationValues;
    Gudhi::vineyard::build_boundary_matrix_from_complex(complex, filtrationValues);
    vineyard_.update(filtrationValues);
  }

  nb::list get_current_vineyard_view() const
  {
    const auto& vy = vineyard_.get_current_vineyard();
    const std::vector<Index>& nberDims = vineyard_.get_number_of_vines_by_dimension();

    if (vy.empty()) return nb::list();

    std::size_t numberOfSteps = vy[0].size() / nberDims[0];

    nb::list ret;
    for (std::size_t d = 0; d < vy.size(); ++d) {
      // read only, ownership is not transferred
      ret.append(nb::ndarray<nb::numpy, const value_type>(vy[d].data(), {numberOfSteps, nberDims[d], 2}));
    }

    return ret;
  }

  nb::list get_latest_representative_cycles()
  {
    nb::list ret;
    for (auto cycle : vineyard_.get_latest_representative_cycles()) {
      ret.append(_wrap_as_numpy_array(std::move(cycle), cycle.size()));
    }
    return ret;
  }

 private:
  Vineyard_builder<value_type, Default_vineyard_options, true> vineyard_;
};

template <class FilteredComplex>
inline nb::tuple build_python_boundary_matrix_from_complex(FilteredComplex& complex)
{
  std::vector<std::vector<typename FilteredComplex::Simplex_key>> boundaries;
  std::vector<int> dimensions;
  std::vector<typename FilteredComplex::Filtration_value> filtrationValues;

  build_boundary_matrix_from_complex(complex, boundaries, dimensions, filtrationValues);

  nb::list b;
  for (auto& boundary : boundaries) {
    b.append(_wrap_as_numpy_array(std::move(boundary), boundary.size()));
  }
  return nb::make_tuple(std::move(b), _wrap_as_numpy_array(std::move(dimensions), dimensions.size()));
}

}  // namespace vineyard
}  // namespace Gudhi

namespace gvy = Gudhi::vineyard;
using gvyi = gvy::Vineyard_interface;

NB_MODULE(_vineyard_ext, m)
{
  m.attr("__license__") = "MIT";

  nb::class_<gvyi>(m, "Vineyard_interface")
      .def(nb::init<bool,typename gvyi::Dimension>())
      .def("_initialize", &gvyi::initialize<double, int>)
      .def("_initialize", &gvyi::initialize<float, int>)
      // TODO: also for Cubical complex
      .def("_initialize_from_complex", &gvyi::initialize_from_complex<Gudhi::Simplex_tree_interface>)
      .def("_update", &gvyi::update)
      // TODO: also for Cubical complex
      .def("_update_from_complex", &gvyi::update_from_complex<Gudhi::Simplex_tree_interface>)
      .def("_get_current_vineyard_view", &gvyi::get_current_vineyard_view)
      .def("get_latest_representative_cycles", &gvyi::get_latest_representative_cycles);

  // TODO: also for Cubical complex
  m.def("_build_boundary_matrix_from_complex",
        &gvy::build_python_boundary_matrix_from_complex<Gudhi::Simplex_tree_interface>);
}
