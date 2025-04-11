/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - 2024/03 Vincent Rouvreau: Renamed Alpha_complex_factory as Delaunay_complex_factory for DelaunayCechComplex.
 *                                  Factorize create_complex
 *      - 2024/10 Vincent Rouvreau: Add square root filtration values interface
 *      - 2025/03 Thibaud Kloczko: Use nanobind instead of Cython for python bindings
 *      - 2025/04 Hannah Schreiber: Re-add possibility of tensors (numpy, torch etc.) as input
 *      - YYYY/MM Author: Description of the modification
 */

#include <vector>
#include <memory>   // for std::unique_ptr
#include <cstddef>  // for std::size_t

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <CGAL/Epeck_d.h>
#include <CGAL/Epick_d.h>

#include <gudhi/Alpha_complex_options.h>
#include <gudhi/MEB_filtration.h>
#include <gudhi/Alpha_complex.h>
#include <python_interfaces/Simplex_tree_interface.h>
#include <python_interfaces/points_utils.h>

namespace Gudhi {
namespace delaunay_complex {

// /////////////////////////////////////////////////////////////////////////////
//  Delaunay_complex_factory
// /////////////////////////////////////////////////////////////////////////////

enum class Delaunay_filtration : char {
  NONE = 'n',   ///< Delaunay complex
  CECH = 'c',   ///< Delaunay Cech complex
  ALPHA = 'a',  ///< Alpha complex
};

// template Functor that transforms a CGAL point to a vector of double as expected by python bindings
template <typename CgalPointType, bool Weighted>
struct Point_cgal_to_python;

// Specialized Unweighted Functor
template <typename CgalPointType>
struct Point_cgal_to_python<CgalPointType, false> {
  std::vector<double> operator()(CgalPointType const& point) const
  {
    std::vector<double> vd;
    vd.reserve(point.dimension());
    for (auto coord = point.cartesian_begin(); coord != point.cartesian_end(); coord++)
      vd.push_back(CGAL::to_double(*coord));
    return vd;
  }
};

// Specialized Weighted Functor
template <typename CgalPointType>
struct Point_cgal_to_python<CgalPointType, true> {
  std::vector<double> operator()(CgalPointType const& weighted_point) const
  {
    const auto& point = weighted_point.point();
    return Point_cgal_to_python<decltype(point), false>()(point);
  }
};

// Function that transforms a vector of double point to a CGAL point
template <typename CgalPointType>
static CgalPointType pt_python_to_cgal(std::vector<double> const& vec)
{
  return CgalPointType(vec.size(), vec.begin(), vec.end());
}

template <typename Delaunay_complex, typename Kernel, bool Weighted, typename Point_cloud>
bool create_complex(Delaunay_complex& delaunay_complex,
                    Simplex_tree_interface* simplex_tree,
                    const Point_cloud& points,
                    double max_alpha_square,
                    bool exact_version,
                    Delaunay_filtration filtration,
                    bool output_squared_values)
{
  if (filtration == Delaunay_filtration::CECH) {
    if (Weighted) throw std::invalid_argument("Weighted Delaunay-Cech complex is not available");
    // Construct the Delaunay complex
    bool result = delaunay_complex.create_complex(
        *simplex_tree, std::numeric_limits<Simplex_tree_interface::Filtration_value>::infinity(), exact_version, true);
    if (result == true) {
      // Construct the Delaunay-Cech complex by assigning filtration values with MEB
      if (!output_squared_values) {
        Gudhi::cech_complex::assign_MEB_filtration<false>(Kernel(), *simplex_tree, points);
        simplex_tree->prune_above_filtration(std::sqrt(max_alpha_square));
      } else {
        Gudhi::cech_complex::assign_MEB_filtration<true>(Kernel(), *simplex_tree, points);
        simplex_tree->prune_above_filtration(max_alpha_square);
      }
    }
    return result;
  } else {
    if (output_squared_values)
      return delaunay_complex.template create_complex<true>(
          *simplex_tree, max_alpha_square, exact_version, filtration == Delaunay_filtration::NONE);
    else
      return delaunay_complex.template create_complex<false>(
          *simplex_tree, max_alpha_square, exact_version, filtration == Delaunay_filtration::NONE);
  }
}

class Abstract_delaunay_complex
{
 public:
  virtual std::vector<double> get_point(int vh) = 0;

  virtual bool create_simplex_tree(Simplex_tree_interface* simplex_tree,
                                   double max_alpha_square,
                                   Delaunay_filtration filtration,
                                   bool output_squared_values = false) = 0;

  virtual std::size_t num_vertices() const = 0;

  virtual ~Abstract_delaunay_complex() = default;
};

// /////////////////////////////////////////////////////////////////////////////
//  Delaunay_complex_interface declaration
// /////////////////////////////////////////////////////////////////////////////

class Delaunay_complex_interface
{
 public:
  Delaunay_complex_interface(const Sequence2D& points,
                             const Sequence1D& weights,
                             bool fast_version,
                             bool exact_version);

  Delaunay_complex_interface(const Tensor2D& points,
                             const Tensor1D& weights,
                             bool fast_version,
                             bool exact_version);

  Delaunay_complex_interface(const Sequence2D& points, bool fast_version, bool exact_version);

  Delaunay_complex_interface(const Tensor2D& points, bool fast_version, bool exact_version);

  std::vector<double> get_point(int vh);

  void create_simplex_tree(Simplex_tree_interface* simplex_tree,
                           double max_alpha_square,
                           Delaunay_filtration filtration,
                           bool output_squared_values);

  static void set_float_relative_precision(double precision);

  static double get_float_relative_precision();

 private:
  std::unique_ptr<Abstract_delaunay_complex> delaunay_ptr_;
};

template <typename Kernel, bool Weighted = false>
class Delaunay_complex_t final : public Abstract_delaunay_complex
{
 private:
  using Bare_point = typename Kernel::Point_d;
  using Point = std::conditional_t<Weighted, typename Kernel::Weighted_point_d, typename Kernel::Point_d>;
  using Delaunay_complex = Gudhi::alpha_complex::Alpha_complex<Kernel, Weighted>;

 public:
  Delaunay_complex_t(const Sequence2D& points, bool exact_version)
      : exact_version_(exact_version),
        points_(boost::begin(boost::adaptors::transform(points, pt_python_to_cgal<Bare_point>)),
                boost::end(boost::adaptors::transform(points, pt_python_to_cgal<Bare_point>))),
        delaunay_complex_(points_)
  {}

  Delaunay_complex_t(const Sequence2D& points, const Sequence1D& weights, bool exact_version)
      : exact_version_(exact_version),
        points_(boost::begin(boost::adaptors::transform(points, pt_python_to_cgal<Bare_point>)),
                boost::end(boost::adaptors::transform(points, pt_python_to_cgal<Bare_point>))),
        delaunay_complex_(points_, weights)
  {}

  virtual std::vector<double> get_point(int vh) override
  {
    // Can be a Weighted or a Bare point in function of Weighted
    return Point_cgal_to_python<Point, Weighted>()(delaunay_complex_.get_point(vh));
  }

  virtual bool create_simplex_tree(Simplex_tree_interface* simplex_tree,
                                   double max_alpha_square,
                                   Delaunay_filtration filtration,
                                   bool output_squared_values) override
  {
    return create_complex<Delaunay_complex, Kernel, Weighted, std::vector<Bare_point>>(
        delaunay_complex_, simplex_tree, points_, max_alpha_square, exact_version_, filtration, output_squared_values);
  }

  virtual std::size_t num_vertices() const override { return delaunay_complex_.num_vertices(); }

 private:
  bool exact_version_;
  std::vector<Bare_point> points_;
  Delaunay_complex delaunay_complex_;
};

// /////////////////////////////////////////////////////////////////////////////
//  Delaunay_complex_interface definition
// /////////////////////////////////////////////////////////////////////////////

Delaunay_complex_interface::Delaunay_complex_interface(const Sequence2D& points,
                                                       const Sequence1D& weights,
                                                       bool fast_version,
                                                       bool exact_version)
{
  // Specific cases for dimensions 2 and 3
  const std::size_t dimension = ((points.size() > 0) ? points[0].size() : 0);
  if (fast_version) {
    if (dimension == 2) {
      delaunay_ptr_ = std::make_unique<Delaunay_complex_t<CGAL::Epick_d<CGAL::Dimension_tag<2>>, true>>(
          points, weights, exact_version);
    } else if (dimension == 3) {
      delaunay_ptr_ = std::make_unique<Delaunay_complex_t<CGAL::Epick_d<CGAL::Dimension_tag<3>>, true>>(
          points, weights, exact_version);
    } else {
      delaunay_ptr_ = std::make_unique<Delaunay_complex_t<CGAL::Epick_d<CGAL::Dynamic_dimension_tag>, true>>(
          points, weights, exact_version);
    }
  } else {
    if (dimension == 2) {
      delaunay_ptr_ = std::make_unique<Delaunay_complex_t<CGAL::Epeck_d<CGAL::Dimension_tag<2>>, true>>(
          points, weights, exact_version);
    } else if (dimension == 3) {
      delaunay_ptr_ = std::make_unique<Delaunay_complex_t<CGAL::Epeck_d<CGAL::Dimension_tag<3>>, true>>(
          points, weights, exact_version);
    } else {
      delaunay_ptr_ = std::make_unique<Delaunay_complex_t<CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>, true>>(
          points, weights, exact_version);
    }
  }
}

Delaunay_complex_interface::Delaunay_complex_interface(const Tensor2D& points,
                                                       const Tensor1D& weights,
                                                       bool fast_version,
                                                       bool exact_version)
    : Delaunay_complex_interface(_get_sequence_from_tensor(points),
                                 _get_sequence_from_tensor(weights),
                                 fast_version,
                                 exact_version)
{}

Delaunay_complex_interface::Delaunay_complex_interface(const Sequence2D& points, bool fast_version, bool exact_version)
{
  // Specific cases for dimensions 2 and 3
  const std::size_t dimension = ((points.size() > 0) ? points[0].size() : 0);
  if (fast_version) {
    if (dimension == 2) {
      delaunay_ptr_ =
          std::make_unique<Delaunay_complex_t<CGAL::Epick_d<CGAL::Dimension_tag<2>>, false>>(points, exact_version);
    } else if (dimension == 3) {
      delaunay_ptr_ =
          std::make_unique<Delaunay_complex_t<CGAL::Epick_d<CGAL::Dimension_tag<3>>, false>>(points, exact_version);
    } else {
      delaunay_ptr_ = std::make_unique<Delaunay_complex_t<CGAL::Epick_d<CGAL::Dynamic_dimension_tag>, false>>(
          points, exact_version);
    }
  } else {
    if (dimension == 2) {
      delaunay_ptr_ =
          std::make_unique<Delaunay_complex_t<CGAL::Epeck_d<CGAL::Dimension_tag<2>>, false>>(points, exact_version);
    } else if (dimension == 3) {
      delaunay_ptr_ =
          std::make_unique<Delaunay_complex_t<CGAL::Epeck_d<CGAL::Dimension_tag<3>>, false>>(points, exact_version);
    } else {
      delaunay_ptr_ = std::make_unique<Delaunay_complex_t<CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>, false>>(
          points, exact_version);
    }
  }
}

Delaunay_complex_interface::Delaunay_complex_interface(const Tensor2D& points, bool fast_version, bool exact_version)
    : Delaunay_complex_interface(_get_sequence_from_tensor(points), fast_version, exact_version)
{}

std::vector<double> Delaunay_complex_interface::get_point(int vh) { return delaunay_ptr_->get_point(vh); }

void Delaunay_complex_interface::create_simplex_tree(Simplex_tree_interface* simplex_tree,
                                                     double max_alpha_square,
                                                     Delaunay_filtration filtration,
                                                     bool output_squared_values)
{
  // Nothing to be done in case of an empty point set
  if (delaunay_ptr_->num_vertices() > 0) {
    delaunay_ptr_->create_simplex_tree(simplex_tree, max_alpha_square, filtration, output_squared_values);
  }
}

void Delaunay_complex_interface::set_float_relative_precision(double precision)
{
  // cf. CGAL::Epeck_d kernel type in Delaunay_complex_interface
  if (precision <= 0 || precision >= 1) {
    throw std::invalid_argument("Precision must be strictly greater than 0 and lower than 0");
  }
  CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>::FT::set_relative_precision_of_to_double(precision);
}

double Delaunay_complex_interface::get_float_relative_precision()
{
  // cf. CGAL::Epeck_d kernel type in Delaunay_complex_interface
  return CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>::FT::get_relative_precision_of_to_double();
}

}  // namespace delaunay_complex
}  // namespace Gudhi

// /////////////////////////////////////////////////////////////////////////////
// Delaunay_complex_interface wrapping
// /////////////////////////////////////////////////////////////////////////////

namespace nb = nanobind;
namespace gdc = Gudhi::delaunay_complex;
using gdci = gdc::Delaunay_complex_interface;

NB_MODULE(_delaunay_complex_ext, m)
{
  m.attr("__license__") = "GPL v3";

  nb::enum_<gdc::Delaunay_filtration>(m, "Filtration")
      .value("NONE", gdc::Delaunay_filtration::NONE, "Default Delaunay Complex")
      .value("CECH", gdc::Delaunay_filtration::CECH, "Delaunay Cech Complex")
      .value("ALPHA", gdc::Delaunay_filtration::ALPHA, "Alpha Complex");

  nb::class_<gdci>(m, "Delaunay_complex_interface")
      .def(nb::init<const Sequence2D&, const Sequence1D&, bool, bool>(), nb::call_guard<nb::gil_scoped_release>())
      .def(nb::init<const Tensor2D&, const Tensor1D&, bool, bool>())
      .def(nb::init<const Sequence2D&, bool, bool>(), nb::call_guard<nb::gil_scoped_release>())
      .def(nb::init<const Tensor2D&, bool, bool>())
      .def("create_simplex_tree", &gdci::create_simplex_tree)
      .def("get_point", &gdci::get_point, R"pbdoc(
This function returns the point corresponding to a given vertex from the :class:`~gudhi.SimplexTree` (the
same as the k-th input point, where `k=vertex`)

Args:
    vertex: The vertex.
Returns:
    the point.

:raises IndexError: In case the point has no associated vertex in the diagram (because of weights or because it
    is a duplicate).
        )pbdoc")
      .def_static("get_float_relative_precision", &gdci::get_float_relative_precision, R"pbdoc(
Returns:
    The float relative precision of filtration values computation when constructing with
    :code:`precision = 'safe'` (the default).
        )pbdoc")
      .def_static("set_float_relative_precision", &gdci::set_float_relative_precision, R"pbdoc(
Args:
    precision: When constructing :func:`~gudhi.delaunay_cech_complex`, :func:`~gudhi.alpha_complex`, or
        :func:`~gudhi.weighted_alpha_complex` with :code:`precision = 'safe'` (the default), one can
        set the float relative precision of filtration values computed. Default is :code:`1e-5` (cf.
        :func:`~gudhi.DelaunayComplex.get_float_relative_precision`). For more details, please refer to
        `CGAL::Lazy_exact_nt<NT>::set_relative_precision_of_to_double <https://doc.cgal.org/latest/Number_types/classCGAL_1_1Lazy__exact__nt.html>`_

:raises ValueError: If precision is not in (0, 1).
        )pbdoc");
}

//
// _delaunay_complex.cc ends here
