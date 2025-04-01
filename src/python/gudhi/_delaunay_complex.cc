/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Thibaud Kloczko
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include "Delaunay_complex_factory.h"

#include <gudhi/Alpha_complex_options.h>
#include <gudhi/MEB_filtration.h>

#include <CGAL/Epeck_d.h>
#include <CGAL/Epick_d.h>

#include <vector>
#include <memory>  // for std::unique_ptr
#include <cstddef>  // for std::size_t

namespace Gudhi {

class Simplex_tree_interface;

namespace delaunay_complex {

// /////////////////////////////////////////////////////////////////////////////
//  Delaunay_complex_interface declaration
// /////////////////////////////////////////////////////////////////////////////

class Delaunay_complex_interface
{
public:
    Delaunay_complex_interface(const std::vector<std::vector<double>>& points,
                               const std::vector<double>& weights,
                               bool fast_version, bool exact_version);

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

// /////////////////////////////////////////////////////////////////////////////
//  Delaunay_complex_interface definition
// /////////////////////////////////////////////////////////////////////////////

Delaunay_complex_interface::Delaunay_complex_interface(const std::vector<std::vector<double>>& points,
                                                       const std::vector<double>& weights,
                                                       bool fast_version, bool exact_version)
{
    const bool weighted = (weights.size() > 0);
    // Specific cases for dimensions 2 and 3
    const std::size_t dimension = ((points.size() > 0) ? points[0].size() : 0);
    if (fast_version) {
        if (weighted) {
            if (dimension == 2) {
                delaunay_ptr_ = std::make_unique<Delaunay_complex_t<CGAL::Epick_d<CGAL::Dimension_tag<2>>, true>>(points, weights, exact_version);
            } else if (dimension == 3) {
                delaunay_ptr_ = std::make_unique<Delaunay_complex_t<CGAL::Epick_d<CGAL::Dimension_tag<3>>, true>>(points, weights, exact_version);
            } else {
                delaunay_ptr_ = std::make_unique<Delaunay_complex_t<CGAL::Epick_d<CGAL::Dynamic_dimension_tag>, true>>(points, weights, exact_version);
            }
        } else {
            if (dimension == 2) {
                delaunay_ptr_ = std::make_unique<Delaunay_complex_t<CGAL::Epick_d<CGAL::Dimension_tag<2>>, false>>(points, exact_version);
            } else if (dimension == 3) {
                delaunay_ptr_ = std::make_unique<Delaunay_complex_t<CGAL::Epick_d<CGAL::Dimension_tag<3>>, false>>(points, exact_version);
            } else {
                delaunay_ptr_ = std::make_unique<Delaunay_complex_t<CGAL::Epick_d<CGAL::Dynamic_dimension_tag>, false>>(points, exact_version);
            }
        }
    } else {
        if (weighted) {
            if (dimension == 2) {
                delaunay_ptr_ = std::make_unique<Delaunay_complex_t<CGAL::Epeck_d<CGAL::Dimension_tag<2>>, true>>(points, weights, exact_version);
            } else if (dimension == 3) {
                delaunay_ptr_ = std::make_unique<Delaunay_complex_t<CGAL::Epeck_d<CGAL::Dimension_tag<3>>, true>>(points, weights, exact_version);
            } else {
                delaunay_ptr_ = std::make_unique<Delaunay_complex_t<CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>, true>>(points, weights, exact_version);
            }
        } else {
            if (dimension == 2) {
                delaunay_ptr_ = std::make_unique<Delaunay_complex_t<CGAL::Epeck_d<CGAL::Dimension_tag<2>>, false>>(points, exact_version);
            } else if (dimension == 3) {
                delaunay_ptr_ = std::make_unique<Delaunay_complex_t<CGAL::Epeck_d<CGAL::Dimension_tag<3>>, false>>(points, exact_version);
            } else {
                delaunay_ptr_ = std::make_unique<Delaunay_complex_t<CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>, false>>(points, exact_version);
            }
        }
    }
}

std::vector<double> Delaunay_complex_interface::get_point(int vh)
{
    return delaunay_ptr_->get_point(vh);
}

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

}

// /////////////////////////////////////////////////////////////////////////////
// Delaunay_complex_interface wrapping
// /////////////////////////////////////////////////////////////////////////////

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;
namespace gdc = Gudhi::delaunay_complex;
using gdci = gdc::Delaunay_complex_interface;

NB_MODULE(_delaunay_complex_ext, m) {
    m.attr("__license__") = "GPL v3";

    nb::enum_<gdc::Delaunay_filtration>(m, "Filtration", "")
            .value("NONE", gdc::Delaunay_filtration::NONE, "Default Delaunay Complex")
            .value("CECH", gdc::Delaunay_filtration::CECH, "Delaunay Cech Complex")
            .value("ALPHA", gdc::Delaunay_filtration::ALPHA, "Alpha Complex");

    nb::class_<gdci>(m, "Delaunay_complex_interface")
            .def(nb::init<const std::vector<std::vector<double>>&, const std::vector<double>&, bool, bool>(), "Constructor")
            .def("create_simplex_tree", &gdci::create_simplex_tree, "")
            .def("get_point", &gdci::get_point, "")
            .def_static("set_float_relative_precision", &gdci::set_float_relative_precision, "")
            .def_static("get_float_relative_precision", &gdci::get_float_relative_precision, "");
}

//
// _delaunay_complex.cc ends here
